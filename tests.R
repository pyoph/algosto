r = 35

contamin = "moyenne_variance"
cluster = FALSE

sigmaSq0 <- (1:d); sigmaSq0 <- sigmaSq0 / mean(sigmaSq0)
SigmaContamin <- diag(sqrt(sigmaSq0)) %*% toeplitz(0.8^(0:(d-1))) %*% diag(sqrt(sigmaSq0))

data <- genererEchantillon(n,d,mu1,mu2 = 20*rep((-1)^(1:d)/sqrt(d), d),p1 = 1- r/100,r/100,Sigma1,Sigma2 = 0.1*SigmaContamin,contamin,cluster)

Z = data$Z

resultats = StreamingOutlierDetection(Z,batch = 1)

Sigma = resultats$Sigma
n_obs = d
evalues_matrix <- matrix(0, nrow = nrow(Z), ncol = d)
evalues1_matrix <- matrix(0, nrow = nrow(Z), ncol = d)

for(t in (1:nrow(Z))){
evalues_matrix[t,] <- eigen(Sigma[t,,])$values
evalues1_matrix[t,] <- eigen(Sigma1)$values


}

#Exclusion des valeurs propres trop grandes
reduce_dimension = function(Sigma){

#Exclusion des valeurs aberrantes pour l'estimation des valeurs propres 

distances = rep(0, nrow(Z))
outliers_labels = rep(0,nrow(Z))
cutoffCorr = rep(0,nrow(Z))

for (i in (1:nrow(Z))){
  
  
  lambda = eigen(Sigma[i,,])$values
  Q <- quantile(lambda, probs = c(0.1, 0.9))
  IQR <- Q[2] - Q[1]
  limites <- c(Q[1] - 1.5 * IQR, Q[2] + 1.5 * IQR)
  indices_non_aberrants <- which(lambda >= max(limites[1],0) & lambda <= limites[2])
  lambda_filtre <- lambda[lambda >= max(limites[1],0) & lambda <= limites[2]]
  #print(paste0("indices non aberrants ", indices_non_aberrants))
  #print(paste0("lambda_filtre ",lambda_filtre))
  distances[i] = mahalanobis_generalizedRcpp(Z[i,indices_non_aberrants],resultats$miter[i,indices_non_aberrants],eigen(Sigma[i,,])$vectors[indices_non_aberrants,indices_non_aberrants], eigen(Sigma[i,,])$values[indices_non_aberrants])
  S = distances[i]
  
         
 
    cutoffCorr[i]  = qchisq(.95,df = d)*median(resultats$distances[1:i])/qchisq(.5,df = d)
    if (distances[i] > cutoffCorr[i]) {outliers_labels[i] = 1}
}
return(outliers_labels)
}

outlier_label = reduce_dimension(resultats$Sigma)

table(data$labelsVrais,outlier_label)


table(data$labelsVrais,resultats$outlier_labels)

# Configuration de la zone de graphique
par(mfrow = c(2, 5))

# Générer les 10 graphiques
for (i in 1:d) {
  # Calculer les différences pour la i-ème composante
  differences <- evalues_matrix[,i] - evalues1_matrix[,i]
  
  # Créer le graphique
  plot(1:nrow(Z), differences, 
       type = "b", 
       main = paste("Composante", i),
       xlab = "Temps/Itération",
       ylab = paste("Différence λ", i, "- λ", i, "'"))
  
  # Ajouter une ligne horizontale à zéro pour référence
  abline(h = 0, col = "red", lty = 2)
}

#eigen(Sigma[nrow(Z),,])$values

table(data$labelsVrais,resultats$outlier_labels)

table(data$labelsVrais,outliers_labels)

lambda= eigen(Sigma)$values
lambda
Q <- quantile(lambda, probs = c(0.1, 0.9))
IQR <- Q[2] - Q[1]

limites <- c(Q[1] - 1.5 * IQR, Q[2] + 1.5 * IQR)

lambda_filtre <- lambda[lambda >= limites[1] & lambda <= limites[2]]
indices_non_aberrants <- which(lambda >= max(limites[1],0) & lambda <= max(limites[2],0))
