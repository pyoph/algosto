#install.packages("microbenchmark")
#install.packages("matlib")
#install.packages("MASS")
#install.packages("reshape2")
#install.packages("corrplot")
#install.packages("dlpyr")
# library(reshape2)
# library(RobRegression)
# library(Gmedian)
# library(ggplot2)
# library(far)
# library(gridExtra)
# library(microbenchmark)
# library("matlib")
# library("MASS")
# library("corrplot")
# library("dplyr")


#Calcule le vecteur distances selon la méthode choisie

calcule_RMSE_FP_AUC_par_methode <- function(data, methode, Sigma1 = Sigma1, mu1 = mu1, cutoff = qchisq(0.95, df = ncol(data$Z))) {
  distances <- NULL
  Z <- data$Z
  if (methode == "Shrinkage") {
    # Méthode Shrinkage
    med <- covComed(Z)$center
    Sigma <- covCor(Z)
    rmseSigma <- norm(Sigma - Sigma1,"F")
    rmseMed <- sqrt(sum((mu1 - med)^2))
    distances <- calcule_vecteur_distances(Z, med, Sigma)
    outliers <- detectionOutliers(distances,cutoff)
    
    tc <- table(data$labelsVrais[1:(nrow(Z))], as.numeric(outliers)[1:(nrow(Z))])
    tc
    tc <- safe_access_tc(tc)
    if((tc["0","0"] + tc["0","1"]) != 0)
    {faux_positifs   <- round((tc["0", "1"]/(tc["0", "1"] + tc["0", "0"]))*100,2)}
    auc <- round(auc(outliers,data$labelsVrais),2)*100
    
    
    
  } else if (methode == "Comédiane") {
    # Méthode Comédiane
    med <- covComed(Z)$center
    Sigma <- covComed(Z)$cov
    distances <- calcule_vecteur_distances(Z, med, Sigma)
    rmseSigma <- norm(Sigma - Sigma1,"F")
    rmseMed <- sqrt(sum((mu1 - med)^2))
    distances <- calcule_vecteur_distances(Z, med, Sigma)
    outliers <- detectionOutliers(distances,cutoff)
    
    tc <- table(data$labelsVrais[1:(nrow(Z))], as.numeric(outliers)[1:(nrow(Z))])
    tc
    tc <- safe_access_tc(tc)
    if((tc["0","0"] + tc["0","1"]) != 0)
    {faux_positifs   <- round((tc["0", "1"]/(tc["0", "1"] + tc["0", "0"]))*100,2)}
    auc <- round(auc(outliers,data$labelsVrais),2)*100
    
  } else if (methode == "OGK") {
    # Méthode OGK
    OGK_result <- covOGK(Z,sigmamu = s_mad)
    med <- OGK_result$center
    Sigma <- OGK_result$cov
    distances <- calcule_vecteur_distances(Z, med, Sigma)
    outliers <- detectionOutliers(distances,cutoff)
    
    tc <- table(data$labelsVrais[1:(nrow(Z))], as.numeric(outliers)[1:(nrow(Z))])
    tc
    tc <- safe_access_tc(tc)
    if((tc["0","0"] + tc["0","1"]) != 0)
    {faux_positifs   <- round((tc["0", "1"]/(tc["0", "1"] + tc["0", "0"]))*100,2)}
    auc <- round(auc(outliers,data$labelsVrais),2)*100
    
    
  } else if (methode == "Cov Empirique") {
    
    med <- colMeans(Z)
    Sigma <- cov(Z)
    distances <- calcule_vecteur_distances(Z, med, Sigma)
    outliers <- detectionOutliers(distances,cutoff)
    
    tc <- table(data$labelsVrais[1:(nrow(Z))], as.numeric(outliers)[1:(nrow(Z))])
    tc
    tc <- safe_access_tc(tc)
    if((tc["0","0"] + tc["0","1"]) != 0)
    {faux_positifs   <- round((tc["0", "1"]/(tc["0", "1"] + tc["0", "0"]))*100,2)}
    auc <- round(auc(outliers,data$labelsVrais),2)*100
    
  } else if (methode == "Offline") {
    # Méthode Offline
    #med <- RobVar(Z)$median
    resultats <- estimation(Z,methodeEstimation = "offline")
    med <- resultats$med
    Sigma <- resultats$SigmaOffline
    distances <- calcule_vecteur_distances(Z, med, Sigma)
    outliers <- detectionOutliers(distances,cutoff)
    
    tc <- table(data$labelsVrais[1:(nrow(Z))], as.numeric(outliers)[1:(nrow(Z))])
    tc
    tc <- safe_access_tc(tc)
    if((tc["0","0"] + tc["0","1"]) != 0)
    {faux_positifs   <- round((tc["0", "1"]/(tc["0", "1"] + tc["0", "0"]))*100,2)}
    auc <- round(auc(outliers,data$labelsVrais),2)*100
    
    
  } else if (methode == "Online") {
    # Méthode Online
    resultats <- estimation(Z,methodeEstimation = "online")
    med <- resultats$med
    Sigma <- resultats$SigmaOnline
    #miter <- results$miter
    #U <- results$U
    #lambda <- results$lambdaIter
    distances <- resultats$distances
    outliers <- detectionOutliers(distances,cutoff)
    
    tc <- table(data$labelsVrais[1:(nrow(Z))], as.numeric(outliers)[1:(nrow(Z))])
    tc
    tc <- safe_access_tc(tc)
    if((tc["0","0"] + tc["0","1"]) != 0)
    {faux_positifs   <- round((tc["0", "1"]/(tc["0", "1"] + tc["0", "0"]))*100,2)}
    auc <- round(auc(outliers,data$labelsVrais),2)*100
    
  } 
  else if (methode == "streaming") {
    # Méthode streaming
    resultats <- estimation(Z,methodeEstimation = "streaming")
    med <- resultats$med
    Sigma <- resultats$SigmaStreaming
    #miter <- results$miter
    #U <- results$U
    #lambda <- results$lambdaIter
    distances <- resultats$distances
    outliers <- detectionOutliers(distances,cutoff)
    
    tc <- table(data$labelsVrais[1:(nrow(Z))], as.numeric(outliers)[1:(nrow(Z))])
    tc
    tc <- safe_access_tc(tc)
    if((tc["0","0"] + tc["0","1"]) != 0)
    {faux_positifs   <- round((tc["0", "1"]/(tc["0", "1"] + tc["0", "0"]))*100,2)}
    auc <- round(auc(outliers,data$labelsVrais),2)*100
    
  } 
  return(distances)
}


#Calcul des distances de Mahalanobis à partir de la médiane géométrique et de la covariance estimées

calcule_vecteur_distances <- function(Z,m,Sigma)

{

  distances <- rep(0,n)

  for  (i in 1:nrow(Z))
  {
    distances[i] = as.numeric(Z[i,] - m) %*% solve(Sigma) %*% (as.numeric(t(Z[i,] - m)))

  }

  return (distances)
}

# Fonction pour calculer le cutoff (deuxième seuil)

calcule_cutoff <- function(distances, d) {

  # Constantes et quantiles
  quantile_95 <- qchisq(0.95, df = d)
  quantile_50 <- qchisq(0.5, df = d)
  median_d <- median(distances)

  # Calcul du cutoff
  cv <- (quantile_95 * median_d) / quantile_50

  return(cv)
}


#Détection des outliers à partir d'un vecteur de distances de Mahalanobis et d'un cutoff


detectionOutliers <- function(distances,cutoff)
  
{
  
  outliers_labels <- rep(0,length(distances))

  
 
  
  for (i in 1:length(distances))
  {
    

    if (distances[i] > cutoff) {outliers_labels[i] <- 1}
    else {outliers_labels[i] <- 0}
  }
  
  return (outliers_labels)
  
}


