# KL divergence btw F0 and F1

rm(list=ls())
library(nlme); library(flexmix); library(fields)

# Functions
KL <- function(parms1, parms2){
  invSigma2 <- solve(parms2$Sigma)
  0.5*(log(det(parms2$Sigma)/det(parms1$Sigma)) - d + sum(diag(invSigma2%*%parms1$Sigma)) +
         t(parms2$mu-parms1$mu)%*%invSigma2%*%(parms2$mu-parms1$mu))[1, 1]
}

FrobeniusNormError = function(Sigmahat,Sigma1){
  return (norm(Sigmahat - Sigma1,"F"))
}


taux_detection <- function(vrais_labels, labels_predits) {
  # S'assurer que les entrées sont binaires
  if (!all(vrais_labels %in% c(0, 1)) || !all(labels_predits %in% c(0, 1))) {
    stop("Les vecteurs doivent contenir uniquement 0 (normal) ou 1 (outlier).")
  }
  
  # Calcul des éléments de la matrice de confusion
  VP <- sum(vrais_labels == 1 & labels_predits == 1)  # Vrai Positifs
  FP <- sum(vrais_labels == 0 & labels_predits == 1)  # Faux Positifs
  FN <- sum(vrais_labels == 1 & labels_predits == 0)  # Faux Négatifs
  VN <- sum(vrais_labels == 0 & labels_predits == 0)  # Vrai Négatifs
  
  # Taux
  TPR <- if ((VP + FN) > 0) VP / (VP + FN) else NA  # Sensibilité
  FPR <- if ((FP + VN) > 0) FP / (FP + VN) else NA  # 1 - Spécificité
  
  return(list(
    TPR = TPR,
    FPR = FPR,
    VP = VP,
    FP = FP,
    FN = FN,
    VN = VN
  ))
}




# Contamination parms: F1
ParmsF1 <- function(m1, k1, l1, rho1){
  d <- length(m1)
  mu1 <- k1*m1
  sigmaSq1 <- l1*sigmaSq0
  Sigma1 <- diag(sqrt(sigmaSq1)) %*% toeplitz(rho1^(0:(d-1))) %*% diag(sqrt(sigmaSq1))
  return(list(mu1=mu1, Sigma1=Sigma1))
}


# Dims
d <- 10

# Null parms: F0
mu0 <- rep(0, d)
sigmaSq0 <- (1:d); sigmaSq0 <- sigmaSq0 / mean(sigmaSq0)
rho0 <- 0.3
Sigma0 <- diag(sqrt(sigmaSq0)) %*% toeplitz(rho0^(0:(d-1))) %*% diag(sqrt(sigmaSq0))
parms0 <- list(mu=mu0, Sigma=Sigma0)

# # Null (null) parms: F0
# mu0 <- rep(0, d)
# sigmaSq0 <- rep(1, d)
# rho0 <- 0
# Sigma0 <- diag(sqrt(sigmaSq0)) %*% toeplitz(rho0^(0:(d-1))) %*% diag(sqrt(sigmaSq0))
# parms0 <- list(mu=mu0, Sigma=Sigma0)

# Contamination
m1 <- rep(1/sqrt(d), d)
KL(parms1=parms0, parms2=ParmsF1(m1, 1, 1, 0.5))

par(mfrow=c(3, 2))
k1grid <- 2^(-4:5); 
plot(k1grid, sapply(k1grid, function(k1){KL(parms0, ParmsF1(m1, k1, 1, rho0))}), log='xy')
l1grid <- 2^(-5:5); 
plot(l1grid, sapply(l1grid, function(l1){KL(parms0, ParmsF1(m1, 0, l1, rho0))}), log='xy')
rho1grid <- seq(0.005, 0.995, by=0.005); 
plot(rho1grid, sapply(rho1grid, function(rho1){KL(parms0, ParmsF1(m1, 0, 1, rho1))}), log='y')

KLk1l1 <- sapply(k1grid, function(k1){sapply(l1grid, function(l1){KL(parms0, ParmsF1(m1, k1, l1, rho0))})})
image.plot(log10(l1grid), log10(k1grid), log10(KLk1l1))
KLk1rho1 <- sapply(k1grid, function(k1){sapply(rho1grid, function(rho1){KL(parms0, ParmsF1(m1, k1, 1, rho1))})})
image.plot(rho1grid, log10(k1grid), log10(KLk1rho1))
KLrho1l1 <- sapply(rho1grid, function(rho1){sapply(l1grid, function(l1){KL(parms0, ParmsF1(m1, 0, l1, rho1))})})
KLrho1l1[KLrho1l1==0] <- NA
image.plot(log10(l1grid), rho1grid, log10(KLrho1l1))


##########################
#Tests effets sur norme de Frobenius variation k
#############################



erreurNormeFrobenius = array(0, dim = c(length(contamin_rate), length(methodes),length(k1grid)))
faux_positifs = array(0, dim = c(length(contamin_rate), length(methodes),length(k1grid)))
faux_negatifs = array(0, dim = c(length(contamin_rate), length(methodes),length(k1grid)))

for (i in seq_along(k1grid))
{
for (j in seq_along(contamin_rate))
{
  
  r = contamin_rate[j]
  print(paste("r ", r))
  k = k1grid[i]
  print(paste("k ",k))
  
  data <- genererEchantillon(n,d,mu1,mu2 = k*m1,p1 = 1- r/100,r/100,Sigma1,Sigma2 = Sigma1,contamin = "moyenne",cluster = FALSE)
  Z = data$Z
  compt = 1
for (m in methodes)  
{
  
  if(m == "online"){
    print (m)
  resultats = StreamingOutlierDetectionCorrectedThreshold(Z,batch = 1)
  Sigma = resultats$Sigma[nrow(Z),,]
  outliers = resultats$outlier_labels
  }
  if(m == "streaming"){
    print (m)
    resultats = StreamingOutlierDetectionCorrectedThreshold(Z,batch = ncol(Z))
    Sigma = resultats$Sigma[nrow(Z),,]
    outliers = resultats$outlier_labels
  }
  
  
  print(paste("i = ",i))
  print(paste("compt = ",compt))
  print(paste("j = ",j))
  #print("m = ",m)
  erreurNormeFrobenius[j,compt,i] = norm(Sigma - Sigma1,"F")  
  
  
  
  tc <- table(data$labelsVrais[1:(nrow(Z))], as.numeric(outliers)[1:(nrow(Z))])
  if ("0" %in% rownames(tc)){
  if((tc["0","0"] + tc["0","1"]) != 0)
  {faux_positifs[j,compt,i]   <-(tc["0", "1"]/(tc["0", "1"] + tc["0", "0"]))*100}
  }
  if ("1" %in% rownames(tc)){
  if((tc["1","0"] + tc["1","1"]) != 0)
  {faux_negatifs[j,compt,i]   <-(tc["1", "1"]/(tc["1", "1"] + tc["1", "0"]))*100}}
  compt = compt + 1
} 
  
}

}

save(erreurNormeFrobenius,file ="erreurNormeFrobeniusKGrid.Rdata") 
save(faux_positifs,file ="faux_positifsKGrid.Rdata") 
save(faux_negatifs,file ="faux_negatifsKGrid.Rdata") 

plot(k1grid,faux_negatifs[9,2,])

###########################################
#Variation de l1grid
###########################################


erreurNormeFrobeniusL1 = array(0, dim = c(length(contamin_rate), length(methodes),length(l1grid)))
faux_positifsL1 = array(0, dim = c(length(contamin_rate), length(methodes),length(l1grid)))
faux_negatifsL1 = array(0, dim = c(length(contamin_rate), length(methodes),length(l1grid)))




for (i in seq_along(l1grid))
{
  for (j in seq_along(contamin_rate))
  {
    
    r = contamin_rate[j]
    print(paste("r ", r))
    k = l1grid[i]
    print(paste("k ",k))
    
    data <- genererEchantillon(n,d,mu1,mu2 = mu1,p1 = 1- r/100,r/100,Sigma1,Sigma2 = k*Sigma1,contamin = "variance",cluster = FALSE)
    Z = data$Z
    compt = 1
    for (m in methodes)  
    {
      
      if(m == "online"){
        print (m)
        resultats = StreamingOutlierDetectionCorrectedThreshold(Z,batch = 1)
        Sigma = resultats$Sigma[nrow(Z),,]
        outliers = resultats$outlier_labels
      }
      if(m == "streaming"){
        print (m)
        resultats = StreamingOutlierDetectionCorrectedThreshold(Z,batch = ncol(Z))
        Sigma = resultats$Sigma[nrow(Z),,]
        outliers = resultats$outlier_labels
      }
      
      
      print(paste("i = ",i))
      print(paste("compt = ",compt))
      print(paste("j = ",j))
      #print("m = ",m)
      erreurNormeFrobeniusL1[j,compt,i] = norm(Sigma - Sigma1,"F")  
      
      
      
      tc <- table(data$labelsVrais[1:(nrow(Z))], as.numeric(outliers)[1:(nrow(Z))])
      if ("0" %in% rownames(tc)){
        if((tc["0","0"] + tc["0","1"]) != 0)
        {faux_positifsL1[j,compt,i]   <-(tc["0", "1"]/(tc["0", "1"] + tc["0", "0"]))*100}
      }
      if ("1" %in% rownames(tc)){
        if((tc["1","0"] + tc["1","1"]) != 0)
        {faux_negatifsL1[j,compt,i]   <-(tc["1", "1"]/(tc["1", "1"] + tc["1", "0"]))*100}}
      compt = compt + 1
    } 
    
  }
  
}



###########################################
#Variation de rho1grid
###########################################


erreurNormeFrobeniusRho = array(0, dim = c(length(contamin_rate), length(methodes),length(rho1grid)))
faux_positifsRho = array(0, dim = c(length(contamin_rate), length(methodes),length(rho1grid)))
faux_negatifsRho = array(0, dim = c(length(contamin_rate), length(methodes),length(rho1grid)))



for (i in seq_along(rho1grid))
{
  for (j in seq_along(contamin_rate))
  {
    
    r = contamin_rate[j]
    print(paste("r ", r))
    k = rho1grid[i]
    print(paste("k ",k))
    Sigma1 <- diag(sqrt(sigmaSq0)) %*% toeplitz(k^(0:(d-1))) %*% diag(sqrt(sigmaSq0))
    data <- genererEchantillon(n,d,mu1,mu2 = mu1,p1 = 1- r/100,r/100,Sigma1,Sigma2 = Sigma1,contamin = "moyenne",cluster = FALSE)
    Z = data$Z
    compt = 1
    for (m in methodes)  
    {
      
      if(m == "online"){
        print (m)
        resultats = StreamingOutlierDetectionCorrectedThreshold(Z,batch = 1)
        Sigma = resultats$Sigma[nrow(Z),,]
        outliers = resultats$outlier_labels
      }
      if(m == "streaming"){
        print (m)
        resultats = StreamingOutlierDetectionCorrectedThreshold(Z,batch = ncol(Z))
        Sigma = resultats$Sigma[nrow(Z),,]
        outliers = resultats$outlier_labels
      }
      
      
      print(paste("i = ",i))
      print(paste("compt = ",compt))
      print(paste("j = ",j))
      #print("m = ",m)
      erreurNormeFrobeniusRho[j,compt,i] = norm(Sigma - Sigma1,"F")  
      
      
      
      tc <- table(data$labelsVrais[1:(nrow(Z))], as.numeric(outliers)[1:(nrow(Z))])
      if ("0" %in% rownames(tc)){
        if((tc["0","0"] + tc["0","1"]) != 0)
        {faux_positifsRho[j,compt,i]   <-(tc["0", "1"]/(tc["0", "1"] + tc["0", "0"]))*100}
      }
      if ("1" %in% rownames(tc)){
        if((tc["1","0"] + tc["1","1"]) != 0)
        {faux_negatifsRho[j,compt,i]   <-(tc["1", "1"]/(tc["1", "1"] + tc["1", "0"]))*100}}
      compt = compt + 1
    } 
    
  }
  
}

#################################################
#Variation de k1, de l1 et de rho
#################################################

erreurNormeFrobeniusk1l1rho1 = array(0, dim = c(length(contamin_rate), length(methodes),length(k1grid),length(l1grid),length(rho1grid)))
faux_positifsk1l1rho1 = array(0, dim = c(length(contamin_rate), length(methodes),length(k1grid),length(l1grid),length(rho1grid)))
faux_negatifsk1l1rho1 = array(0, dim = c(length(contamin_rate), length(methodes),length(k1grid),length(l1grid),length(rho1grid)))


for (t in seq_along(rho1grid)){
for (i in seq_along(k1grid))
{
  for (s in seq_along(l1grid)){
    
  for (j in seq_along(contamin_rate))
  {
    
    r = contamin_rate[j]
    print(paste("r ", r))
    k = k1grid[i]
    g = l1grid[s]
    w = rho1grid[t]
    print(paste("k ",k))
    Sigma1 <- diag(sqrt(sigmaSq0)) %*% toeplitz(w^(0:(d-1))) %*% diag(sqrt(sigmaSq0))
    data <- genererEchantillon(n,d,mu1,mu2 = m1*k,p1 = 1- r/100,r/100,Sigma1,Sigma2 = g*Sigma1,contamin = "moyenne_variance",cluster = FALSE)
    Z = data$Z
    compt = 1
    for (m in methodes)  
    {
      
      if(m == "online"){
        print (m)
        resultats = StreamingOutlierDetectionCorrectedThreshold(Z,batch = 1)
        Sigma = resultats$Sigma[nrow(Z),,]
        outliers = resultats$outlier_labels
      }
      if(m == "streaming"){
        print (m)
        resultats = StreamingOutlierDetectionCorrectedThreshold(Z,batch = ncol(Z))
        Sigma = resultats$Sigma[nrow(Z),,]
        outliers = resultats$outlier_labels
      }
      
      
      print(paste("i = ",i))
      print(paste("compt = ",compt))
      print(paste("j = ",j))
      #print("m = ",m)
      erreurNormeFrobeniusk1l1rho1[j,compt,i,s,t] = norm(Sigma - Sigma1,"F")  
      
      
      
      tc <- table(data$labelsVrais[1:(nrow(Z))], as.numeric(outliers)[1:(nrow(Z))])
      if ("0" %in% rownames(tc)){
        if((tc["0","0"] + tc["0","1"]) != 0)
        {faux_positifsk1l1rho1[j,compt,i,s,t]   <-(tc["0", "1"]/(tc["0", "1"] + tc["0", "0"]))*100}
      }
      if ("1" %in% rownames(tc)){
        if((tc["1","0"] + tc["1","1"]) != 0)
        {faux_negatifsk1l1rho1[j,compt,i,s,t]   <-(tc["1", "1"]/(tc["1", "1"] + tc["1", "0"]))*100}}
      compt = compt + 1
    } 
  }
  }
}
}





#########################################################
#Seuil
########################################################

calcule_fp = function(corr = TRUE){

  fp_seuil =array(0,dim = c(length(contamin_rate), length(methodes)))
  
for (k in seq_along(contamin_rate))
{
 r = 0 
 compt = 1
  data <- genererEchantillon(n,d,mu1,mu2 = k*mu1,p1 = 1- r/100,r/100,Sigma1,Sigma2 = Sigma1,contamin = "moyenne",cluster = FALSE)
  Z = data$Z
  for (m in methodes){
    
    
    
    
    if(m == "offline"){
      
      resultats = OfflineOutlierDetection(Z)
      
      mu_hat = resultats$median
      Sigma = resultats$variance
      distances = resultats$distances
      if (corr == TRUE){
      cutoffCorr = qchisq(.95,df = d)*median(resultats$distances)/qchisq(.5,df = d)
      
      resultats = OfflineOutlierDetection(Z,cutoff = cutoffCorr)
      }
      outliers = resultats$outlier_labels
      
      
    }
    
  if(m == "online"){
    print (m)
    resultats = StreamingOutlierDetectionCorrectedThreshold(Z,batch = 1)
    Sigma = resultats$Sigma[nrow(Z),,]
    distance = resultats$distances
    outliers = resultats$outlier_labels
    
    if(corr == TRUE){
    cutoffCorr = rep(0,nrow(Z))
    outliers = rep(0,nrow(Z))
    for (s in (1 : nrow(Z)))
    {
      cutoffCorr[s]  = qchisq(.95,df = d)*median(resultats$distances[1:s])/qchisq(.5,df = d)
      if (distance[s] > cutoffCorr[s]) {outliers[s] = 1}
    }}
    #cutoffCorr = qchisq(.95,df = d)*median(resultats$distances)/qchisq(.5,df = d)
    #resultats = StreamingOutlierDetection(Z,batch = 1,cutoff = cutoffCorr)
    #outliers = resultats$outlier_labels
   
  }
  if(m == "streaming"){
    print (m)
    resultats = StreamingOutlierDetectionCorrectedThreshold(Z,batch = ncol(Z))
    Sigma = resultats$Sigma[nrow(Z),,]
    outliers = resultats$outlier_labels
    #cutoffCorr = qchisq(.95,df = d)*median(resultats$distances)/qchisq(.5,df = d)
    #resultats = StreamingOutlierDetection(Z,batch = 1,cutoff = cutoffCorr)
    distance = resultats$distances
    
    outliers = resultats$outlier_labels
    if(corr == TRUE){
    cutoffCorr = rep(0,nrow(Z))
    outliers = rep(0,nrow(Z))
    for (s in (1 : nrow(Z)))
    {
      cutoffCorr[s]  = qchisq(.95,df = d)*median(resultats$distances[1:s])/qchisq(.5,df = d)
      if (distance[s] > cutoffCorr[s]) {outliers[s] = 1}
    }}
    
  }
  
  
  tc <- table(data$labelsVrais[1:(nrow(Z))], as.numeric(outliers)[1:(nrow(Z))])
  if ("0" %in% rownames(tc)){
    if((tc["0","0"] + tc["0","1"]) != 0)
    {fp_seuil[k,compt]   <-(tc["0", "1"]/(tc["0", "1"] + tc["0", "0"]))*100
    print(fp_seuil[k,compt])}
  }
  # if ("1" %in% rownames(tc)){
  #   if((tc["1","0"] + tc["1","1"]) != 0)
  #   {fp_seuil[k,compt]   <-(tc["1", "1"]/(tc["1", "1"] + tc["1", "0"]))*100
  #   print(fp_seuil[k,compt])}}
   compt = compt + 1}
}
return(fp_seuil = fp_seuil)
}  
  
fp_seuil = calcule_fp(corr = FALSE)


# Paramètres
alpha <- 0.05
conf_level <- 1 - alpha
p_seq <- seq(0.6, 1, length.out = 200)
N <- nrow(Z)  # N doit être défini au préalable

# Calcul des n0 et des x théoriques
n0 <- floor(p_seq * N)
x <- round(n0 * alpha)

# Calcul de l’intervalle de confiance exact (Clopper-Pearson)
conf <- binom.confint(x = x, n = n0, conf.level = conf_level, methods = "exact")


# Tracé des bornes
plot(p_seq, conf$lower, type = "l", col = "blue", lwd = 2,
     ylim = c(0, 0.15), xlab = "p (proportion pour n0 = p * N)",
     ylab = "Intervalle de confiance de la proportion",
     main = "IC à 95% pour Bin(n0, alpha = 0.05)")
lines(p_seq, conf$upper, col = "red", lwd = 2)

fp_positions = c(1.00, 0.95, 0.9,  0.85, 0.80, 0.75, 0.70, 0.65,0.60)
length(fp_positions)
fp_online_values = fp_seuil[,1]/100
points(fp_positions, fp_online_values, col = "darkgreen", pch = 19)


fp_streaming_values = fp_seuil[,2]/100
points(fp_positions, fp_streaming_values, col = "orange", pch = 19)


fp_offline_values = fp_seuil[,3]/100
points(fp_positions, fp_offline_values, col = "black", pch = 19)


# Configuration de la fenêtre graphique 
par(mfrow = c(2, 2))

# --- Graphique 1 : Online ---
plot(p_seq, conf$lower, type = "l", col = "blue", lwd = 2,
     ylim = c(0, 0.15), 
     xlab = "p (proportion pour n0 = p * N)",
     ylab = "Taux de faux positifs",
     main = "Online Streaming")
lines(p_seq, conf$upper, col = "red", lwd = 2)
points(fp_positions, fp_online_values, col = "darkgreen", pch = 19)
#legend("topright", legend = c("IC inférieur", "IC supérieur", "Observé Online"),
 #      col = c("blue", "red", "darkgreen"), lwd = c(2, 2, NA), pch = c(NA, NA, 19))

# --- Graphique 2 : Streaming ---
plot(p_seq, conf$lower, type = "l", col = "blue", lwd = 2,
     ylim = c(0, 0.15),
     xlab = "p (proportion pour n0 = p * N)",
     ylab = "Taux de faux positifs",
     main = "Streaming")
lines(p_seq, conf$upper, col = "red", lwd = 2)
points(fp_positions, fp_streaming_values, col = "orange", pch = 19)
# legend("topright", legend = c("IC inférieur", "IC supérieur", "Observé Streaming"),
#        col = c("blue", "red", "orange"), lwd = c(2, 2, NA), pch = c(NA, NA, 19))

# --- Graphique 3 : Offline ---
plot(p_seq, conf$lower, type = "l", col = "blue", lwd = 2,
     ylim = c(0, 0.15),
     xlab = "p (proportion pour n0 = p * N)",
     ylab = "Taux de faux positifs",
     main = "Offline")
lines(p_seq, conf$upper, col = "red", lwd = 2)
points(fp_positions, fp_offline_values, col = "black", pch = 19)
# legend("topright", legend = c("IC inférieur", "IC supérieur", "Observé Offline"),
#        col = c("blue", "red", "black"), lwd = c(2, 2, NA), pch = c(NA, NA, 19))

# Titre global
mtext("Comparaison des taux de faux positifs - Streaming Online et Offline après correction", 
      outer = TRUE, cex = 1.2, font = 2)

