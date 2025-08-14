
######################################
################Packages nécessaires#####
#########################################

packages = c("Rcpp","Gmedian","MASS","ggplot2","cowplot","patchwork")

for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
  library(p, character.only = TRUE)
}



library(Rcpp)
library(RMM)
library(Gmedian)
library(MASS)
library(ggplot2)
library(cowplot)
library(patchwork)

############################
#############FIchiers nécessaires
#################################
source("~/algosto/algorithmes.R")
sourceCpp("~/algosto/src/CodesRCpp.cpp")


#######################################################################
#Simulation of a dataset with different parameters####################
######################################################################

genererEchantillon <- function(n, d, mu1, mu2,Sigma1, Sigma2,r) {
  # InitialSisation
  n1 <- floor((1 -r/100) * n)  # Taille du groupe non contaminé
  n2 <- n - n1         # Taille du groupe contaminé
  
  labels_mu1 <- rep(0, n1)  # Labels pour les vecteurs avec moyenne mu1
  labels_mu2 <- rep(1, n2)  # Labels pour les vecteurs avec moyenne mu2
  
  if (r > 0) {
    # Générer les vecteurs selon le type de contamination
    vecteurs_mu1 <- mvrnorm(n1, mu1, Sigma1)
    vecteurs_mu2 <- mvrnorm(n2, mu2, Sigma2)
    
    # Combinaison des vecteurs
    Z <- rbind(vecteurs_mu1, vecteurs_mu2)
    labelsVrais <- c(labels_mu1, labels_mu2)  
  }
  
  
  
  else {
    # Pas de contamination
    Z <- mvrnorm(n, mu1, Sigma1)
    labelsVrais <- rep(0, n)
  } # Mélanger aléatoirement les données
  set.seed(123)  # Pour garantir la reproductibilité
  indices <- sample(nrow(Z))
  Z <- Z[indices, ]
  labelsVrais <- labelsVrais[indices]
  
  
  return(list(Z = Z,labelsVrais = labelsVrais))
}



# Functions
KL <- function(parms1, parms2){
  invSigma2 <- solve(parms2$Sigma)
  0.5*(log(det(parms2$Sigma)/det(parms1$Sigma)) - d + sum(diag(invSigma2%*%parms1$Sigma)) +
         t(parms2$mu-parms1$mu)%*%invSigma2%*%(parms2$mu-parms1$mu))[1, 1]
}

#
# # Contamination parms: F1
ParmsF1 <- function(m1, k1, l1, rho1){
  d <- length(m1)
  mu1 <- k1*m1
  sigmaSq1 <- l1*sigmaSq0
  Sigma1 <- diag(sqrt(sigmaSq1)) %*% toeplitz(rho1^(0:(d-1))) %*% diag(sqrt(sigmaSq1))
  return(list(mu1=mu1, Sigma1=Sigma1))
}


# Functions
KL <- function(parms1, parms2){
  invSigma2 <- solve(parms2$Sigma)
  0.5*(log(det(parms2$Sigma)/det(parms1$Sigma)) - d + sum(diag(invSigma2%*%parms1$Sigma)) +
         t(parms2$mu-parms1$mu)%*%invSigma2%*%(parms2$mu-parms1$mu))[1, 1]
}

test_outliers = function(distances,dim = 10,cutoff)
{
  outlab = rep(0,length(distances))
  
  
  for (i in (1:length(distances)))
  {
    if(distances[i] > cutoff){outlab[i] = 1}
  }
  return(outlab)
}
