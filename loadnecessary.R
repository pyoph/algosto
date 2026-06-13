
######################################
################Packages nécessaires#####
#########################################
setwd("~/algosto")
packages = c("Rcpp","Gmedian","MASS","DescTools" ,"capushe","checkmate", "doFuture", "future",'mclust', 'LaplacesDemon', 'genieclust', 'reshape2','cowplot','scales',"bookdown")
#
for (p in packages) {
   if (!requireNamespace(p, quietly = TRUE)) {
     install.packages(p)
   }
   library(p, character.only = TRUE)
 }
#
setwd("~/algosto")
install.packages("RMM_1.0.tar.gz",repos = NULL,type = "source")
install.packages("binom")
install.packages("STARRS_1.0.tar.gz")

library(Rcpp)
library(RMM)
library(Gmedian)
library(MASS)
library(ggplot2)
library(mvtnorm)
library(binom)
library(reshape2)
library(dplyr)
library(tidyr)
library(scales)
library("robustbase")
library(STARRS)
library(mclust)
library(pROC)

############################
#############FIchiers nécessaires
#################################
source("~/algosto/algorithmes.R")
sourceCpp("~/algosto/src/CodesRCpp.cpp")
source('~/algosto/FunctionsKLgauss.R')

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


genererEchantillon_new <- function(n, d, mu1, mu2,Sigma1, Sigma2,r,id_outliers = NULL,dirac = FALSE) {
  # InitialSisation
  Z = matrix(0,n,d)
  
  n1 <- floor((1 -r/100) * n)  # Taille du groupe non contaminé
  n2 <- n - n1         # Taille du groupe contaminé
  labelsVrais = rep(0,n)
  # labels_mu1 <- rep(0, n1)  # Labels pour les vecteurs avec moyenne mu1
  # labels_mu2 <- rep(1, n2)  # Labels pour les vecteurs avec moyenne mu2
  # 
  #print(paste0("n1 = ",n1))
 if(r > 0){ 
  # inliers
  idx_in <- setdiff(1:n, id_outliers)
 # print(paste0("length idx_in = ",length(idx_in)))
#  print(paste0("dim(Z) = ",dim(Z)[2]))
  
  Z[idx_in, ] <- mvrnorm(n1,mu1, Sigma1)
  
  labelsVrais[idx_in] <- 0
  if(dirac == FALSE){
  # outliers
  Z[id_outliers, ] <- mvrnorm(n2,mu2, Sigma2)}
  else if (dirac == TRUE){Z[id_outliers, ] = mu2}
  labelsVrais[id_outliers] <- 1
 }
  # 
  # if (r > 0) {
  #   # Générer les vecteurs selon le type de contamination
  #   for (k in (1:n)){
  #     if (k %in% id_outliers)
  #       {
  #      Z[k,] = mvrnorm(1, mu2, Sigma2)
  #      labelsVrais[k] = 1
  #     }
  #     else {
  #       Z[k,] = mvrnorm(1, mu1, Sigma1)
  #     }
  #   }
  #   
  #   
    # # Combinaison des vecteurs
    # Z <- rbind(vecteurs_mu1, vecteurs_mu2)
    # labelsVrais <- c(labels_mu1, labels_mu2)  
  
  
  
  
  else if (r == 0){
    # Pas de contamination
    Z <- mvrnorm(n, mu1, Sigma1)
    #labelsVrais <- rep(0, n)
  } # Mélanger aléatoirement les données
#  set.seed(123)  # Pour garantir la reproductibilité
  #indices <- sample(nrow(Z))
  #Z <- Z[indices, ]
  #labelsVrais <- labelsVrais[indices]
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

#Inverse with Shermann Morisson formula

inverse_Schermann_Morisson = function(A,Z) 
{
  nbcol = ncol(A)
  invA = diag(nbcol)
  
  
  niter = nrow(Z)
  
  
  for (i in (1:(niter-1)))
  {
    scal = 1 + t(Z[i+1,])%*%invA%*%Z[i+1,]
    if(scal !=0){
      invA = invA - 1/scal[1]*invA%*%(Z[i+1,]%*%t(Z[i+1,]))%*%invA}
    
  }
  
  
  return (invA)
  
}


calcule_cutoff = function(distances,c_m = 2, n,cutinit = 0.6,cutoff = .9)
{
  cutoffcor=quantile(distances[1:n],probs=cutinit)
  c_med=median(distances[1:n])
  for (i in 1 :n) {
    cutoffcor <- cutoffcor - c_m*c_med*(i+1)^(-0.75)*( as.numeric(distances[i] <= cutoffcor) - cutoff)
  }
  return(cutoffcor)
}
inv_safe <- function(S, eps = 1e-4) {
  
  S <- as.matrix(S)
  invS = diag(ncol(S))
  # sécurité NA / Inf
  if (any(!is.finite(S))) {
    stop("Matrix contains NA/Inf")
  }
  
  tryCatch({
    
  invS =   solve(S)
    
  }, error = function(e) {
    
    cat("Singular matrix -> ridge correction applied\n")
    
    invS = solve(S + eps * diag(ncol(S)))
    
  })
  
  return (invS)
}
  # }
# test_outliers = function(distances,dim = 10,cutoff)
# {
#   outlab = rep(0,length(distances))
#   
#   
#   for (i in (1:length(distances)))
#   {
#     if(distances[i] > cutoff){outlab[i] = 1}
#   }
#   return(outlab)
# }


##############################################
####################Online sample covariance
###############################################SA

SampleCovOnline = function(Z,quantcutoff = FALSE,nDataInit=100,cutoffquant =.9,cutinit = .6,c_m = 2)
{
  nblignes = nrow(Z)
  n0 = nDataInit #Nombre données init
  
  # Initialisation
  mean = colMeans(Z[1:n0, , drop = FALSE])
  meanOld = mean
  Sigma = cov(Z[1:n0, , drop = FALSE])
  
  M = n0*mean
  
  A = matrix(0,ncol(Z),ncol(Z))
  
  for (j in (1:n0))
  {
    A = A + (Z[j,] - mean)%*%t(Z[j,] - mean)
  }
  
  B <- NULL
  
  for (lambda in c(0, 1e-10, 1e-8, 1e-6, 1e-4)) {
    
    B <- tryCatch(
      solve(A + lambda * diag(ncol(A))),
      error = function(e) NULL
    )
    
    if (!is.null(B)) break
  }
  
  if (is.null(B))
    stop("Unable to invert covariance matrix.")
  
  # B = solve(A)
  
  invSigma = n0*B
  
  meanIter = matrix(0, nrow(Z), ncol(Z))
  SigmaIter = array(0, dim = c(nrow(Z), ncol(Z), ncol(Z)))
  
  meanIter[1:n0,] = mean
  
  distances = rep(0, nrow(Z))
  outliers_labels = rep(0, nrow(Z))
  
  meanOld = mean
  
  distances = rep(0, nrow(Z))
  outliers_labels = rep(0,nrow(Z))
  
  cutoff = qchisq(.95,df = ncol(Z))
  
  for(j in 1:n0)
  {
    distances[j] = t(Z[j,] - meanIter[j,])%*%invSigma%*%(Z[j,] - meanIter[j,])
  }
  
  if(quantcutoff == TRUE)
  {
    cutoffcor = quantile(distances[1:n0], probs = cutinit)
    c_med=median(distances[1:n0])
    
  }
  
  for (j in 1:n0)
  {
    
    if(quantcutoff == FALSE)
    {
      if(distances[j] > cutoff)
      {
        outliers_labels[j] = 1
      }
    }
    else
    {  
      cutoffcor <- cutoffcor - c_m*c_med*(i+1)^(-0.75)*
        (as.numeric(distances[j] <= cutoffcor) - cutoffquant)
      
      if (distances[j] > cutoffcor)
      {
        outliers_labels[j] <- 1
      }
    }  
  }
  
  niterr = 0 
  
  for (i in ((n0 + 1):(nblignes-1)))
  {
    niterr = niterr + 1
    
    M = M + Z[i,]
    Xbar = M/i
    
    U = Z[i,] - Xbar
    
    scal = 1 + t(U)%*%B%*%U
    
    if(scal !=0)
    {
      B = B - 1/scal[1]*(B%*%U)%*%(t(B%*%U))
    }
    
    invSigma = i*B
    
    mean = mean + (1.0 / (i + 1)) * (Z[i+1,] - mean)
    
    Sigma = (i - 1)/i*Sigma +
      1/(i + 1)*((Z[i+1,] - meanOld)%*%t(Z[i+1,] - meanOld))
    
    meanOld = mean
    
    meanIter[i,] = mean
    SigmaIter[i,,] = Sigma
    
    distances[i] = t(Z[i,] - meanIter[i,])%*%invSigma%*%(Z[i,] - meanIter[i,])
    
    S = distances[i]
    
    
    cutoff = qchisq(.95,df = ncol(Z))
    
    if(quantcutoff == FALSE)
    {
      if (S > cutoff)
      {
        outliers_labels[i] = 1
      }
    }
    else
    {
      cutoffcor <- cutoffcor -
        c_m*c_med*(n0 + (niterr-1) + j)^(-0.75)*
        (as.numeric(S <= cutoffcor) - cutoffquant)
      
      if(S > cutoffcor)
      {
        outliers_labels[i] = 1
      }
    }
  }
  
  SigmaIter[nrow(Z),,] = Sigma
  meanIter[nrow(Z),] = mean
  
  return(list(
    mean = mean,
    Sigma = Sigma,
    meanIter = meanIter,
    SigmaIter = SigmaIter,
    distances = distances,
    outliers_labels = outliers_labels
  ))
}

# --- Transformation pseudo-log  ---
pseudo_log <- function(y) {
  return(log10(1+y))  
}


###############Calcule taux pour trajectoires###############################

compute_rates <- function(pred, labels) {
  
  n <- length(pred)
  
  FP_rate <- rep(0, n)
  FN_rate <- rep(0, n)
  fp <- 0
  tn <- 0
  tp <- 0
  fn <- 0
  
  for (t in 1:n) {
    
    if (pred[t] == 1 && labels[t] == 0) fp <- fp + 1
    if (pred[t] == 0 && labels[t] == 1) fn <- fn + 1
    if (pred[t] == 1 && labels[t] == 1) tp <- tp + 1
    if (pred[t] == 0 && labels[t] == 0) tn <- tn + 1
    
    # vrai taux streaming
    FP_rate[t] <- fp / (fp + tn + 1e-10)
    FN_rate[t] <- fn / (fn + tp + 1e-10)
  }
  
  list(FP_rate = FP_rate,
       FN_rate = FN_rate)
}