
######################################
################Packages nécessaires#####
#########################################
setwd("~/algosto")
packages = c("Rcpp","Gmedian","MASS","DescTools" ,"capushe","checkmate", "doFuture", "future",'mclust', 'LaplacesDemon', 'genieclust', 'reshape2','cowplot','scales',"bookdown","xfun","dplyr","binom","pROC","mclust")
#
for (p in packages) {
   if (!requireNamespace(p, quietly = TRUE)) {
     install.packages(p)
   }
   library(p, character.only = TRUE)
 }
#
setwd("~/algosto")

packages_us = c("STARRS_1.0.tar.gz")

for (p in packages_us) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
  
}

library(Rcpp)
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
  # Initialisation
  Z = matrix(0,n,d)
  idx_out <- id_outliers
  idx_in  <- setdiff(1:n, idx_out)
  
  n1 <- length(idx_in)
  n2 <- length(idx_out)
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


calcule_cutoff = function(distances,c_m = 2, n = 1e4,cutinit = 0.6,cutoff = .95,type = "raw",nbcol = 10)
{
  
  if(type == "raw"){return (qchisq(.95,df = nbcol))}
  
  if(type == "quantcorr"){
    cutoffcor=quantile(distances[1:n],probs=cutinit)
    c_med=median(distances[1:n])
    for (i in 1 :n) {
      cutoffcor <- cutoffcor - c_m*c_med*(i+1)^(-0.75)*( as.numeric(distances[i] <= cutoffcor) - cutoff)
    }
    return(cutoffcor)}
  
  if(type == "rescale_dist"){
    
    return(qchisq(.95,df = nbcol) * median(distances)/qchisq(0.5, df = nbcol))
    
  }
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
      cutoffcor <- cutoffcor - c_m*c_med*(j+1)^(-0.75)*
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
    else{
      cutoffcor <- cutoffcor -
        c_m*c_med*(n0 + (niterr-1) + j)^(-0.75)*
        (as.numeric(S <= cutoffcor) - cutoffquant)
      
      if(S > cutoffcor){
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


#################SMD data preprocessing functions###########################"


load_smd_machine <- function(test_path, label_path) {
  
  Z <- as.matrix(
    read.table(test_path, sep = ",")
  )
  
  labels <- scan(label_path)
  
  list(
    Z = Z,
    labels = labels,
    name = basename(test_path)
  )
}





remove_high_corr <- function(Z, threshold = 0.9999)
{
  C <- cor(Z)
  
  C[lower.tri(C, diag = TRUE)] <- 0
  
  remove <- which(apply(abs(C) > threshold, 2, any))
  
  if(length(remove) > 0)
  {
    cat("Removing cols:", remove, "\n")
  }
  
  Z[, -remove, drop = FALSE]
  return(Z)
}

sequential_clean_columns <- function(Z, begin = 1, cut = NULL) {
  
  Z <- as.matrix(Z)
  p <- ncol(Z)
  
  if (is.null(cut)) cut <- 0.95
  
  start <- NULL
  
  # =========================
  # 1) trouver prefix stable
  # =========================
  
  for (k in (begin + 1):p) {
    
    ok <- tryCatch({
      
      onlineRobustVariance(
        Z[, begin:k, drop = FALSE],
        computeOutliers = TRUE,
        cutoff = cut,
        cutinit = 0.6,
        nDataInit = 1e3,
        c_m = 1,
        batch = 3
      )
      
      TRUE
      
    }, error = function(e) {
      FALSE
    })
    
    if (ok) {
      start <- k
      cat("✔ first stable prefix found at:", k, "\n")
      break
    }
  }
  
  if (is.null(start)) {
    stop("No stable prefix found")
  }
  
  keep <- c(begin, start)
  
  for (j in (start + 1):p) {
    
    cat("j =", j, "keep =", keep, "\n")
    
    ok <- tryCatch({
      # 
      onlineRobustVariance(
        Z[, c(keep, j), drop = FALSE],
        computeOutliers = TRUE,
        cutoff = cut,
        cutinit = 0.6,
        nDataInit = 1e3,
        c_m = 1,
        batch = 3
      )
      TRUE
      
    }, error = function(e) FALSE)
    
    if (ok) {
      keep <- c(keep, j)
      cat("✔ col:", j, "OK\n")
    } else {
      cat("✘ col:", j, "KO\n")
    }
  }
  
  return(list(
    Z = Z[, keep, drop = FALSE],
    keep = keep
  ))
}


KL <- function(parms1, parms2){
  invSigma2 <- solve(parms2$Sigma)
  0.5*(log(det(parms2$Sigma)/det(parms1$Sigma)) - d + sum(diag(invSigma2%*%parms1$Sigma)) +
         t(parms2$mu-parms1$mu)%*%invSigma2%*%(parms2$mu-parms1$mu))[1, 1]
}

FrobeniusNormError = function(Sigmahat,Sigma){
  return (norm(Sigmahat - Sigma,"F"))
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

ARI_manual <- function(true_labels, pred_labels) {
  
  n <- length(true_labels)
  
  # matrices de comparaison paire-à-paire
  comb <- function(x) outer(x, x, FUN = "==")
  
  TP <- sum(comb(true_labels) & comb(pred_labels)) / 2
  TN <- sum(!comb(true_labels) & !comb(pred_labels)) / 2
  FP <- sum(!comb(true_labels) & comb(pred_labels)) / 2
  FN <- sum(comb(true_labels) & !comb(pred_labels)) / 2
  
  index <- TP + TN
  expected <- ((TP + FP) * (TP + FN) + (FN + TN) * (FP + TN)) / choose(n, 2)
  max_index <- (TP + FP + FN + TN + index) / 2
  
  ARI <- (index - expected) / (max_index - expected)
  
  return(ARI)
}


auc_manual <- function(score, labels) {
  
  thresholds <- sort(unique(score))
  
  TPR <- numeric(length(thresholds))
  FPR <- numeric(length(thresholds))
  
  for (i in seq_along(thresholds)) {
    
    th <- thresholds[i]
    pred <- ifelse(score >= th, 1, 0)
    
    TP <- sum(pred == 1 & labels == 1)
    TN <- sum(pred == 0 & labels == 0)
    FP <- sum(pred == 1 & labels == 0)
    FN <- sum(pred == 0 & labels == 1)
    
    TPR[i] <- TP / (TP + FN + 1e-12)
    FPR[i] <- FP / (FP + TN + 1e-12)
  }
  
  ord <- order(FPR, TPR)
  
  FPR <- c(0, FPR[ord], 1)
  TPR <- c(0, TPR[ord], 1)
  
  auc <- sum(diff(FPR) * (head(TPR, -1) + tail(TPR, -1)) / 2)
  
  return(list(
    FPR = FPR,
    TPR = TPR,
    auc = auc
  ))
}

#####Calcule critères

compute_criteres = function(variance,outlab,distances,labels_vrais,SigmaTrue = Sigma0,r){
  
  erreurFrob <- norm(variance - SigmaTrue, "F")
  FP <- sum(outlab == 1 & labels_vrais == 0)
  FN <- sum(outlab == 0 & labels_vrais == 1)
  if(r != 0){
    
    # ARI
     ari <- adjustedRandIndex(
       labels_vrais,
       outlab
     )
  
    
    #ari = ARI_manual(labels_vrais,outlab)  
    print(paste0("ARI ",ari))
      # AUC
    auc_val <- as.numeric(auc(roc(labels_vrais, distances, direction='<')))
    
    #auc_val <- auc_manual(as.numeric(distances),labels_vrais)$auc
    
    print(paste0("AUC ",auc_val))
    
  }
  else{
    auc_val = .5
    ari = 0
  }
  
  return(list(
    erreurFrob = erreurFrob,
    FP = FP,
    FN = FN,
    ARI = ari,
    AUC = auc_val  ))
}





