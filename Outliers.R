#install.packages("microbenchmark")
#install.packages("matlib")
#install.packages("MASS")
#install.packages("reshape2")
#install.packages("corrplot")
#install.packages("dlpyr")
library(reshape2)
library(RobRegression)
library(Gmedian)
library(ggplot2)
library(far)
library(gridExtra)
library(microbenchmark)
library("matlib")
library("MASS")
library("corrplot")
library("dplyr")


#Calcule le vecteur distances selon la méthode choisie
<<<<<<< HEAD

=======
>>>>>>> 661875862712b4f7546a18cb504c4e554d78b4e1
calcule_distances_par_methode <- function(Z, methode) {
  distances <- NULL
  
  if (methode == "Shrinkage") {
    # Méthode Shrinkage
    med <- covComed(Z)$center
    SigmaShrink <- covCor(Z)
    distances <- calcule_vecteur_distances(Z, med, SigmaShrink)
    
  } else if (methode == "Comédiane") {
    # Méthode Comédiane
    med <- covComed(Z)$center
    SigmaComed <- covComed(Z)$cov
    distances <- calcule_vecteur_distances(Z, med, SigmaComed)
    
  } else if (methode == "OGK") {
    # Méthode OGK
    OGK_result <- covOGK(Z,sigmamu = s_mad)
    med <- OGK_result$center
    SigmaOGK <- OGK_result$cov
    distances <- calcule_vecteur_distances(Z, med, SigmaOGK)
    
  } else if (methode == "Cov Empirique") {
    # Méthode Empirique
    med <- colMeans(Z)
    SigmaEmp <- cov(Z)
    distances <- calcule_vecteur_distances(Z, med, SigmaEmp)
    
  } else if (methode == "Offline") {
    # Méthode Offline
    med <- RobVar(Z)$median
    SigmaOffline <- RobVar(Z)$variance
    distances <- calcule_vecteur_distances(Z, med, SigmaOffline)
    
  } else if (methode == "Online") {
    # Méthode Online
    results <- estimMV(Z)
    #miter <- results$miter
    #U <- results$U
    #lambda <- results$lambdaIter
    distances <- results$distances
<<<<<<< HEAD
    
=======
  
>>>>>>> 661875862712b4f7546a18cb504c4e554d78b4e1
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


