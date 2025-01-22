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


