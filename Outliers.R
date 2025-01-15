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


#Calcul des distances de Mahalanobis à partir de la vraie covariance 
calcule_vecteur_distancesEmpiriqueVrai <- function(Z,Sigma)
  
{
  
  distances <- rep(0,nrow(Z))
  CovEmp <- Sigma1
  #CovEmp <- solve(cov(Z))%*%CovEmp
  
  for  (i in 1:nrow(Z))
  {
    #CovEmp <- (t(Z[i,])%*%Z[i,])
    distances[i] = as.numeric(t(Z[i,])) %*% solve(CovEmp) %*% (as.numeric(Z[i,]))
    #distances[i] = as.numeric(Z[i,])%*%as.numeric(Z[i,])
  }
  
  return (distances)
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

#Calcul des distances de Mahalanobis à partir de la médiane géométrique des vecteurs propres et des valeurs propres estimées
calcule_vecteur_distancesOnline <- function(Z,miter,U,lambda)
{

  distances <- rep(0,nrow(Z))

for (i in 1:(n-1))
{
  U[i,,] <- U[i,,] %*% diag(1/sqrt(colSums(U[i,,]^2)))

  # for (j in (1:d))
  #  {
  #   SigmaEstim <- 1/lambda[i,j]*U[i,,j]%*%t(U[i,,j])
  #}
  Sigma <- (U[i,,]) %*% diag(lambda[i,])%*% t(U[i,,])
  #SigmaEstim <- diag(lambda[i,])
  distances[i] = as.numeric(Z[i,] - miter[i]) %*% solve(Sigma) %*% (as.numeric(t(Z[i,] - miter[i])))

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


detectionOutliers <- function(distances, n,d,cutoff =qchisq(p = 0.95,df = d))
  
{
  
  outliers_labels <- rep(0,n)
  
  #cutoff <- qchisq(p = 0.975, df = ncol(Z))
  
  
  
  
  
  #m = rep(0,d)
  for (i in 1:n)
  {
    
    #print(i)
    
    if (distances[i] > cutoff) {outliers_labels[i] <- 1}
    else {outliers_labels[i] <- 0}
    #}
  }
  
  #print(which(outliers_labels == 1))
  return (outliers_labels)
  
}


