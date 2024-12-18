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


n <- 1e4

d <- 10

rho <- 0.7

seuil_p_value <- 0.05

Sigma1 <- creerMatriceToeplitz(rho,d)

Sigma2 <- creerMatriceToeplitz(0.4,d)
#Sigma1 <- diag(sqrt(1:d))
#Sigma1 <- diag(d)
mu1 <- rep(0,d)

mu2 <- 5*rep(1,d)




p1 <- 0.98
p2 <- 1 - p1


resultsSimul <- genererEchantillon(n, d, mu1, mu2, p1, p2, Sigma1 = Sigma1, Sigma2 = Sigma2)

#cov(Z)

Z <- resultsSimul$Z




calcule_vecteur_distances <- function(Z,m,Sigma)
  
{
  
  distances <- rep(0,n)
  
  for  (i in 1:nrow(Z))
  {
    distances[i] = as.numeric(Z[i,] - m) %*% solve(Sigma) %*% (as.numeric(t(Z[i,] - m)))
    
  }
  
  return (distances)
}



# Fonction pour calculer le cutoff
calcule_cutoff <- function(distances, d) {
  
  # Constantes et quantiles
  quantile_95 <- qchisq(0.95, df = d) 
  quantile_50 <- qchisq(0.5, df = d)    
  median_d <- median(distances)        
  
  # Calcul du cutoff
  cv <- (quantile_95 * median_d) / quantile_50
  
  return(cv)
}


