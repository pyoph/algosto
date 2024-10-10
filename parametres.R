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
#Dimension



initialiser_parametres <- function(d,q,n,c,m0,Sigma1,Sigma2,p1,p2,mu1,mu2,r,k) {
  # Initialisation des paramètres
  params <- list(
    d = d   ,              # Dimension des données mettre ncol(Y)
    n = n,             # Taille de l'échantillon nrow(Y)
    q = q,              # Nombre de plus grandes valeurs propres
    c = c,        # Taux d'apprentissage
    m = m0,      # Moyenne réelle (vecteur de zéros)
    Sigma1 = Sigma1,        # Covariance des non outliers
    Sigma2 = Sigma2,    # Covariance des outliers
    p1 = p1,            # Proportion de non outliers
    p2 = p2,            # Proportion d'outliers
    mu1 = mu1,    # Moyenne des non outliers
    mu2 = mu2, #a Moyenne des outliers
    r = r, #Paramètre r pour générer les variables aléatoires
    k = k #Taille des blocs
      )

  return(params)
}

