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



initialiser_parametres <- function(Y,c,r,k,Sigma = diag(ncol(Y))) {
  # Initialisation des paramètres
  params <- list(
    d = ncol(Y)   ,              # Dimension des données mettre ncol(Y)
    n = nrow(Y),# Taille de l'échantillon nrow(Y)
    Sigma = Sigma,
    q = d,              # Nombre de plus grandes valeurs propres
    c = c,        # Taux d'apprentissage
    #m = m0,      # Moyenne réelle (vecteur de zéros)
    r = r, #Paramètre r pour générer les variables aléatoires
    k = k #Taille des blocs
      )

  return(params)
}

