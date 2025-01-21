#install.packages("microbenchmark")
#install.packages("matlib")
#install.packages("MASS")
#install.packages("reshape2")
#install.packages("corrplot")
#install.packages("dlpyr")
#install.packages("devtools")
#install.packages("usethis")

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


###Initialisation des paramètres communs à toutes les simulations


#Permuter lignes colonnes Toeplitz

# Fonction pour permuter des lignes et des colonnes dans une matrice
permuterLignesColonnes <- function(matrice, lignes_a_permuter = c(1, 2), colonnes_a_permuter = c(1, 2)) {
  
  # Permuter les lignes
  matrice[lignes_a_permuter, ] <- matrice[rev(lignes_a_permuter), ]
  
  # Permuter les colonnes
  matrice[, colonnes_a_permuter] <- matrice[, rev(colonnes_a_permuter)]
  
  # Retourner la matrice modifiée
  return(matrice)
}


#Estimation de Sigma critère et détection outliers

#Mettre l'identité dans la diagonale de Toeplitz
creerMatriceToeplitz <- function(rho,d)
{
  
  # Construire la matrice de Toeplitz
  toeplitz_matrix <- matrix(0, nrow = d, ncol = d)
  
  # Remplir la matrice avec rho^{|i - j|}
  for (i in 1:d) {
    for (j in 1:d) {
      toeplitz_matrix[i, j] <- rho^abs(i - j)
    }
  }
  
  
  
    return(diag(sqrt(1:d))%*%toeplitz_matrix%*%diag(sqrt(1:d)))
  #return(toeplitz_matrix)
}





#Initialisation des paramètres communs à toutes les simulations

nbruns = 1
n = 1e4
d = 10
mu1 = rep(0,d)
mu2 = 2*rep(1,d)
rho = 0.8
Sigma1 = creerMatriceToeplitz(rho,d)

lignes_a_permuter = c(1, 2)
colonnes_a_permuter = c(1, 2)
Sigma2 = permuterLignesColonnes(Sigma1,lignes_a_permuter , colonnes_a_permuter)
#Sigma2 = 0.2* Sigma1
