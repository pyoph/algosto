#install.packages("microbenchmark")
#install.packages("matlib")
#install.packages("MASS")
#install.packages("reshape2")
#install.packages("corrplot")
#install.packages("dlpyr")
#install.packages("devtools")
#install.packages("usethis")

#library(reshape2)
#library(RobRegression)
#library(Gmedian)
#library(ggplot2)
#library(far)
#library(gridExtra)
#library(microbenchmark)
#library("matlib")
#library("MASS")
#library("corrplot")
#library("dplyr")


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
  
   
  
  #return(diag(sqrt(1:d))%*%toeplitz_matrix%*%diag(sqrt(1:d)))
  return(toeplitz_matrix)
}

  

#Faire un boxplot des résultats enregistrer les résultats à chaque run

#Initialisation des paramètres communs à toutes les simulations
#d = 100 000 
nbruns = 1
n = 1e4
d = 1e2
mu1 = rep(0,d)
mu2 = rep(1/sqrt(d), d)
rho = 0.3
#Sigma1 = creerMatriceToeplitz(rho,d) 
sigmaSq0 <- (1:d); sigmaSq0 <- sigmaSq0 / mean(sigmaSq0)

Sigma1 <- diag(sqrt(sigmaSq0)) %*% toeplitz(rho^(0:(d-1))) %*% diag(sqrt(sigmaSq0))

#Sigma1 = diag(d)
#Sigma1 = diag(c((1:d)))
lignes_a_permuter = c(1, 2)
colonnes_a_permuter = c(1, 2)
#Sigma2 = permuterLignesColonnes(Sigma1,lignes_a_permuter , colonnes_a_permuter)
Sigma2 = 0.1* Sigma1
