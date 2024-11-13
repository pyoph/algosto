#install.packages("microbenchmark")
#install.packages("matlib")
#install.packages("MASS")
#install.packages("reshape2")
#install.packages("corrplot")
#install.packages("dlpyr")
#install.packages("fbm")
#install.packages("mvnfast")
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
library("mvnfast")

genererEchantillon <- function(n,d,mu1,mu2,p1,p2,Sigma1,Sigma2)
{


  #S'il y a des outliers
  if (p2 > 0) {
  # Calculer les tailles des sous-échantillons
  n1 <- floor(p1 * n)
  n2 <- n - n1


  # Calculer les tailles des sous-échantillons
  n1 <- floor(p1 * n)
  n2 <- n - n1

  # Générer les vecteurs gaussiens
  vecteurs_mu1 <- mvrnorm(n1, mu1, Sigma1)
  #vecteurs_mu2 <- mvrnorm(n2, mu2, Sigma1)
  vecteurs_mu2 <- rmvt(n2, mu1,diag(d),df = 1, ncores = 1, A = NULL)
  
  # Ajouter des labels
  labels_mu1 <- rep(0, n1)  # Labels pour les vecteurs avec moyenne mu1
  labels_mu2 <- rep(1, n2)  # Labels pour les vecteurs avec moyenne mu2
  labels <- c(labels_mu1, labels_mu2)  # Combinaison des labels

  # Créer une matrice Z combinant les vecteurs de mu1 et mu2
  Z1 <- rbind(vecteurs_mu1, vecteurs_mu2)

  # Créer un data frame en associant Z et les labels
  data <- data.frame(Z = Z1, label = labels)

  # Mélanger aléatoirement les lignes tout en gardant les labels associés
  set.seed(123)  # Pour garantir la reproductibilité
  data_melangee <- data[sample(nrow(data)), ]

  # Récupérer les données et les labels mélangés
  Z <- as.matrix(data_melangee[, 1:d])  # Récupérer uniquement les colonnes Z (d dimensions)
  labelsVrais <- data_melangee$label    # Récupérer les labels mélangés

  }
  else {
    Z <- mvrnorm(n, mu1, Sigma1)
      # Labels pour les vecteurs avec moyenne mu
      labelsVrais <- rep(0, n)
    }

  #Génération des vraies valeurs de la MCM
  Vvrai <- WeiszfeldCov(Z,nitermax = 1000000)$covmedian

  VcovVrai <- GmedianCov(Z, scores = 10)
  #Génération des vraies valeurs propres de la MCM
  VpvraiesV <- eigen(VcovVrai$covmedian)$values

  return(list(Z = Z,Vvrai = Vvrai,VcovVrai = VcovVrai,VpvraiesV = VpvraiesV,labelsVrais = labelsVrais))

}
