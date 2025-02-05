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


#Génération d'un échantillon de n vecteurs de taille d selon différents modes de contamination

genererEchantillon <- function(n, d, mu1, mu2, p1, p2, Sigma1, Sigma2, contamin = "moyenne") {
  # Initialisation
  n1 <- floor(p1 * n)  # Taille du groupe majoritaire
  n2 <- n - n1         # Taille du groupe minoritaire
  
  labels_mu1 <- rep(0, n1)  # Labels pour les vecteurs avec moyenne mu1
  labels_mu2 <- rep(1, n2)  # Labels pour les vecteurs avec moyenne mu2
  
  if (p2 > 0) {
    # Générer les vecteurs selon le type de contamination
    if (contamin == "moyenne") {
      vecteurs_mu1 <- mvrnorm(n1, mu1, Sigma1)
      vecteurs_mu2 <- mvrnorm(n2, mu2, Sigma1)
    } else if (contamin == "variance") {
      vecteurs_mu1 <- mvrnorm(n1, mu1, Sigma1)
      vecteurs_mu2 <- mvrnorm(n2, mu1, Sigma2)
    } else if (contamin == "student") {
      vecteurs_mu1 <- mvrnorm(n1, mu1, Sigma1)
      vecteurs_mu2 <- rmvt(n2, mu1, sigma = Sigma1, df = 1, ncores = 1, A = NULL)
    } else if (contamin == "uniforme") {
      vecteurs_mu1 <- mvrnorm(n1, mu1, Sigma1)
      vecteurs_mu2 <- matrix(runif(n2 * d, min = -50, max = 50), nrow = n2, ncol = d)
    } 
    else if (contamin == "zero") {
      vecteurs_mu1 <- mvrnorm(n1, mu1, Sigma1)
      vecteurs_mu2 <- matrix(0, nrow = n2, ncol = d)
    } 
    
    else {
      stop("Type de contamination non reconnu.")
    }
    
    # Combinaison des vecteurs
    Z <- rbind(vecteurs_mu1, vecteurs_mu2)
    labelsVrais <- c(labels_mu1, labels_mu2)
  } else {
    # Pas de contamination
    Z <- mvrnorm(n, mu1, Sigma1)
    labelsVrais <- rep(0, n)
  }
  
  # Mélanger aléatoirement les données
  set.seed(123)  # Pour garantir la reproductibilité
  indices <- sample(nrow(Z))
  Z <- Z[indices, ]
  labelsVrais <- labelsVrais[indices]
  
  # Calcul des vraies valeurs pour la médiane géométrique
  Vvrai <- WeiszfeldCov(Z, nitermax = 1000000)$covmedian
  VcovVrai <- GmedianCov(Z, scores = 10)
  VpvraiesV <- eigen(VcovVrai$covmedian)$values
  
  # Retourner les résultats
  return(list(Z = Z, Vvrai = Vvrai, VcovVrai = VcovVrai, VpvraiesV = VpvraiesV, labelsVrais = labelsVrais))
}


