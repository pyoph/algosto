#install.packages("microbenchmark")
#install.packages("matlib")
#install.packages("MASS")
#install.packages("reshape2")
#install.packages("corrplot")
#install.packages("dlpyr")
#install.packages("fbm")
#install.packages("mvnfast")
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
#library("mvnfast")


#Génération d'un échantillon de n vecteurs de taille d selon différents modes de contamination
#Génération d'un échantillon de n vecteurs de taille d selon différents modes de contamination

genererEchantillon <- function(n, d, mu1, mu2, p1, p2, Sigma1, Sigma2, contamin = "moyenne",deltacarre = 1.2, k = 5,rhomult = 0.9) {
  # Initialisation
  n1 <- floor(p1 * n)  # Taille du groupe non contaminé
  n2 <- n - n1         # Taille du groupe contaminé
  
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
      vecteurs_mu2 <- matrix(runif(n2 * d, min = -2, max = 2), nrow = n2, ncol = d)
    } 
    else if (contamin == "zero") {
      vecteurs_mu1 <- mvrnorm(n1, mu1, Sigma1)
      vecteurs_mu2 <- matrix(0, nrow = n2, ncol = d)
    } 
    else if (contamin == "MaronnaZamar")
    {
      vecteurs_mu1 <- mvrnorm(n1, rep(0,d), diag(d)) #Non outliers
      
      #Calcul de a0
      a0 <- rep(0,d)
      for(i in 1:d)
      {
        a0[i] <- runif(1,0,1)
      }
      #Normalisation de a0
      a0 <- (a0 - mean(a0))/sqrt(sum((a0 - mean(a0))^2))
      vecteurs_mu2 <- mvrnorm(n2,k*a0,deltacarre*diag(d)) #Outliers
      
      #Création de la matrice de Toeplitz pour introduire de la corrélation : 
      
      R <- matrix(0,d,d)
      
      for (i in 1:d) 
      {
        for (j in 1:d)
        {
          if (i == j) {R[i,i] = 1}
          else R[i,j] = rhomult
        }
      }
      
      #R <- creerMatriceToeplitz(rhomult,d)
      #Calcul des nouvelles matrices corrélées
      vecteurs_mu1 <- R %*% t(vecteurs_mu1)
      vecteurs_mu1 <- t(vecteurs_mu1)
      vecteurs_mu2 <- R %*% t(vecteurs_mu2)
      vecteurs_mu2 <- t(vecteurs_mu2)
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
  #Vvrai <- WeiszfeldCov(Z, nitermax = 1000000)$covmedian
  #VcovVrai <- GmedianCov(Z, scores = 10)
  #VpvraiesV <- eigen(VcovVrai$covmedian)$values
  
  # Retourner les résultats
  return(list(Z = Z,labelsVrais = labelsVrais))
}




