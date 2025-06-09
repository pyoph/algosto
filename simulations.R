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

genererEchantillon <- function(n, d, mu1, mu2, p1, p2, Sigma1, Sigma2, contamin = "moyenne",deltacarre = 1.2, k = 5,rhomult = 0.9, cluster = FALSE,nclust = d,cutoff = qchisq(0.95, df = d)) {
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
    } 
    else if (contamin == "studentTronquee") {
      vecteurs_mu1 <- mvrnorm(n1, mu1, Sigma1)
      # vecteurs_mu2 <- rmvt(n2, mu1, sigma = Sigma1, df = 1, ncores = 1, A = NULL)
      # vecteurs_mu2 <- ifelse(vecteurs_mu2 > -2 & vecteurs_mu2 < 2, 
      #                        ifelse(vecteurs_mu2 < 0, -2, 2), 
      #                        vecteurs_mu2)
      compt = 0
      vecteurs_mu2 = matrix(0,0, d)
      #Inverser Sigma1 avant
      invSigma1 = solve(Sigma1)
      while (compt < n2) 
      {
        student <- rmvt(1, mu1, sigma = Sigma1, df = 1, ncores = 1, A = NULL)
        
        if(((student - mu1) %*% invSigma1 %*% t(student - mu1)) > cutoff)
        {
          vecteurs_mu2 = rbind(vecteurs_mu2,student)
          compt = compt + 1
        }
      }
      
    }
    else if (contamin == "UniformeTronquee") {
      vecteurs_mu1 <- mvrnorm(n1, mu1, Sigma1)
      # vecteurs_mu2 <- rmvt(n2, mu1, sigma = Sigma1, df = 1, ncores = 1, A = NULL)
      # vecteurs_mu2 <- ifelse(vecteurs_mu2 > -2 & vecteurs_mu2 < 2, 
      #                        ifelse(vecteurs_mu2 < 0, -2, 2), 
      #                        vecteurs_mu2)
      compt = 0
      vecteurs_mu2 = matrix(0,0, d)
      invSigma1 = solve(Sigma1)
      
      
      while (compt < n2) 
      {
        uniforme = runif(d,-2,2)
        
        if((t(uniforme - mu1) %*% invSigma1 %*% (uniforme - mu1)) > cutoff)
        {
          vecteurs_mu2 = rbind(vecteurs_mu2,uniforme)
          compt = compt + 1
        }
      }
      
    }
    else if (contamin == "uniforme") {
      vecteurs_mu1 <- mvrnorm(n1, mu1, Sigma1)
      vecteurs_mu2 <- matrix(runif(n2 * d, min = -10, max = 10), nrow = n2, ncol = d)
      #Contamination avec des données hors de l'ellipse X^T Sigma^(-1) X > q_{1 - a}
      # #Suppression des valeurs comprises entre -2 : remplacement par -2  pour les valeurs et 2
      # vecteurs_mu2 <- ifelse(vecteurs_mu2 > -2 & vecteurs_mu2 < 2, 
      #                        ifelse(vecteurs_mu2 < 0, -2, 2), 
      #                        vecteurs_mu2)
      # 
    } 
    else if (contamin == "zero") {
      vecteurs_mu1 <- mvrnorm(n1, mu1, Sigma1)
      vecteurs_mu2 <- matrix(0, nrow = n2, ncol = d)
    } 
    else if (contamin == "MaronnaZamar")
    {
      vecteurs_mu1 <- mvrnorm(n1, mu1, diag(d)) #Non outliers
      print(dim(vecteurs_mu1))
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
      #print(dim(R))
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
  if(cluster == FALSE){
    # Mélanger aléatoirement les données
    set.seed(123)  # Pour garantir la reproductibilité
    indices <- sample(nrow(Z))
    Z <- Z[indices, ]
    labelsVrais <- labelsVrais[indices]
  }
  else if (cluster == TRUE)
  {
    points_par_cluster <- d
    nclust = n2/d
    # Créer grappes d'outliers bien espacées
    Z <- matrix(NA, nrow = n, ncol = d)  # assez de place
    labelsVrais <- rep(NA, nrow(Z))
    
    # Placer les outliers à des indices spécifiques
    debut_positions <- seq(5, by = 50, length.out = nclust)  # positions 1, 51, 201, ...
    print(paste0("début positions",debut_positions))
    #print(length(debut_positions))
    compt = 1
    for (i in 1:(nclust)) {
      debut_indices = (compt - 1) * points_par_cluster + 1
      indices <- debut_indices : (debut_indices + (points_par_cluster - 1))
      # print(paste0("longueur indices ",length(indices)))
      # print(paste0("indices ", indices))
      # #print(debut_positions[i]:(debut_positions[i] + points_par_cluster))
      # print(paste0("positions ",debut_positions[i]:(debut_positions[i] + points_par_cluster - 1)))
      # print(paste0("longueur positions "), length(debut_positions[i]:(debut_positions[i] + points_par_cluster - 1)))
      Z[debut_positions[i]:(debut_positions[i] + points_par_cluster - 1), ] <- vecteurs_mu2[indices, ]
      labelsVrais[debut_positions[i]:(debut_positions[i] + points_par_cluster - 1)] <- 1
      compt = compt + 1
    }
    
    # Remplir les trous restants avec les points "propres"
    remaining_rows <- which(is.na(labelsVrais))
    print(length(remaining_rows))
    Z[remaining_rows, ] <- vecteurs_mu1
    labelsVrais[remaining_rows] <- 0    }
  # Calcul des vraies valeurs pour la médiane géométrique
  #Vvrai <- WeiszfeldCov(Z, nitermax = 1000000)$covmedian
  #VcovVrai <- GmedianCov(Z, scores = 10)
  #VpvraiesV <- eigen(VcovVrai$covmedian)$values
  
  # Retourner les résultats
  return(list(Z = Z,labelsVrais = labelsVrais))
}




