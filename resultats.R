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
source("C://Users/Paul/Downloads/parametres.R")
source("C://Users/Paul/Downloads/simulations.R")
source("C://Users/Paul/Downloads/algorithmes.R")

#Calcule les erreurs d'estimation de la médiane géométrique
calculErreursM <- function(miter,mvrai,n)
{

  erreursM <- rep(0,n)

  for (i in 1:(n-1))
  {
    erreursM[i] <- sqrt(sum((mvrai - miter[i,])^2))
    #print(erreursM[i])
  }
  return (erreursM)

}

#Calcule les erreurs d'estimation de la matrice de covariance médiane comparer avec I_d

calculErreursV <- function(Viter,Vvrai,n)
{
  erreursV <- rep(0,n-1)

  for (i in 1:(n-1))
  {
    erreursV[i] <- norm(Vvrai - Viter[,,i],type = "F")
  }
return(erreursV)
}

#Calcule les erreurs d'estimation des valeurs propres de la matrice de covariance médiane

calculErreursVpMCM <- function(vpMCM,vpMCMvrai,n)
{

  erreursVpMCM <- rep(0,n - 1)

  for (i in 1:(n-1))
  {
    erreursVpMCM[i] <- sqrt(sum((vpMCMvrai - vpMCM[i,])^2))
  }
  return(erreursVpMCM)
}

#Calcul norme de Frobenius matrice de covariance


#Calcule les erreurs d'estimation des valeurs propres de la matrice de covariance médiane

calculErreursVpMCM <- function(vpMCM,vpMCMvrai,n)
{

  erreursVpMCM <- rep(0,n - 1)

  for (i in 1:(n-1))
  {
    erreursVpMCM[i] <- sqrt(sum((vpMCMvrai - vpMCM[i,])^2))
  }
  return(erreursVpMCM)
}



#Calcule les erreurs d'estimation des vecteurs propres de la MCM

calculErreursVectpMCM <- function(UMCM,UMCMvrai,n,d)
{

  erreursUMCM <- rep(0,n - 1)

  normes_colonnes <- apply(UMCMvrai, 2, function(col) sqrt(sum(col^2)))

  UMCMvrai <- 1/normes_colonnes*UMCMvrai


  for (i in 1:(n-1))
  {
    normes_colonnes <- apply(UMCM[i,,], 2, function(col) sqrt(sum(col^2)))

    UMCM[i,,] <- 1/normes_colonnes*abs(UMCM[i,,])

    erreursUMCM[i] <- norm(UMCMvrai - UMCM[i,,],type = "F")
  }
  return(erreursUMCM)
}


#Calcule les erreurs d'estimation de la matrice Sigma

calculErreursSigma <- function(SigmaVrai,U,lambda,n,d)
{

  erreursSigmaEstim <- rep(0,n - 1)


  for (i in 1:(n-1))
  {

    SigmaEstim <- solve(U[i,,]) %*% diag(lambda) %*% U[i,,]

    erreursSigmaEstim[i] <- norm(SigmaEstim - SigmaVrai,type = "F")
  }
  return(erreursSigmaEstim)
}


#Calcule les erreurs d'estimation des valeurs propres de la matrice de covariance

calculErreursValPropreCov <- function(lambdaIter,lambdaVrai,n)
{

  erreursLambda <- rep(0,n - 1)

  for (i in 1:(n-1))
  {
    erreursLambda[i] <- sqrt(sum((lambdaIter[i,] - lambdaVrai)^2))
  }
  return(erreursLambda)
}



#Affiche les erreurs d'estimation de la médiane géométrique
affiche_erreurs_mediane_geometrique <- function(erreursM,n)
{
  # Créer le data frame
  data <- data.frame(
    x = 1:(n - 1),
    erreursM = erreursM[1:(n - 1)],
    type = rep("Erreur d'estimation de m", n - 1)
  )


  # Créer le graphique avec ggplot2
  p <- ggplot(data, aes(x = x, y = erreursM, color = type)) +
    geom_line(size = 1) +
    scale_x_log10() +   # Échelle logarithmique pour l'axe x
    scale_y_log10() +   # Échelle logarithmique pour l'axe y
    scale_color_manual(values = c("Erreur d'estimation de m" = "red"),
                       labels = c("Erreur d'estimation de m")) +
    labs(x = "Taille de l'échantillon", y = "Erreur Quadratique",
         title = "Erreur d'estimation de m",
         color = "Type d'Erreur") +
    theme_minimal()


  # Afficher le graphique
  print(p)

}


#Affiche les erreurs d'estimation de la médiane géométrique
affiche_erreurs_mcm <- function(erreursV,n,c)
{

  # Créer le data frame
  data <- data.frame(
    x = 1:(n - 1),
    erreursV = erreursV[1:(n - 1)],
    type = rep("Erreur d'estimation de V", n - 1)
  )


  p <- ggplot(data, aes(x = x, y = erreursV, color = type)) +
    geom_line(size = 1) +
    scale_x_log10() +   # Échelle logarithmique pour l'axe x
    scale_y_log10() +   # Échelle logarithmique pour l'axe y
    scale_color_manual(values = c("Erreur V" = "red"),
                       labels = c(paste("Erreur d'estimation de V c = ", c))) +
    labs(x = "Taille de l'échantillon", y = "Erreur norme de Frobenius",
         title = "Erreur d'estimation de V",
         color = "Type d'Erreur") +
    theme_minimal()

  print(p)


}

#Affiche les erreurs d'estimation des valeurs propres de la médiane géométrique
affiche_erreurs_vp_mcm <- function(erreursVpMCM,n)
{

  # Créer le data frame
  data <- data.frame(
    x = 1:(n - 1),
    erreursVpMCM = erreursVpMCM[1:(n - 1)],
    type = rep("Erreur d'estimation de m", n - 1)
  )


  # Créer le graphique avec ggplot2
  p <- ggplot(data, aes(x = x, y = erreursVpMCM, color = type)) +
    geom_line(size = 1) +
    scale_x_log10() +   # Échelle logarithmique pour l'axe x
    scale_y_log10() +   # Échelle logarithmique pour l'axe y
    scale_color_manual(values = c("Erreur d'estimation de m" = "red"),
                       labels = c("Erreur d'estimation de m")) +
    labs(x = "Taille de l'échantillon", y = "Erreur Quadratique",
         title = "Erreur d'estimation des valeurs propres de la MCM",
         color = "Type d'Erreur") +
    theme_minimal()


  # Afficher le graphique
  print(p)

}


#Affiche les erreurs d'estimation des vecteurs propres de la MCM
affiche_erreurs_vectp_mcm <- function(erreursVectPMCM,n,c)
{
  # Créer le data frame
  data <- data.frame(
    x = 1:(n - 1),
    erreursVectPMCM = erreursVectPMCM[1:(n - 1)],
    type = rep("Erreur d'estimation des valeurs propres de la MCM", n - 1)
  )


  p <- ggplot(data, aes(x = x, y = erreursVectPMCM, color = type)) +
    geom_line(size = 1) +
    scale_x_log10() +   # Échelle logarithmique pour l'axe x
    scale_y_log10() +   # Échelle logarithmique pour l'axe y
    scale_color_manual(values = c("Erreur V" = "red"),
                       labels = c(paste("Erreur d'estimation des vecteurs propres de la MCM  c = ", c))) +
    labs(x = "Taille de l'échantillon", y = "Erreur norme de Frobenius",
         title = "Erreur d'estimation des vecteurs propres de la MCM",
         color = "Type d'Erreur") +
    theme_minimal()

  print(p)

}

#Affiche les erreurs d'estimation des valeurs propres de la matrice de covariance
affiche_erreurs_cov <- function(erreursCov,n)
{
  # Créer le data frame
  data <- data.frame(
    x = 1:(n - 1),
    erreursCov = erreursCov[1:(n - 1)],
    type = rep("Erreur d'estimation des valeurs propres de la matrice de covariance", n - 1)
  )


  # Créer le graphique avec ggplot2
  p <- ggplot(data, aes(x = x, y = erreursCov, color = type)) +
    geom_line(size = 1) +
    scale_x_log10() +   # Échelle logarithmique pour l'axe x
    scale_y_log10() +   # Échelle logarithmique pour l'axe y
    scale_color_manual(values = c("Erreur d'estimation des valeurs propres de la matrice de covariance" = "red"),
                       labels = c("Erreur d'estimation des valeurs propres de la matrice de covariance")) +
    labs(x = "Taille de l'échantillon", y = "Erreur Quadratique",
         title = "Erreur d'estimation des valeurs propres de la matrice de covariance",
         color = "Type d'Erreur") +
    theme_minimal()


  # Afficher le graphique
  print(p)

}


#Affiche les erreurs d'estimation de la médiane géométrique
affiche_erreurs_Sigma <- function(erreursSigma,n,c)
{

  # Créer le data frame
  data <- data.frame(
    x = 1:(n - 1),
    erreursSigma = erreursSigma[1:(n - 1)],
    type = rep("Erreur d'estimation de Sigma", n - 1)
  )


  p <- ggplot(data, aes(x = x, y = erreursSigma, color = type)) +
    geom_line(size = 1) +
    scale_x_log10() +   # Échelle logarithmique pour l'axe x
    scale_y_log10() +   # Échelle logarithmique pour l'axe y
    scale_color_manual(values = c("Erreur Sigma" = "red"),
                       labels = c(paste("Erreur d'estimation de Sigma c = ", c))) +
    labs(x = "Taille de l'échantillon", y = "Erreur norme de Frobenius",
         title = "Erreur d'estimation de Sigma",
         color = "Type d'Erreur") +
    theme_minimal()

  print(p)


}



#Affiche table de contingence

table_contingence <- function(labels,outliers_labels)
{
  table(labels, outliers_labels)
}


#Affiche l'histogramme et la densité du Khi deux

densiteHistKhi2 <- function(S,d)
{

  # Créer un dataframe pour ggplot
  df <- data.frame(S = S)

  # Densité théorique du Chi2
  x_vals <- seq(0, max(S), length.out = 500)
  chi2_density <- data.frame(
    x = x_vals,
    density = dchisq(x_vals, df = d)
  )

  # Représentation graphique : Histogramme des valeurs de S et densité théorique du Chi2
  ggplot(df, aes(x = S)) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
    geom_line(data = chi2_density, aes(x = x, y = density), color = "red", size = 1) +
    labs(title = "Comparaison des valeurs de S et de la densité théorique du Chi2",
         x = "Valeurs de S", y = "Densité") +
    theme_minimal()


}

