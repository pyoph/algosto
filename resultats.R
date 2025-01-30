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
library(ROCR)
#install.packages("ROCR")

#Calcule les erreurs d'estimation de la médiane géométrique critère RMSE

calculErreursM <- function(miter,mvrai)
{
  nb_iterations <- dim(miter)[2]
  erreursM <- rep(0,nb_iterations)
  
  for (i in 1:(nb_iterations-1))
  {
    erreursM[i] <- sqrt(sum((mvrai - miter[i,])^2))
    
  }
  return (erreursM)
  
}

#Calcul des erreurs en norme de Frobenius

calculErreursNormeFrobenius <- function (Estim,Vrai)
{
  nbiterations <- dim(Estim)[1]
  
  erreurs <- rep(0,nbiterations)
  
  for (i in (1:nbiterations))
  {
    erreurs[i] <-  norm(Vrai - Estim[i,,],type = "F")
  }
  return (erreurs)
}



# Fonction pour afficher les erreurs d'estimation de m online et offline pour différents taux de contamination. 

#Les erreurs offline seront par la suite ajoutées
affiche_erreursM <- function(erreurs_online, erreurs_offline = NULL, contamination = 2) {
  
  nberreurs_online <- length(erreurs_online)
  
  if (!is.null(erreurs_offline)) {
    nberreurs_offline <- length(erreurs_offline)
    
    if (nberreurs_online != nberreurs_offline) {
      stop("Les vecteurs 'erreurs_online' et 'erreurs_offline' doivent avoir la même longueur.")
    }
    
    # Créer le data frame combiné
    data <- data.frame(
      x = rep(1:(nberreurs_online - 1), 2),
      erreurs = c(erreurs_online[1:(nberreurs_online - 1)], erreurs_offline[1:(nberreurs_offline - 1)]),
      type = c(rep("Online", nberreurs_online - 1), rep("Offline", nberreurs_offline - 1))
    )
    
    # Définir les couleurs et labels
    colors <- c("Online" = "blue", "Offline" = "green")
    labels <- c(
      paste("Online (contamination = ", contamination , "%)", sep = ""),
      paste("Offline (contamination = ", contamination , "%)", sep = "")
    )
  } else {
    # Créer le data frame uniquement pour les erreurs online
    data <- data.frame(
      x = 1:(nberreurs_online - 1),
      erreurs = erreurs_online[1:(nberreurs_online - 1)],
      type = rep("Online", nberreurs_online - 1)
    )
    
    # Définir les couleurs et labels pour online uniquement
    colors <- c("Online" = "blue")
    labels <- c(paste("Online (contamination = ", contamination , "%)", sep = ""))
  }
  
  # Créer le graphique
  p <- ggplot(data, aes(x = x, y = erreurs, color = type, linetype = type)) +
    geom_line(size = 1) +
    scale_x_log10() +   # Échelle logarithmique pour l'axe x
    scale_y_log10() +   # Échelle logarithmique pour l'axe y
    scale_color_manual(values = colors, labels = labels) +
    labs(x = "Taille de l'échantillon", 
         y = "Erreur RMSE",
         title = "Comparaison des erreurs d'estimation de m",
         color = "Méthode", 
         linetype = "Méthode") +
    theme_minimal()
  
  print(p)
}


# Fonction pour afficher les erreurs d'estimation de Sigma online et offline pour différents taux de contamination. 

#Les erreurs offline seront par la suite ajoutées
affiche_erreursSigma <- function(erreurs_online, erreurs_offline = NULL, contamination = 0) {
  
  nberreurs_online <- length(erreurs_online)
  
  if (!is.null(erreurs_offline)) {
    nberreurs_offline <- length(erreurs_offline)
    
    if (nberreurs_online != nberreurs_offline) {
      stop("Les vecteurs 'erreurs_online' et 'erreurs_offline' doivent avoir la même longueur.")
    }
    
    # Créer le data frame combiné
    data <- data.frame(
      x = rep(1:(nberreurs_online - 1), 2),
      erreurs = c(erreurs_online[1:(nberreurs_online - 1)], erreurs_offline[1:(nberreurs_offline - 1)]),
      type = c(rep("Online", nberreurs_online - 1), rep("Offline", nberreurs_offline - 1))
    )
    
    # Définir les couleurs et labels
    colors <- c("Online" = "blue", "Offline" = "green")
    labels <- c(
      paste("Online (contamination = ", contamination , "%)", sep = ""),
      paste("Offline (contamination = ", contamination , "%)", sep = "")
    )
  } else {
    # Créer le data frame uniquement pour les erreurs online
    data <- data.frame(
      x = 1:(nberreurs_online - 1),
      erreurs = erreurs_online[1:(nberreurs_online - 1)],
      type = rep("Online", nberreurs_online - 1)
    )
    
    # Définir les couleurs et labels pour online uniquement
    colors <- c("Online" = "blue")
    labels <- c(paste("Online (contamination = ", contamination , "%)", sep = ""))
  }
  
  # Créer le graphique
  p <- ggplot(data, aes(x = x, y = erreurs, color = type)) +
    geom_line(size = 1) +
    scale_x_log10() +   # Échelle logarithmique pour l'axe x
    scale_y_log10() +   # Échelle logarithmique pour l'axe y
    scale_color_manual(values = colors, labels = labels) +
    labs(x = "Taille de l'échantillon", 
         y = "Erreur norme de Frobenius",
         title = "Comparaison des erreurs d'estimation de Sigma",
         color = "Méthode") +   # Une seule rubrique pour "Méthode"
    theme_minimal()
  
  print(p)
}


#Affichage comparaison Sigma

plot_comparaison_sigma <- function(Sigma1, SigmaEstimOnline, SigmaEstimOffline, delta) {
  par(mfrow = c(1, 2))  # Afficher deux graphiques côte à côte
  
  # Comparaison Online
  plot(Sigma1, SigmaEstimOnline, col = 2, 
       main = paste("Comparaison Online (Contamination:", delta, "%)"), 
       xlab = "Sigma réelle", ylab = "Estimations Online")
  abline(0, 1)
  points(Sigma1, SigmaEstimOnline, col=2)
  
  # Comparaison Offline
  plot(Sigma1, SigmaEstimOffline, col = 4, 
       main = paste("Comparaison Offline (Contamination:", delta, "%)"), 
       xlab = "Sigma réelle", ylab = "Estimation Offline")
  abline(0, 1)
  
  points(Sigma1, SigmaEstimOffline, col=2)
  par(mfrow = c(1, 1))  # Réinitialisation du paramètre graphique
}


#Affichage boxplot des erreurs

creer_boxplot_erreurs <- function(erreursSigmaBoxplot, taux_contamination, methode) {
  # Vérifier que la matrice d'erreurs a bien des lignes correspondant aux taux de contamination
  if (length(taux_contamination) != nrow(erreursSigmaBoxplot)) {
    stop("Le nombre de lignes de erreursSigmaBoxplot doit correspondre à la longueur de taux_contamination")
  }
  
  # Transformer la matrice en un format long pour ggplot
  library(ggplot2)
  library(reshape2)
  
  # Conversion de la matrice en data frame long
  df_erreurs <- data.frame(erreursSigmaBoxplot)
  df_erreurs$taux_contamination <- factor(taux_contamination)  # Facteur pour l'axe X
  df_long <- melt(df_erreurs, id.vars = "taux_contamination")
  
  # Création du boxplot avec ggplot2
  p <- ggplot(df_long, aes(x = taux_contamination, y = value)) +
    geom_boxplot(fill = "lightblue", color = "blue", outlier.color = "red") +
    labs(
      title = paste("Boxplots des erreurs - Méthode :", methode),
      x = "Taux de contamination (%)",
      y = "Erreur (norme de Frobenius)"
    ) +
    theme_minimal()
  
  print(p)
}