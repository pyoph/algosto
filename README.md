---
editor_options: 
  markdown: 
    wrap: 72
---

# Algorithmes stochastiques pour la statistique robuste

## Fichiers du code

-   parametres.R : initialisation des différents paramètres communs à
    toutes les simulations. Fonctions de création de matrices de
    Toeplitz, et de permutation de lignes et de colonnes pour les
    matrices de variance et de covariance.

-   simulations.R : Génère $n$ vecteurs de dimension $d$, avec ou sans
    outliers selon 3 modes de contamination (en moyenne, en variance,
    avec une loi de Student à 1 degré de liberté). Fonction principale :
    genererEchantillon() : utilisée pour générer un échantillon non
    contaminé ou contaminé

-   algorithmes.R : contient les différents algorithmes qui visent à

    -   estimer la médiane géométrique, et la matrice de covariance ;
        Fonctionne à partir des données générées (notamment le vecteur
        Y) par la fonction generer_echantillon du fichier simulations.R
        2 fonctions principales :
        
          estimV (estimation en ligne) et streaming (pour une estimation par blocs)

-   resultats.R : a comme fonctions principales :

    -   calculer les erreurs ;
    -   afficher les erreurs d'estimation des différentes quantités
        d'intérêt

-   Outliers.R : 
    - fonction qui calcule les distances de Mahalanobis pour différentes méthodes ;
    - détecte les outliers en focntion d'une distance et d'un seuil

-   seuils.R : 
    - fonction courbeROC qui calcule les AUC pour différentes méthodes ;

-   computeOutliers.R : - applique plusieurs algorithmes basés sur
    différentes estimations de la matrice de covariance pour détecter
    les outliers

-   main.R :

    -   lancement des fonctions dans les précédents fichiers

## Exemple d'utilisation

Dans le fichier parametres.R modifier les paramètres souhaités \`\`\`r

nbruns = 1 n = 1e4 d = 10 mu1 = rep(0,d) mu2 = 5\*rep(1,d) rho = 0.8
Sigma1 = creerMatriceToeplitz(rho,d)

lignes_a_permuter = c(1, 2) colonnes_a_permuter = c(1, 2) Sigma2 =
permuterLignesColonnes(Sigma1,lignes_a_permuter , colonnes_a_permuter)

Dans le fichier main.R, lancer l'intégralité des lignes de code du
fichier pour la détection des outliers

\`\`\`r

library(reshape2) library(RobRegression) library(Gmedian)
library(ggplot2) library(far) library(gridExtra) library(microbenchmark)
library("matlib") library("MASS") library("corrplot") library("dplyr")
source("\~/algosto/parametres.R") source("\~/algosto/simulations.R")
source("\~/algosto/algorithmes.R") source("\~/algosto/resultats.R")
source("\~/algosto/Outliers.R") source("\~/algosto/computeOutliers.R")


#Estimation de la matrice de covariance pour différents taux de contamination



#Utilisation de la fonction d'estimation online et affichage des erreurs


taux_contamination <- c(0,2, 5, 10, 15, 20, 25, 30, 40)



for (i in seq_along(taux_contamination)) 
{
  delta <- taux_contamination[i]
 
 #Proportion dans l'échantillon non contaminée et contaminée
  p1 <- 1 - delta / 100
  
  p2 <- 1 - p1
  #mu1
 
  
  resultsSimul <- genererEchantillon(n, d, mu1, mu2, p1, p2, Sigma1 = Sigma1, Sigma2 = Sigma2,contamin = "moyenne")
  Z <- resultsSimul$Z
  
  #Estimation online
  results <- estimMV(Z)
  
  
  #Calcul des erreurs
  
  erreursMiter <- calculErreursM(results$miter,rep(0,d))
  erreursSigma <- calculErreursNormeFrobenius(results$Sigma,Sigma1)
  
  #Affichage des erreurs
  
  #affiche_erreursM(erreursMiter,contamination = delta)
  
  affiche_erreursSigma(erreursSigma,contamination= delta)
  
  
}



#Calcul des outliers

for (i in (1 : nbruns)){ results_outliers \<-
calcule_outliers(n,d,c,rho,mu1,mu2,Sigma1,Sigma2,contamin = "student") }


