# Algorithmes stochastiques pour la statistique robuste 

## Fichiers du code

* simulations.R : génère une matrice de variance de Toeplitz à partir d'un paramètre $\rho$. Génère $n$ vecteurs de dimension $d$, avec ou sans outliers selon 3 modes de contamination (en moyenne, en variance, avec une loi de Student à 1 degré de liberté).
      

* algorithmes.R : contient les différents algorithmes qui visent à
    - estimer la médiane géométrique, et la matrice de covariance ;
     Fonctionne à partir des données générées (notamment le vecteur Y) par la fonction generer_echantillon du fichier simulations.R

* resultats.R : a trois fonctions principales :
    - calculer les erreurs ;
    - afficher les erreurs d'estimation des différentes quantités d'intérêt :

* Outliers.R : fonctions qui : 
      - calculent les distances de Mahalanobis ;
      - detectent les outliers

* computeOutliers.R : 
      - appliquent plusieurs algorithmes pour détecter les outliers

* main.R :
  - lancement des fonctions dans les précédents fichiers

## Exemple d'utilisation 

Exemple d'utilisation pour une estimation des différents paramètres et une détection des outliers : 
```r
source("~/algosto/Outliers.R")
source("~/algosto/simulations.R")
source("~/algosto/algorithmes.R")
taux_contamination <- 2
p1 <- 1 - delta / 100
p2 <- 1 - p1
#Génération de l'échantillon à partir de la fonction genererEchantillon de parametres.R  
resultsSimul <- genererEchantillon(n, d, mu1, mu2, p1, p2, Sigma1 = Sigma1, Sigma2 = Sigma2,contamin = "variance")
Z <- resultsSimul$Z

#Calcul des quantités d'intérêt
results <- estimMVOutliers(Z, c = sqrt(d), n = nrow(Z), d = ncol(Z), d = ncol(Z), r = 1.5, aa = 0.75, niter = 1e4,niterRMon = d ,methode = "eigen")
#Médiane géométrique
miter <- results$miter
#Vecteurs propres
U <- results$U
#valeurs propres
lambda <- results$lambdaIter
#Distances de Mahalanobis
distances <- calcule_vecteur_distancesOnline(Z,miter,U,lambda)
#Recherche des outliers
outliers_labels <- detectionOutliers(distances, cutoff =  qchisq(p = 0.95, df = ncol(Z)))
#Tableau de contingence
tc <- table_contingence(resultsSimul$labelsVrais[1:9999], outliers_labels[1:9999])


