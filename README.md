# Algorithmes stochastiques pour la statistique robuste 

## Fichiers du code

* simulations.R : génère une matrice de variance de Toeplitz à partir d'un paramètre $\rho$. Génère $n$ vecteurs de dimension $d$, avec ou sans outliers selon 3 modes de contamination (en moyenne, en variance, avec une loi de Student à 1 degré de liberté).
      
* Parametres.R : renvoie les différents paramètres communs et spécifiques aux différentes méthodes, notamment : 
  - Y : matrice dont les lignes sont générées selon une loi donnée ;
  - k : taille des blocs (spécifique au streaming) ;
  - c : learning rate (spécifique à la méthode online et à la méthode offline) ;

* algorithmes.R : contient les différents algorithmes qui visent à
    - estimer la médiane géométrique, et la matrice de covariance ;
    - calculer les seuils de détection des outliers ;
    - détecter les outliers ;
 Fonctionne à partir des données générées (notamment le vecteur Y) par la fonction generer_echantillon du fichier simulations.R
      
* resultats.R : a trois fonctions principales :
    - calculer les erreurs ;
    - afficher les erreurs d'estimation des différentes quantités d'intérêt :
    - calculer les boucles de détection des outliers ;
 
* main.R : exécution des fonctions précédentes

## Exemple d'utilisation 

Exemple d'utilisation pour une estimation des différents paramètres et une détection des outliers : 
```r
source("~/algosto/parametres.R")
source("~/algosto/simulations.R")
source("~/algosto/algorithmes.R")
taux_contamination <- 2
p1 <- 1 - delta / 100
p2 <- 1 - p1
#Génération de l'échantillon à partir de la fonction genererEchantillon de parametres.R  
resultsSimul <- genererEchantillon(n, d, mu1, mu2, p1, p2, Sigma1 = Sigma1, Sigma2 = Sigma2,contamin = "variance")
Z <- resultsSimul$Z

#Calcul des quantités d'intérêt
results <- estimMVOutliers(Z, params$c, params$n, params$d, params$d, params$r, aa = 0.75, niter = 1e4,niterRMon = d ,methode = "eigen")
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


