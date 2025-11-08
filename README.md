# Online Robust Estimation of the Variance and Outlier Detection



Ce projet contient un pipeline de simulation en R destiné à estimer online de manière robuste la matrice de variance, et à détecter les outliers.

## Structure du projet

Le projet s'articule autour de quatre fichiers principaux :

### 0. `loadnecessary.R`

- **Rôle** : charge les packages et fonctions nécessaires et les packages nécessaires

### 1. `parametres.R`

- **Rôle** : Génère les fichiers de paramètres pour les différents scénarios de simulation.
- **Entrée** : Aucune.
- **Sortie** :
  - `SimParmsGrid-n10000-d10.Rdata` : Grille de paramètres pour $d = 10$.
  - Plusieurs fichiers nommés selon le format :
    ```
    Parms-d<d>-n<n>-k<k>-l<l>-rho1<rho1>-r<r>-sim<sim>.Rdata
    ```
    où chaque variable représente un paramètre de simulation.


### 2. `algorithmes.R`

- **Rôle** : Lance les méthodes d'estimation et de détection d'outlier concurrentes : cov online naïf, shrinkage selon les idées de Wolf et de Cabana.
- **Entrée** :
  - Dataset $Z \in \mathcal{M}_{n,d}(\mathbb{R})$.
- **Sortie** :
  - Estimateur de $\Sigma$ : $\widehat{\Sigma}$ final et à chaque itération pour méthodes online, outliers_labels (vecteurs où l'entrée $i$ si la donnée $i$ traitée est un outlier, $0$ sinon.

  
### 3. `affiche_resultats.R`

- **Rôle** : Graphiques (scénarios de contamination, erreurs norme de Frobenius matrice de covariance, faux positifs, faux négatifs)
- **Entrée** :
  - Résultats des runs 
- **Sortie** :
  - Graphiques.
  
