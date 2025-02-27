#Shrinkage pour la médiane géométrique 

shrinkage_med <- function(Y,muVrai = mu1)
{
  
  #Calcul de la médiane composante par composante
  muCCm <- colMedians(Y)
  
  #Calcul de target
  
  nu <- (muCCm %*% rep(1,ncol(Y)))/ncol(Y)
  
  target <- rep(nu,ncol(Y))
  
  VarTot <- 0
  
  for (j in 1:ncol(Y))
  {
  VarTot <- VarTot + var(Y[,j])  
  }
  
  #Calcul de eta 
  
  eta <- pi/(2*ncol(Y))*VarTot
  
  eta <- eta/(eta + sum(diag(t(muVrai - target)%*%(muVrai - target))))
  
  muShrink <- eta*target + (1 - eta)*muCCm
  
  return (list(muCCm = muCCm, nu = nu,eta = eta,muShrink = muShrink))
}


#Reprise fonction package Michael Wolf changement de matrice de covariance empirique par SCCM

shrinkage_SCCM <- function(Y, k = -1) {
  dim.Y <- dim(Y)
  N <- dim.Y[1]  
  p <- dim.Y[2]  
  
  if (k < 0) {    
    Y <- scale(Y, scale = FALSE)  # Centrage des données
    k <- 1
  }
  n <- N - k    # Taille d'échantillon effective
  c <- p / n    # Ratio de concentration
  
  #  Calcul de la Matrice de Comédiane (COM)
  COM <- matrix(0, p, p)
  for (j in 1:p) {
    for (t in 1:p) {
      COM[j, t] <- median((Y[, j] - median(Y[,j])) * (Y[, t]  - median(Y[,t]))) 
    }
  }
  
  # Ajustement SCCM
  SCCM <- 2.198 * COM
  
  #  Calcul de la cible de shrinkage
  SCCM_diag <- diag(SCCM)
  sqrtvar <- sqrt(SCCM_diag)
  rBar <- (sum(SCCM / outer(sqrtvar, sqrtvar)) - p) / (p * (p - 1))
  #target <- rBar * outer(sqrtvar, sqrtvar)
  #diag(target) <- SCCM_diag  # Conservation de la diagonale
  target <- diag(ncol(SCCM))
  nu <- sum(diag(SCCM))/ncol(SCCM) 
  
  target <- nu*target
  
  # Estimation de pi (variance des éléments hors-diagonale)

  # SCCM2 <- 2.198 * COM2
  piMat <- SCCM - SCCM^2
  pihat <- sum(piMat)
  
  #  Estimation de gamma (distance entre SCCM et la cible)
  gammahat <- norm((SCCM - target), type = "F")^2
  
  # ⃣ Partie diagonale de rho 
  rho_diag <- sum(diag(piMat))
  
  #  Partie hors-diagonale de rho 
  thetaMat <- matrix(0, p, p)
  for (j in 1:p) {
    for (t in 1:p) {
      if (j != t) {
        thetaMat[j, t] <- median(Y[, j] * Y[, t]) - SCCM_diag[j] * SCCM_diag[t]
      }
    }
  }
  rho_off <- rBar * sum(outer(1/sqrtvar, sqrtvar) * thetaMat)
  
  #  Calcul de l'intensité du shrinkage
  rhohat <- rho_diag + rho_off
  kappahat <- (pihat - rhohat) / gammahat
  shrinkage <- max(0, min(1, kappahat / n))
  
  #  Calcul de l'estimateur shrinké
  sigmahat <- shrinkage * target + (1 - shrinkage) * SCCM
  
  return(list(SCCM_shrinked = sigmahat,SCCM = SCCM ,shrinkage_intensity = shrinkage))
}

