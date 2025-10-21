



#Estimation de la matrice de covariance par shrinkage

covCor <- function(Y, k = -1) {
  dim.Y <- dim(Y)
  N <- dim.Y[1]
  p <- dim.Y[2]
  if (k < 0) {    # demean the data and set k = 1
    Y <- scale(Y, scale = F)
    k <- 1
  }
  n <- N - k    # effective sample size
  c <- p / n    # concentration ratio
  #sample <- (t(Y) %*% Y) / n   
  sample <- covComed(Y)$cov
  # compute shrinkage target
  samplevar <- diag(sample)
  sqrtvar <- sqrt(samplevar)
  rBar <- (sum(sample / outer(sqrtvar, sqrtvar)) - p) / (p * (p - 1))
  target <- rBar * outer(sqrtvar, sqrtvar)
  diag(target) <- samplevar
  
  # estimate the parameter that we call pi in Ledoit and Wolf (2003, JEF)
  Y2 <- Y^2
  sample2 <- (t(Y2) %*% Y2) / n   
  #sample2 <- covComed(Y2)
  piMat <- sample2 - sample^2
  pihat <- sum(piMat)
  
  # estimate the parameter that we call gamma in Ledoit and Wolf (2003, JEF)
  gammahat <- norm(c(sample - target), type = "2")^2
  
  # diagonal part of the parameter that we call rho 
  rho_diag <- sum(diag(piMat))
  
  # off-diagonal part of the parameter that we call rho 
  term1 <- (t(Y^3) %*% Y) / n;
  term2 <- rep.row(samplevar, p) * sample;
  term2 <- t(term2)
  thetaMat <- term1 - term2
  diag(thetaMat) <- 0
  rho_off <- rBar * sum(outer(1/sqrtvar, sqrtvar) * thetaMat)
  
  # compute shrinkage intensity
  rhohat <- rho_diag + rho_off
  kappahat <- (pihat - rhohat) / gammahat
  shrinkage <- max(0, min(1, kappahat / n))
  
  # compute shrinkage estimator
  sigmahat <- shrinkage * target + (1 - shrinkage) * sample
  
  
  return (sigmahat)
  
}

rep.row <- function(x, n){
  matrix(rep(x, each = n), nrow = n)
}


##############################################
####################Online sample covariance
###############################################SA

SampleCovOnline = function(Z)
{
  nblignes = nrow(Z)
  
  # Initialisation avec les 2 premières données
  mean = colMeans(Z[1:2, , drop = FALSE])        # moyenne empirique
  meanOld = mean
  Sigma = cov(Z[1:2, , drop = FALSE])            # covariance empirique
  
  meanIter = matrix(0, nrow(Z), ncol(Z))
  SigmaIter = array(0, dim = c(nrow(Z), ncol(Z), ncol(Z)))
  
  
  meanIter[1,] = mean
  
  inSigma = diag(ncol(Z))  
  
  

  distances = rep(0, nrow(Z))
  outliers_labels = rep(0, nrow(Z))
  cutoff = qchisq(.95, df = ncol(Z))
  
  meanOld = mean
  #Sigma = 1.5*diag(ncol(Z))
  
  #SigmaIter = array(0, dim = c(nrow(Z),ncol(Z),ncol(Z)))
  distances = rep(0, nrow(Z))
  outliers_labels = rep(0,nrow(Z))
  cutoff = qchisq(.95,df = ncol(Z))
  #nb_out= 0
  for (i in (2:(nblignes-1)))
  {

    mean   = mean   + (1.0 / (i + 1)) * (Z[i+1,] - mean);
    Sigma = (i - 1)/i*Sigma + 1/(i + 1)*((Z[i+1,] - meanOld)%*%t(Z[i+1,] - meanOld))
    meanOld = mean
    meanIter[i,] = mean
    SigmaIter[i,,] = Sigma
    #distances[i] = mahalanobis_generalizedRcpp(Z[i,],meanIter[i,],eigen(SigmaIter[i,,])$vectors, eigen(SigmaIter[i,,])$values)
    #distances[i] = mahalanobis(Z[i,],meanIter[i,],SigmaIter[i,,],inverted = FALSE)
    #distances[i] = t(Z[i,] - meanIter[i,])%*%solve(SigmaIter[i,,])%*%(Z[i,] - meanIter[i,])
    
    scal = 1 + t(Z[i+1,])%*%invSigma%*%Z[i+1,]
    if(scal !=0){
      invA = invA - 1/scal[1]*invSigma%*%(Z[i+1,]%*%t(Z[i+1,]))%*%invA}
    distances[i] = t(Z[i,] - meanIter[i,])%*%invA%*%(Z[i,] - meanIter[i,])
    
    S = distances[i]
    
    if (S > cutoff) {outliers_labels[i] = 1
    }  
  }
    
    
  
    # if(!is.nan(S)){
    # 
    # if (S > cutoff) {outliers_labels[i] = 1}}else{print("S = Nan ")}      
    # }
    # 
    #print(paste0("cutoff =",cutoff))
    
  SigmaIter[nrow(Z),,] = Sigma
  meanIter[nrow(Z),] = mean
  return(list(mean = mean, Sigma = Sigma, meanIter = meanIter, SigmaIter = SigmaIter,distances = distances, outliers_labels = outliers_labels))
  
  }

  


#Exclusion des valeurs propres trop grandes
reduce_dimension = function(Sigma){
  
  #Exclusion des valeurs aberrantes pour l'estimation des valeurs propres 
  
  distances = rep(0, nrow(Z))
  outliers_labels = rep(0,nrow(Z))
  cutoffCorr = rep(0,nrow(Z))
  
  for (i in (1:nrow(Z))){
    
    
    lambda = eigen(Sigma[i,,])$values
    Q <- quantile(lambda, probs = c(0.1, 0.9))
    IQR <- Q[2] - Q[1]
    limites <- c(Q[1] - 1.5 * IQR, Q[2] + 1.5 * IQR)
    indices_non_aberrants <- which(lambda >= max(limites[1],0) & lambda <= limites[2])
    lambda_filtre <- lambda[lambda >= max(limites[1],0) & lambda <= limites[2]]
    #print(paste0("indices non aberrants ", indices_non_aberrants))
    #print(paste0("lambda_filtre ",lambda_filtre))
    distances[i] = mahalanobis_generalizedRcpp(Z[i,indices_non_aberrants],resultats$miter[i,indices_non_aberrants],eigen(Sigma[i,,])$vectors[indices_non_aberrants,indices_non_aberrants], eigen(Sigma[i,,])$values[indices_non_aberrants])
    S = distances[i]
    
    
    
    cutoffCorr[i]  = qchisq(.95,df = d)*median(resultats$distances[1:i])/qchisq(.5,df = d)
    if (distances[i] > cutoffCorr[i]) {outliers_labels[i] = 1}
  }
  return(outliers_labels)
}
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
