#' @title Generation of Gaussian mixture distribution with Student contamination
#'
#' @description Function to generate a sample of a gaussian mixture model with Student contamination
#' \eqn{X\sim \sum_{k=1}^K \pi_k Y_k}, where \eqn{Y_k\sim \mathcal N (\mu_k, \Sigma_k)}.
#'
#' @param nk Vector of the number of observations in each component, giving the proportions \eqn{\pi_k}
#' @param muG Matrix of the means of the Gaussian components \eqn{\mu = (\mu_1, \ldots, \mu_K)} TODO plutot array
#' @param sigmaG Array of the variance-covariance matrices of the Gaussian components
#' @param delta Proportion of contamination
#' @param dfT Degrees of freedom of the Student contamination
#' @param muT Matrix of the means of the Student contamination
#' @param sigmaT Array of the variances of the Student contamination
#'
#' @return A list containing the sample X, the indicator of contamination C and the latent variable Z, giving the classification
#'
#' @examples
#' \donttest{
#' # Dimensions
#' K <- 3; nk <- rep(500, K); d <- 5
#'
#' # Gaussian parameters
#' mu <- rbind(mu1=rep(0, d), mu2=rep(3, d), mu3=rep(-3, d))
#' sigma <- array(dim=c(K, d, d))
#' sigma[1, , ] <- 2*diag(d); sigma[2, , ] <- diag(1:d); sigma[3, , ] <- diag(1/(1:d))
#'
#' # Student contamination
#' fraction <- 0.3
#' dfT <- 1; dfT <- 1; muT <- mu; sigmaT <- sigma
#'
#' # Simulation
#' GMMcontStudent <- rGMMcontStudent(nk=nk, muG=mu, sigmaG=sigma,
#'                                   delta=fraction,
#'                                   dfT=dfT, muT=muT, sigmaT=sigmaT)
#' }
#'
#' @export
rGMMcontStudent <- function(nk, muG, sigmaG, delta, dfT, muT, sigmaT){
  K <- length(nk); p <- ncol(muG); n <- sum(nk)
  nkCum <- c(0, cumsum(nk))
  nkG <- round(nk*(1-delta)); nkA <- nk - nkG
  Z <- C <- rep(0, n); X <- matrix(NA, n, p)
  for(k in 1:K){
    if(nkA[k] > 0){
      X[nkCum[k]+(1:nk[k]), ] <-
        rbind(mvtnorm::rmvnorm(n=nkG[k], mean=muG[k, ], sigma=sigmaG[k, , ]),
              mvtnorm::rmvt(n=nkA[k], delta=muT[k, ], sigma=sigmaT[k, , ], df=dfT))
      C[nkCum[k]+nkG[k]+(1:nkA[k])] <- 1
    }else{
      X[nkCum[k]+(1:nk[k]), ] <-
        mvtnorm::rmvnorm(n=nk[k], mean=muG[k, ], sigma=sigmaG[k, , ])
    }
    Z[nkCum[k]+(1:nk[k])] <- k
  }
  randOrder <- rank(runif(n)); Z <- Z[randOrder]; C <- C[randOrder]; X <- X[randOrder, ];
  return(list(Z=Z, C=C, X=X))
}

#' @title Generation of Gaussian mixture distribution with uniform contamination
#'
#' @description Function to generate a sample of a gaussian mixture model with uniform contamination
#' \eqn{X\sim \sum_{k=1}^K \pi_k Y_k}, where \eqn{Y_k\sim \mathcal N (\mu_k, \Sigma_k)}.
#'
#' @param nk Vector of the number of observations in each component, giving the proportions \eqn{\pi_k}
#' @param muG Matrix of the means of the Gaussian components \eqn{\mu = (\mu_1, \ldots, \mu_K)} TODO plutot array
#' @param sigmaG Array of the variance-covariance matrices of the Gaussian components
#' @param delta Proportion of contamination
#' @param lUnif Vector of lower bounds of the uniform contamination
#' @param uUnif Vector of upper bounds of the uniform contamination
#'
#' @return A list containing the sample X, the indicator of contamination C and the latent variable Z, giving the classification
#' @examples
#' \donttest{
#' # Dimensions
#' K <- 3; nk <- rep(500, K); d <- 5
#'
#' # Gaussian parameters
#' mu <- rbind(mu1=rep(0, d), mu2=rep(3, d), mu3=rep(-3, d))
#' sigma <- array(dim=c(K, d, d))
#' sigma[1, , ] <- 2*diag(d); sigma[2, , ] <- diag(1:d); sigma[3, , ] <- diag(1/(1:d))
#'
#' # Student contamination
#' fraction <- 0.3
#' lower <- rep(0,d); upper <- rep(1,d)
#'
#' # Simulation
#' GMMcontUniform <- rGMMcontUniform(nk=nk, muG=mu, sigmaG=sigma,
#'                                   delta=fraction,
#'                                   lUnif=lower, uUnif=upper)
#' }
#'
#' @export
rGMMcontUniform <- function(nk, muG, sigmaG, delta, lUnif, uUnif){
  K <- length(nk); p <- ncol(muG); n <- sum(nk)
  nkCum <- c(0, cumsum(nk))
  nkG <- round(nk*(1-delta)); nkA <- nk - nkG
  Z <- C <- rep(0, n); X <- matrix(NA, n, p)
  for(k in 1:K){
    if(nkA[k] > 0){
      X[nkCum[k]+(1:nk[k]), ] <-
        rbind(mvtnorm::rmvnorm(n=nkG[k], mean=muG[k, ], sigma=sigmaG[k, , ]),
              sapply(1:p, function(j){runif(n=nkA[k], min=lUnif[j], max=uUnif[j])}))
      C[nkCum[k]+nkG[k]+(1:nkA[k])] <- 1
    }else{
      X[nkCum[k]+(1:nk[k]), ] <-
        mvtnorm::rmvnorm(n=nk[k], mean=muG[k, ], sigma=sigmaG[k, , ])
    }
    Z[nkCum[k]+(1:nk[k])] <- k
  }
  randOrder <- rank(runif(n)); Z <- Z[randOrder]; C <- C[randOrder]; X <- X[randOrder, ];
  return(list(Z=Z, C=C, X=X))
}


#' @title Generation of Student mixture distribution with Student contamination
#'
#' @description Function to generate a sample of a Student mixture model with Student contamination
#' \eqn{X\sim \sum_{k=1}^K \pi_k Y_k}, where \eqn{Y_k} follow a Student distribution with means \eqn{\mu = (\mu_1, \ldots, \mu_K)}
#' and variance-covaraince matrices \eqn{\Sigma = (\Sigma_1, \ldots, \Sigma_K)}.
#'
#' @param nk Vector of the number of observations in each component, giving the proportions \eqn{\pi_k}
#' @param dfT0 Degrees of freedom of the Student distribution
#' @param muT0 Matrix of the means of the Student distribution
#' @param sigmaT0 Array of the variances of the Student distribution
#' @param delta Proportion of contamination
#' @param dfT1 Degrees of freedom of the Student contamination
#' @param muT1 Matrix of the means of the Student contamination
#' @param sigmaT1 Array of the variances of the Student contamination
#'
#' @return A list containing the sample X, the indicator of contamination C and the latent variable Z, giving the classification
#'
#' @examples
#' \donttest{
#' # Dimensions
#' K <- 3; nk <- rep(500, K); d <- 5
#'
#' # Gaussian parameters
#' muT <- rbind(mu1=rep(0, d), mu2=rep(3, d), mu3=rep(-3, d))
#' sigmaT <- array(dim=c(K, d, d))
#' sigmaT[1, , ] <- 2*diag(d); sigmaT[2, , ] <- diag(1:d); sigmaT[3, , ] <- diag(1/(1:d))
#'
#' # Student contamination
#' fraction <- 0.3
#' dfT <- 1;
#'
#' # Simulation
#' TMMcontStudent <- rTMMcontStudent(nk=nk, dfT0=dfT, muT0=muT, sigmaT0=sigmaT,
#'                                   delta=fraction,
#'                                   dfT1=dfT, muT1=muT, sigmaT1=sigmaT)
#' }
#'
#' @export
rTMMcontStudent <- function(nk, dfT0, muT0, sigmaT0, delta, dfT1, muT1, sigmaT1){
  K <- length(nk); p <- ncol(muT0); n <- sum(nk)
  nkCum <- c(0, cumsum(nk))
  nkT0 <- round(nk*(1-delta)); nkA <- nk - nkT0
  Z <- C <- rep(0, n); X <- matrix(NA, n, p)
  for(k in 1:K){
    if(nkA[k] > 0){
      X[nkCum[k]+(1:nk[k]), ] <-
        rbind(mvtnorm::rmvt(n=nkT0[k], delta=muT0[k, ], sigma=sigmaT0[k, , ], df=dfT0),
              mvtnorm::rmvt(n=nkA[k], delta=muT1[k, ], sigma=sigmaT1[k, , ], df=dfT1))
      C[nkCum[k]+nkT0[k]+(1:nkA[k])] <- 1
    }else{
      X[nkCum[k]+(1:nk[k]), ] <-
        mvtnorm::rmvt(n=nk[k], delta=muT0[k, ], sigma=sigmaT0[k, , ], df=dfT0)
    }
    Z[nkCum[k]+(1:nk[k])] <- k
  }
  randOrder <- rank(runif(n)); Z <- Z[randOrder]; C <- C[randOrder]; X <- X[randOrder, ];
  return(list(Z=Z, C=C, X=X))
}


#' @title Generation of Student mixture distribution with uniform contamination
#'
#' @description Function to generate a sample of a Student mixture model with uniform contamination
#' \eqn{X\sim \sum_{k=1}^K \pi_k Y_k}, where \eqn{Y_k} follow a Student distribution with means \eqn{\mu = (\mu_1, \ldots, \mu_K)}
#' and variance-covaraince matrices \eqn{\Sigma = (\Sigma_1, \ldots, \Sigma_K)}.
#'
#' @param nk Vector of the number of observations in each component, giving the proportions \eqn{\pi_k}
#' @param dfT Degrees of freedom of the Student distribution
#' @param muT Matrix of the means of the Student distribution
#' @param sigmaT Array of the variances of the Student distribution
#' @param delta Proportion of contamination
#' @param lUnif Vector of lower bounds of the uniform contamination
#' @param uUnif Vector of uppers bounds of the uniform contamination
#'
#' @return A list containing the sample X, the indicator of contamination C and the latent variable Z, giving the classification
#'
#' @examples
#' \donttest{
#' # Dimensions
#' K <- 3; nk <- rep(500, K); d <- 5
#'
#' # Gaussian parameters
#' muT <- rbind(mu1=rep(0, d), mu2=rep(3, d), mu3=rep(-3, d))
#' sigmaT <- array(dim=c(K, d, d))
#' sigmaT[1, , ] <- 2*diag(d); sigmaT[2, , ] <- diag(1:d); sigmaT[3, , ] <- diag(1/(1:d))
#'
#' # Student contamination
#' fraction <- 0.3
#' dfT <- 1;
#' lower <- rep(0,d); upper <- rep(1,d)
#'
#' # Simulation
#' TMMcontUniform <- rTMMcontUniform(nk=nk, dfT=dfT, muT=muT, sigmaT=sigmaT,
#'                                   delta=fraction,
#'                                   lUnif=lower, uUnif=upper)
#' }
#'
#' @export
rTMMcontUniform <- function(nk, dfT, muT, sigmaT, delta, lUnif, uUnif){
  K <- length(nk); p <- ncol(muT); n <- sum(nk)
  nkCum <- c(0, cumsum(nk))
  nkT <- round(nk*(1-delta)); nkA <- nk - nkT
  Z <- C <- rep(0, n); X <- matrix(NA, n, p)
  for(k in 1:K){
    if(nkA[k] > 0){
      X[nkCum[k]+(1:nk[k]), ] <-
        rbind(mvtnorm::rmvt(n=nkT[k], delta=muT[k, ], sigma=sigmaT[k, , ], df=dfT),
              sapply(1:p, function(j){runif(n=nkA[k], min=lUnif[j], max=uUnif[j])}))
      C[nkCum[k]+nkT[k]+(1:nkA[k])] <- 1
    }else{
      X[nkCum[k]+(1:nk[k]), ] <-
        mvtnorm::rmvt(n=nk[k], delta=muT[k, ], sigma=sigmaT[k, , ], df=dfT)
    }
    Z[nkCum[k]+(1:nk[k])] <- k
  }
  randOrder <- rank(runif(n)); Z <- Z[randOrder]; C <- C[randOrder]; X <- X[randOrder, ];
  return(list(Z=Z, C=C, X=X))
}


#' Gaussian mixture on sphere sample generation
#'
#' @description Generate a sample of a Gaussian Mixture Model whose centers are generate randomly on a sphere of a given radius.
#'
#' @param n Number of points per cluster
#' @param d Dimension of the space
#' @param K Number of clusters
#' @param pcont Proportion of contamination
#' @param df Degrees of freedom of the Student contamination
#' @param cont Type of contamination: can be "Student" (default) or "Unif"
#' @param min Minimum value of the contamination
#' @param max Maximum value of the contamination
#' @param radius Radius of the sphere on which the centers are generated
#'
#' @examples
#' \donttest{
#' n <- 500
#' K <- 3
#' pcont <- 0.2
#'
#' sphereGM <- rSphereGM(n=n,K=K,pcont=pcont)
#' }
#'
#'@export
rSphereGM <- function(n=500,d=5,K=3,pcont=0,df=1,cont="Student",min=-5,max=5,radius=5)
{
  Sigma <- (diag(d)*2)
  X <- c()
  Tclassif <- c()
  mean <- rep(0,d)
  for (k in 1:K)
  {
    Z <- rnorm(d)
    Z <- radius*Z/sqrt(sum(Z^2))
    X <- rbind(X,mvtnorm::rmvnorm(n ,mean=Z))
    Tclassif <- c(Tclassif,rep(k,n))
  }

  if (pcont>0)
  {
    if(cont=='Unif')
    {
      Z <- matrix(runif(pcont*n*K*d,min,max),ncol=d)
    }
    if (cont=='Student')
    {
      Z <- matrix(rt(pcont*n*K*d,df=df),ncol=d)
    }
    I <- sample(1:(K*n),size=pcont*n*K)
    X[I,] <- Z
    Tclassif[I] <- "outliers"
  }

  mel <- sample.int(K*n)
  Tclassif <- Tclassif[mel]
  X <- X[mel,]
  return(list(X=X,cluster=Tclassif))
}
