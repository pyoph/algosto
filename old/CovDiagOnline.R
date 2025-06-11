covDiagOnline <- function(Y, k = -1) {
  dim.Y <- dim(Y)
  N <- dim.Y[1]
  p <- dim.Y[2]
  sample= matrix(0,p,p)
  sample2= matrix(0,p,p)
  if (k < 0) {    # demean the data and set k = 1
    Y <- scale(Y, scale = F)
    k <- 1
  }
  n <- N - k    # effective sample size
  c <- p / n    # concentration ratio
  Y2 = Y^2
  for (i in N){
    #print((Y[i,]%*%(t(Y[i,]))))
    sample = sample  + (1.0 / (n + 1)) * (Y[i,] %*% t(Y[i,]) - sample )
    sample2 = sample2  + (1.0 / (n + 1)) * (Y2[i,] %*% t(Y2[i,]) - sample2 )
    piMat <- sample2 - sample^2
    pihat <- sum(piMat)
    
    target <- diag(diag(sample))
    gammahat <- norm(c(sample - target), type = "2")^2
    
    # diagonal part of the parameter that we call rho 
    rho_diag <- sum(diag(piMat))
    
    # off-diagonal part of the parameter that we call rho 
    rho_off <- 0
    
    # compute shrinkage intensity
    rhohat <- rho_diag + rho_off
    kappahat <- (pihat - rhohat) / gammahat
    shrinkage <- max(0, min(1, kappahat / n))
    
    # compute shrinkage estimator
    sigmahat <- shrinkage * target + (1 - shrinkage) * sample   
  }
 return(sigmahat = sigmahat) 
 
}