# Functions
KL <- function(parms1, parms2){
  invSigma2 <- solve(parms2$Sigma)
  0.5*(log(det(parms2$Sigma)/det(parms1$Sigma)) - d + sum(diag(invSigma2%*%parms1$Sigma)) +
         t(parms2$mu-parms1$mu)%*%invSigma2%*%(parms2$mu-parms1$mu))[1, 1]
}

# Contamination parms: F1
ParmsF1 <- function(m1, k1, l1, rho1){
  d <- length(m1)
  mu1 <- k1*m1
  sigmaSq1 <- l1*sigmaSq0
  Sigma1 <- diag(sqrt(sigmaSq1)) %*% toeplitz(rho1^(0:(d-1))) %*% diag(sqrt(sigmaSq1))
  return(list(mu=mu1, Sigma=Sigma1))
}

# Simulations
SimSample <- function(n, rate, parms0, parms1){
  x <- rmvnorm(n, mean=parms0$mu, sigma=parms0$Sigma)
  c <- rbinom(n, 1, prob=rate)
  if(sum(c) > 0){
    x1 <- rmvnorm(sum(c), mean=parms1$mu, sigma=parms1$Sigma)
    x[which(c==1), ] <- x1
  }
  return(list(x=x, c=c))
}
