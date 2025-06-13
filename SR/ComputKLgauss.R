# KL divergence btw F0 and F1

rm(list=ls())
library(nlme); library(flexmix); library(fields)

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
  return(list(mu1=mu1, Sigma1=Sigma1))
}


# Dims
d <- 10

# Null parms: F0
mu0 <- rep(0, d)
sigmaSq0 <- (1:d); sigmaSq0 <- sigmaSq0 / mean(sigmaSq0)
rho0 <- 0.3
Sigma0 <- diag(sqrt(sigmaSq0)) %*% toeplitz(rho0^(0:(d-1))) %*% diag(sqrt(sigmaSq0))
parms0 <- list(mu=mu0, Sigma=Sigma0)

# # Null (null) parms: F0
# mu0 <- rep(0, d)
# sigmaSq0 <- rep(1, d)
# rho0 <- 0
# Sigma0 <- diag(sqrt(sigmaSq0)) %*% toeplitz(rho0^(0:(d-1))) %*% diag(sqrt(sigmaSq0))
# parms0 <- list(mu=mu0, Sigma=Sigma0)

# Contamination
m1 <- rep(1/sqrt(d), d)
KL(parms1=parms0, parms2=ParmsF1(m1, 1, 1, 0.5))

par(mfrow=c(3, 2))
k1grid <- 2^(-4:5); 
plot(k1grid, sapply(k1grid, function(k1){KL(parms0, ParmsF1(m1, k1, 1, rho0))}), log='xy')
l1grid <- 2^(-5:5); 
plot(l1grid, sapply(l1grid, function(l1){KL(parms0, ParmsF1(m1, 0, l1, rho0))}), log='xy')
rho1grid <- seq(0.005, 0.995, by=0.005); 
plot(rho1grid, sapply(rho1grid, function(rho1){KL(parms0, ParmsF1(m1, 0, 1, rho1))}), log='y')

KLk1l1 <- sapply(k1grid, function(k1){sapply(l1grid, function(l1){KL(parms0, ParmsF1(m1, k1, l1, rho0))})})
image.plot(log10(l1grid), log10(k1grid), log10(KLk1l1))
KLk1rho1 <- sapply(k1grid, function(k1){sapply(rho1grid, function(rho1){KL(parms0, ParmsF1(m1, k1, 1, rho1))})})
image.plot(rho1grid, log10(k1grid), log10(KLk1rho1))
KLrho1l1 <- sapply(rho1grid, function(rho1){sapply(l1grid, function(l1){KL(parms0, ParmsF1(m1, 0, l1, rho1))})})
KLrho1l1[KLrho1l1==0] <- NA
image.plot(log10(l1grid), rho1grid, log10(KLrho1l1))
