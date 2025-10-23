# KL divergence btw F0 and F1

rm(list=ls())
source('FunctionsKLgauss.R')

# Dims
d <- 100; KLval <- c(0, 1, 10, 100)
d <- 10; KLval <- c(0, 1, 10, 33, 100)
d <- 10; KLval <- c(0, 1, 5, 10, 25, 50, 100)

# Null parms: F0
mu0 <- rep(0, d)
sigmaSq0 <- (1:d); sigmaSq0 <- sigmaSq0 / mean(sigmaSq0)
rho0 <- 0.3
Sigma0 <- diag(sqrt(sigmaSq0)) %*% toeplitz(rho0^(0:(d-1))) %*% diag(sqrt(sigmaSq0))
parms0 <- list(mu=mu0, Sigma=Sigma0)

# Contamination
m1 <- (-1)^(1:d)/sqrt(d)
#KL(parms1=parms0, parms2=ParmsF1(m1, 1, 1, 0.5))

# Setting simulation parms
k1grid <- c(0, 2^seq(-4, 5, by=.001)); # k1val <- c(0, 2, 5, 10)
KLk1 <- sapply(k1grid, function(k1){KL(parms0, ParmsF1(m1, k1, 1, rho0))})
k1val <- k1grid[sapply(KLval, function(kl){which.min(abs(KLk1 - kl))})]
k1val

l1grid <- 2^seq(-5, 0, .01); # l1val <- c(1/20, .2, .5, 1, 2, 20)
KLl1 <- sapply(l1grid, function(l1){KL(parms0, ParmsF1(m1, 0, l1, rho0))})
l1val <- l1grid[sapply(KLval, function(kl){which.min(abs(KLl1 - kl))})]
l1grid <- 2^seq(0, 40, .01); # l1val <- c(1/20, .2, .5, 1, 2, 20)
KLl1 <- sapply(l1grid, function(l1){KL(parms0, ParmsF1(m1, 0, l1, rho0))})
l1val <- round(c(l1val, l1grid[sapply(KLval, function(kl){which.min(abs(KLl1 - kl))})]), 3)
l1val <- unique(l1val)
l1val
l1grid <- 2^seq(-5, 10, 0.01)

rho1grid <- seq(rho0, 0.995, by=0.005); #= rho1val <- c(0, .3, .6, 0.8, .95)
KLrho1 <- sapply(rho1grid, function(rho1){KL(parms0, ParmsF1(m1, 0, 1, rho1))})
rho1val <- rho1grid[sapply(KLval, function(kl){which.min(abs(KLrho1 - kl))})]
rho1val
rho1grid <- seq(-0.995, rho0, by=0.005); #= rho1val <- c(0, .3, .6, 0.8, .95)
KLrho1 <- sapply(rho1grid, function(rho1){KL(parms0, ParmsF1(m1, 0, 1, rho1))})
rho1val <- round(c(rho1val, rho1grid[sapply(KLval, function(kl){which.min(abs(KLrho1 - kl))})]), 3)
rho1val <- unique(rho1val)
rho1val
rho1grid <- seq(-0.995, 0.995, by=0.005); #= rho1val <- c(0, .3, .6, 0.8, .95)

# Export
save(d, KLval, 
     mu0, sigmaSq0, Sigma0, rho0, 
     m1, 
     k1val, l1val, rho1val, 
     file=paste0('SimParmsGrid-d', d, '.Rdata'))

# # Plot
# par(mfrow=c(3, 1))
# plot(k1grid, sapply(k1grid, function(k1){KL(parms0, ParmsF1(m1, k1, 1, rho0))}), log='xy')
# abline(h=KLval); abline(v=k1val)
# plot(l1grid, sapply(l1grid, function(l1){KL(parms0, ParmsF1(m1, 0, l1, rho0))}), log='xy')
# abline(h=KLval); abline(v=l1val)
# plot(rho1grid, sapply(rho1grid, function(rho1){KL(parms0, ParmsF1(m1, 0, 1, rho1))}), log='y')
# abline(h=KLval); abline(v=rho1val)
# 
