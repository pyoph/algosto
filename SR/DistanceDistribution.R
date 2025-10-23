# KL divergence btw F0 and F1

rm(list=ls()); par(pch=20, lwd=2)
library(mvtnorm)
source('FunctionsKLgauss.R')

# Dims
d <- 10; KLval <- c(0, 1, 5, 10, 25, 50, 100)

# Parms
load(paste0('SimParmsGrid-d', d, '.Rdata'))
l1val <- sort(l1val[which(l1val >= 1)])
rho1val <- sort(rho1val[which(rho1val >= rho0)])
rbind(KLval, k1val, l1val, rho1val)
klNb <- length(KLval)

# Sigma0
Sigma0inv <- solve(Sigma0)
Sigma0eig <- eigen(Sigma0)
Sigma0invHalf <- Sigma0eig$vectors%*%diag(1/sqrt(Sigma0eig$values))%*%t(Sigma0eig$vectors)
Sigma0half <- Sigma0eig$vectors%*%diag(sqrt(Sigma0eig$values))%*%t(Sigma0eig$vectors)

par(mfrow=c(2, 1))
B <- 1e6

# Distribution 
U <- matrix(rnorm(d*B), B, d)
X0 <- rmvnorm(B, mean=mu0, sigma=Sigma0)
d0 <- rowSums((X0%*%Sigma0invHalf)^2)
f0 <- density((d0))
xMin <- 1e-1; xMax <- 1e2
plot(f0, xlim=c((xMin), (xMax)), ylim=c(0, .1), main='Mahalonobis')  
# curve(dchisq(x, df=d), from=1e-2, to=max(d0), add=TRUE, lty=2)
text(0.98*f0$x[which.max(f0$y)], 1.02*max(f0$y), labels='0', col=1)
abline(v=(quantile(d0, probs=.95)), lty=2)
for(kl in 2:klNb){
  parms1 <- ParmsF1(m1=m1, k1=k1val[kl], l1=l1val[kl], rho1=rho1val[kl])
  X1 <- rmvnorm(B, mean=parms1$mu, sigma=parms1$Sigma)
  d1 <- rowSums((X1%*%Sigma0invHalf)^2)
  f1 <- density((d1))
  lines(f1, col=kl)  
  text(0.98*f1$x[which.max(f1$y)], 1.02*max(f1$y), labels=KLval[kl], col=kl)
}

# Distribution (log scale)
U <- matrix(rnorm(d*B), B, d)
X0 <- rmvnorm(B, mean=mu0, sigma=Sigma0)
d0 <- rowSums((X0%*%Sigma0invHalf)^2)
f0 <- density(log10(d0))
xMin <- 1e-1; xMax <- 1e11
plot(f0, xlim=c(log10(xMin), log10(xMax)), ylim=c(0, 2.2), main='log10-Mahalonobis')  
# curve(dchisq(x, df=d), from=1e-2, to=max(d0), add=TRUE, lty=2)
text(0.98*f0$x[which.max(f0$y)], 1.02*max(f0$y), labels='0', col=1)
abline(v=log10(quantile(d0, probs=.95)), lty=2)
for(kl in 2:klNb){
  parms1 <- ParmsF1(m1=m1, k1=k1val[kl], l1=l1val[kl], rho1=rho1val[kl])
  X1 <- rmvnorm(B, mean=parms1$mu, sigma=parms1$Sigma)
  d1 <- rowSums((X1%*%Sigma0invHalf)^2)
  f1 <- density(log10(d1))
  lines(f1, col=kl)  
  text(0.98*f1$x[which.max(f1$y)], 1.02*max(f1$y), labels=KLval[kl], col=kl)
}

