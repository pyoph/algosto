# KL divergence btw F0 and F1

rm(list=ls()); par(pch=20, lwd=2)
setwd('/home/robin/Bureau/Paul/algosto/SR')
library(mvtnorm)
source('FunctionsKLgauss.R')
library('RGMM')

# Dims
d <- 10; KLval <- c(0, 1, 5, 10, 25)

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
xMin <- 1e-1; xMax <- 1e5
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

# Correction of the Mahalanobis
# Correction with the median does not work : 
# Due to the outliers, both the variance matrix and the median are upwarldy biased so
# 1 - the Mahalanobis distance is under-estimated
# 2 - the correction coeficient is smaller than one (which makes things worse) 
par(mfrow=c(1, 1))
n <- 1e5
kl <- 4; rate <- .4
n0 <- round((1-rate)*n); n1 <- n - n0
X0 <- rmvnorm(n0, mean=mu0, sigma=Sigma0)
parms1 <- ParmsF1(m1=m1, k1=k1val[kl], l1=l1val[kl], rho1=rho1val[kl])
if(n1 > 0){
  X1 <- rmvnorm(n1, mean=parms1$mu, sigma=parms1$Sigma)
  X <- rbind(X0, X1)
}else{X <- X0}
SigmaHat <- RobVar(X)$variance
plot(Sigma0, SigmaHat); abline(0, 1)
SigmaHatInv <- solve(SigmaHat)
Dmaha <- sapply(1:nrow(X), function(i){X[i, ]%*%SigmaHatInv%*%X[i, ]})
hist(Dmaha[1:n0], breaks=sqrt(nrow(X)), freq=FALSE)
curve(dchisq(x, df=d), from=0, to=100, col=2, add=TRUE)
alpha <- .5
coefCor <- qchisq(p=alpha, df=d) / quantile(Dmaha, probs=alpha)
Dcor <- Dmaha * coefCor
lines(density(Dcor[1:n0]), col=4)
abline(v=qchisq(p=.95, df=d), col=2)
c(mean(Dmaha[1:n0] > qchisq(p=.95, df=d)), 
  mean(Dcor[1:n0] > qchisq(p=.95, df=d)))
qqplot(qchisq(ppoints(n0), df=d), Dmaha[1:n0])
abline(0, 1)

U <- sapply(1:d, function(j){qnorm(rank(X[, j])/(n+1))})
Ucov <- cov(U)
sum(diag(Ucov%*%SigmaHatInv))
sum(diag(Ucov%*%SigmaHat))
