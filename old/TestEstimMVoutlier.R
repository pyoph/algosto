# Test algo robustes 

rm(list=ls())
source('LoadAll.R')
# source("../algosto/parametres.R")
# source("../algosto/simulations.R")
# source("../algosto/algorithmes.R")
# source("../algosto/resultats.R")
# source("../algosto/Outliers.R")
# source("../algosto/computeOutliers.R")
# source("../algosto/seuils.R")
library(mvtnorm)
library(ROCR)

# Parms gene
# seed <- 1; set.seed(seed)
par(mfrow=c(1, 1), pch=20); palette('R3')

# Parms 
n <- 5e3; d <- 10; p1 <- .01
sigmaTrue <- diag(d); omegaTrue <- solve(sigmaTrue)
sTrue <- sqrt(rgamma(d, 0.5, 0.5))
corTrue <- exp(-as.matrix(dist(matrix(rnorm(2*d), d, 2)), diag=TRUE))
sigmaTrue <- diag(sTrue)%*%corTrue%*%diag(sTrue); omegaTrue <- solve(sigmaTrue)
eigen(sigmaTrue)$values

# Simuls
Y <- rmvnorm(n, sigma=sigmaTrue)
Z <- rbinom(n, 1, prob=p1); n1 <- sum(Z)
Y[which(Z==1), ] <- rmvt(n1, sigma=2*sigmaTrue, df=1)
# Y[which(Z==1), ] <- 5*matrix(runif(n1*d), n1, d)
thresh <- qchisq(0.95, df=d)

# Fit
sigmaNaive <- cov(Y); muNaive <- colMeans(Y); omegaNaive <- solve(sigmaNaive)
resOn <- estimMV(Y=Y)
sigmaOn <- resOn$Sigma[dim(resOn$Sigma)[1]-1, , ]; muOn <- resOn$moyennem; omegaOn <- solve(sigmaOn)
resOff <- RobVar(X=Y)
sigmaOff <- resOff$variance; muOff <- as.vector(resOff$median); omegaOff <- solve(sigmaOff)
plot(sigmaTrue, sigmaOn, col=2); abline(0, 1)
points(sigmaTrue, sigmaOff, col=4)
plot(omegaTrue, omegaOn, col=2); abline(0, 1)
points(omegaTrue, omegaOff, col=4)

# Distance
dist2True <- sapply(1:n, function(i){(Y[i, ]%*%omegaTrue%*%Y[i, ])[1, 1]})
dist2Naive <- sapply(1:n, function(i){(Y[i, ] - muNaive)%*%omegaNaive%*%(Y[i, ] - muNaive)})
dist2On <- sapply(1:n, function(i){(Y[i, ] - muOn)%*%omegaOn%*%(Y[i, ] - muOn)})
dist2Off <- sapply(1:n, function(i){(Y[i, ] - muOff)%*%omegaOff%*%(Y[i, ] - muOff)})
plot(dist2True, dist2On, col=2); abline(0, 1)
points(dist2True, dist2Off, col=4)

# Outliers
Performance <- function(score, class, thresh, add=FALSE, col=1){
  perf <- performance(prediction(score, class), "tpr", "fpr")
  tab <- cbind(unlist(perf@alpha.values), unlist(perf@x.values), unlist(perf@y.values))
  numThresh <- min(which(unlist(perf@alpha.values) < thresh))
  perfThresh <- tab[numThresh, -1]; 
  names(perfThresh) <- c("tpr", "fpr")
  plot(perf, add=add, col=col); abline(0, 1, col=8)
  points(perfThresh[1], perfThresh[2], col=col, cex=2)
  return(list(perf=perf, tab=tab, perfThresh=perfThresh))
}

perfTrue <- Performance(score=dist2True, class=Z, thresh=thresh)
perfNaive <- Performance(score=dist2Naive, class=Z, thresh=thresh, add=TRUE, col=3)
perfOn <- Performance(score=dist2On, class=Z, thresh=thresh, add=TRUE, col=2)
perfOff <- Performance(score=dist2Off, class=Z, thresh=thresh, add=TRUE, col=4)

rbind(perfTrue$perfThresh, perfOn$perfThresh, perfOff$perfThresh, perfNaive$perfThresh)
