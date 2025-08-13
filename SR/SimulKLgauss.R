# Simul : contamination of F0 with F1

rm(list=ls())
library(mvtnorm); library(fields)
simDir <- 'Simuls/'
figDir <- 'Figures/'
exportFig <- TRUE

# Functions
source('FunctionsKLgauss.R')

# Dims
d <- 10

# Parms
k1List <- c(0, 2, 5, 10); k1Nb <- length(k1List)
l1List <- c(1/20, .2, .5, 1, 2, 20); l1Nb <- length(l1List)
rho1List <- c(-.95, -.7, -.15, .3, .6, .85, .975); rho1Nb <- length(rho1List)
rList <- seq(0, .5, by=.05); rNb <- length(rList)

# Parms
mu0 <- rep(0, d)
sigmaSq0 <- (1:d); sigmaSq0 <- sigmaSq0 / mean(sigmaSq0)
rho0 <- 0.3
Sigma0 <- diag(sqrt(sigmaSq0)) %*% toeplitz(rho0^(0:(d-1))) %*% diag(sqrt(sigmaSq0))
parms0 <- list(mu=mu0, Sigma=Sigma0)
m1 <- rep(1/sqrt(d), d)

# RMSE
B <- 100; n <- 1e4
parmList <- k1List; parmNb <- k1Nb; parmName = 'k1'
# parmList <- l1List; parmNb <- l1Nb; parmName = 'lambda'
# parmList <- rho1List; parmNb <- rho1Nb; parmName = 'rho1'
resFile <- paste0(simDir, 'SimContamination-', parmName, '-B', B, '.Rdata')
if(!file.exists(resFile)){
  kl <- rmse <- array(NA, dim=c(B, rNb, parmNb), dimnames=list(1:B, rList, parmList))
  klTheo  <- rmseTheo <- rep(NA, parmNb)
  for(pp in 1:parmNb){
    cat(parmList[pp], '')
    if(parmName=='k1'){parms1 <- ParmsF1(m1, parmList[pp], 1, rho0)}
    if(parmName=='lambda'){parms1 <- ParmsF1(m1, 0, parmList[pp], rho0)}
    if(parmName=='rho1'){parms1 <- ParmsF1(m1, 0, 1, parmList[pp])}
    rmseTheo[pp] <- FrobeniusNormError(parms1$Sigma1, Sigma0)           
    klTheo[pp] <- KL(parms1=parms0, parms2=parms1)           
    for(rr in 1:rNb){
      n1 <- round(rList[rr]*n); n0 <- n - n1
      for(b in 1:B){
        if(n1 > 0){
          x <- rbind(rmvnorm(n0, mean=mu0, sigma=Sigma0), 
                     rmvnorm(n1, mean=parms1$mu1, sigma=parms1$Sigma1))
        }else{
          x <- rmvnorm(n, mean=mu0, sigma=Sigma0)
        }
        sigmaHat <- cov(x)
        muHat <- mean(x)
        rmse[b, rr, pp] <- FrobeniusNormError(sigmaHat, Sigma0)           
        parmsHat <- list(mu=muHat, Sigma=sigmaHat)
        kl[b, rr, pp] <- KL(parms1=parms0, parms2=parmsHat)           
      }
    }
  }
  cat('\n')
  save(parmName, parmList, parmNb, klTheo, kl, rmseTheo, rmse, 
       file=resFile)
}else{load(resFile)}

figFile <- paste0(figDir, 'SimContamination-', parmName, '-B', B, '.png')
if(exportFig){png(figFile)}
par(mfcol=c(2, parmNb))
for(pp in 1:parmNb){
  boxplot(rmse[, , pp], main=paste(parmName, '=', parmList[pp]), 
          ylim=c(min(rmse), max(max(rmse), max(rmseTheo))), , log='y', ylab='RMSE')
  abline(h=rmseTheo[pp], lty=2, col=2, lwd=2)
  boxplot(kl[, , pp], main=paste(parmName, '=', parmList[pp]), 
          ylim=c(min(kl), max(max(kl), max(klTheo))), log='y', ylab='KL')
  abline(h=klTheo[pp], lty=2, col=2, lwd=2)
}
if(exportFig){dev.off()}
