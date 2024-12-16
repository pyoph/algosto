#install.packages("microbenchmark")
#install.packages("matlib")
#install.packages("MASS")
#install.packages("reshape2")
#install.packages("corrplot")
#install.packages("dlpyr")
#install.packages("robustbase")
#install.packages("remotes")
#install.packages("rrcov")
library(rrcov)
library(reshape2)
library("robustbase")
library(RobRegression)
library(Gmedian)
library(ggplot2)
library(far)
library(gridExtra)
library(microbenchmark)
library("matlib")
library("MASS")
library("corrplot")
library("dplyr")
source("~/algosto/parametres.R")
source("~/algosto/simulations.R")
source("~/algosto/algorithmes.R")
source("~/algosto/resultats.R")
source("~/algosto/seuils.R")



n <- 1e4

d <- 10

rho <- 0.8

seuil_p_value <- 0.05

distances <- rep(0,n)

Sigma1 <- creerMatriceToeplitz(rho,d)

Sigma2 <- permuterLignesColonnes(Sigma1,lignes_a_permuter = c(1,d),colonnes_a_permuter = c(1,d))

#Sigma2 <- creerMatriceToeplitz(0.4,d)
#Sigma1 <- diag(sqrt(1:d))
#Sigma1 <- diag(d)
mu1 <- rep(0,d)

mu2 <- 5*rep(1,d)




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
  sample <- (t(Y) %*% Y) / n   
  #sample <- covComed(Z)
  # compute shrinkage target
  samplevar <- diag(sample)
  sqrtvar <- sqrt(samplevar)
  rBar <- (sum(sample / outer(sqrtvar, sqrtvar)) - p) / (p * (p - 1))
  target <- rBar * outer(sqrtvar, sqrtvar)
  diag(target) <- samplevar
  
  # estimate the parameter that we call pi in Ledoit and Wolf (2003, JEF)
  Y2 <- Y^2
  sample2 <- (t(Y2) %*% Y2) / n   
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


p1 <- 0.98
p2 <- 1 - p1

Est <- covComed(Z)

Z

#Est2 <- cov.comed(Z)


m <- covComed(Z)


m$center

resultsSimul <- genererEchantillon(n, d, mu1, mu2, p1, p2, Sigma1 = Sigma1, Sigma2 = Sigma2)

#cov(Z)

Z <- resultsSimul$Z

sigmahat <- covCor(Z)


Com <- Est$cov

Est$



erreursSigmaMoyenne <- rep(0, n - 1)
erreursSigmaEstim <- rep(0, n - 1)


for (i in 1:(n-1))
{
  
  #SigmaEstim <- diag(lambda[i,])
  erreursSigmaEstim[i] <- (norm(Com - Sigma1,type = "F"))/i
  erreursSigmaMoyenne[i] <- mean(erreursSigmaEstim[1:i])
  
}

affiche_erreurs_Sigma(erreursSigmaMoyenne, n = 1e4)
