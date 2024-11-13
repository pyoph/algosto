#install.packages("microbenchmark")
#install.packages("matlib")
#install.packages("MASS")
#install.packages("reshape2")
#install.packages("corrplot")
#install.packages("dlpyr")
#install.packages("devtools")
#install.packages("usethis")
#install.packages("REPPlab")
#install.packages("rJava")

library(reshape2)
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
library("rJava")
library("REPPlab")
source("~/algosto/parametres.R")
source("~/algosto/simulations.R")
source("~/algosto/algorithmes.R")
source("~/algosto/resultats.R")
#usethis::create_package("C:/Users/Paul/Documents/codeRTheseLPSM")
#devtools::load_all("C:/Users/Paul/Documents/codeRTheseLPSM")


n <- 1e4

d <- 10

mu1 <- rep(0,d)

mu2 <- 5*rep(1,d)

Sigma1 <- diag(sqrt(1:10))


Sigma2 <- 10*diag(1:10)

set.seed(123)

#Création d'une matrice de covariance non nécessairement diagonale
A <- matrix(runif(100, min=-1, max=1), nrow=10)
A <- (A + t(A)) / 2  
Sigma1 <- A %*% t(A)
Sigma1=(Sigma1+diag(10))/10
eigen(Sigma1)$values

#Pourcentage de non outliers
p1 <- 1

#Pourcentage d'outliers
p2 <- 1-p1

resultsSimul <- genererEchantillon(n,d,mu1,mu2,p1,p2,Sigma1 =Sigma1 ,Sigma2)


Z <- resultsSimul$Z

set.seed(4567)

res.KurtM.Tribe <- EPPlab(Z, PPalg = "Tribe", n.simu = 100,maxiter = 200, sphere = TRUE)

OUTms <- EPPlabOutlier(res.KurtM.Tribe, k = 1, location = median, scale = sd)
help(EPPlabOutlier)
summary(OUTms)

plot(OUTms, las = 1)
