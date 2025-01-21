#install.packages("easystats")
#install.packages("bigutilsr")
#install.packages("xtable")
#install.packages("RGMM")
#install.packages("rrcov")
#install.packages("mvnfast")
library(xtable)
library(reshape2)
library(RobRegression)
library("robustbase")
library(Gmedian)
library(ggplot2)
library(far)
library(gridExtra)
library(microbenchmark)
library("matlib")
library("MASS")
library("corrplot")
library("dplyr")
library("easystats")
library("bigutilsr")
#library("rJava")
#library("REPPlab")
library("RGMM")
source("~/work/algosto/parametres.R")
source("~/work/algosto/simulations.R")
source("~/work/algosto/algorithmes.R")
source("~/work/algosto/resultats.R")
source("~/work/algosto/Outliers.R")
source("~/work/algosto/computeOutliers.R")

for (i in 1:nbruns)
{
  resultats_outliers <- calcule_outliers()
}


taux_contamination <- c(0,2, 5, 10, 15, 20, 25, 30, 40)


for (i in seq_along(taux_contamination)) 
{
  contamin = "moyenne"
  
  delta <- taux_contamination[i]
  
  
  #Génération des échantillons
  p1 <- 1 - delta / 100
  
  p2 <- 1 - p1

  resultsSimul <- genererEchantillon(n, d, mu1, mu2, p1, p2, Sigma1 = Sigma1, Sigma2 = Sigma2,contamin = contamin)
  Z <- resultsSimul$Z
  
  
  #Estimation online
  results <- estimMV(Z, c = sqrt(d),niterRMon = d ,methode = "eigen")
  SigmaEstimOnline <- results$Sigma
  
  #Estimation offline
  resultsOffline <- RobVar(Z)
  SigmaEstimOffline <- resultsOffline$variance
  
  erreursonline <- calculErreursNormeFrobenius(SigmaEstimOnline,Sigma1)
  
  erreursoffline <- calculErreursNormeFrobenius(SigmaEstimOffline,Sigma1)
  affiche_erreursSigmav2(erreursonline,erreursoffline,delta)
  
}