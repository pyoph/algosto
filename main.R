#install.packages("microbenchmark")
#install.packages("matlib")
#install.packages("MASS")
#install.packages("reshape2")
#install.packages("corrplot")
#install.packages("dlpyr")
#install.packages("devtools")
#install.packages("usethis")

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
source("~/work/algosto/parametres.R")
source("~/work/algosto/simulations.R")
source("~/work/algosto/algorithmes.R")
source("~/work/algosto/resultats.R")
#usethis::create_package("C:/Users/Paul/Documents/codeRTheseLPSM")
#devtools::load_all("C:/Users/Paul/Documents/codeRTheseLPSM")


params <- initialiser_parametres(
  d = 10,
  q = 10,
  n = 1e6,
  c = sqrt(10),
  m0 = rep(0, 10),
  Sigma2 = diag(1:10),#Tester une matrice de covariance non diagonale
  Sigma1 =  matrix(c(2, 0.8, 0.5, 0.3, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005,
                     0.8, 3, 0.7, 0.4, 0.3, 0.2, 0.1, 0.05, 0.02, 0.01,
                     0.5, 0.7, 1.5, 0.6, 0.5, 0.3, 0.2, 0.1, 0.05, 0.02,
                     0.3, 0.4, 0.6, 2.5, 0.7, 0.5, 0.3, 0.2, 0.1, 0.05,
                     0.2, 0.3, 0.5, 0.7, 1.8, 0.8, 0.5, 0.3, 0.2, 0.1,
                     0.1, 0.2, 0.3, 0.5, 0.8, 2.2, 0.7, 0.5, 0.3, 0.2,
                     0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.7, 0.6, 0.5, 0.3,
                     0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.6, 2.1, 0.7, 0.5,
                     0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.9, 0.8,
                     0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 2.5),
                   nrow = 10, ncol = 10, byrow = TRUE),
  #Sigma3 = 10 * diag((1:10)),
  p1 = 1,
  p2 = 0,
  
  mu1 = rep(0,10 ),
  mu2 = 5 * rep(1, 10),
  r = 1,
  k = 1
)
params$n


# Creation d'une matrice de covariance
A <- matrix(runif(100, min=-1, max=1), nrow=10)
A <- (A + t(A)) / 2 
cov_matrix <- A %*% t(A)

cov_matrix

#eigen(resultsSimul$Vvrai)$vectors

resultsSimul <- genererEchantillon(params$n,params$d,params$mu1,params$mu2,params$p1,params$p2,cov_matrix,params$Sigma2)

#resultsSimul$labelsVrais

Z <- resultsSimul$Z
#Z
results <- estimMVOutliers(Z,params$c ,params$n,params$d,params$d,params$r,aa = 0.75)



  
resultsStr <- streaming(Z,params$n/params$k,params$k,params$c,params$n,params$d,params$q,params$r)

#lambdaResultat <- RobbinsMC2(1e4,results$vpMCM[99,],10^(-8),alpha = 0.75,c= 2,w=2,samp=1e4, initbarre = results$vpMCM[99,],cbarre=0)

#lambdaResultat$lambdalist

#eigen(params$Sigma1)$values

#Affichage de l'estimation des erreurs pour l'estimation de mbarre

affiche_erreurs_mediane_geometrique(calculErreursM(results$miter,params$m,params$n/params$k),params$n/params$k)

#Affichage de l'estimation des erreurs pour l'estimation de Vbarre

affiche_erreurs_mcm(calculErreursV(results$VIter,resultsSimul$Vvrai,params$n/params$k),params$n/params$k,params$c)

#Affichage de l'erreur d'estimation des valeurs propres de la MCM

affiche_erreurs_vp_mcm(calculErreursVpMCM(results$vpMCM,resultsSimul$VpvraiesV,params$n/params$k),params$n/params$k)

#results$U[99,,]

#Affichage de l'erreur d'estimation des vecteurs propres de la MCM

affiche_erreurs_vectp_mcm(calculErreursVectpMCM(results$U,eigen(cov_matrix)$vectors,params$n/params$k,params$d),params$n/params$k,params$c)
#resultsSimul$VcovVrai
#Affichage de l'erreur d'estimation des valeurs propres de la matrice de covariance
affiche_erreurs_cov(calculErreursValPropreCov(results$lambdaIter,eigen(cov_matrix)$values,params$n/params$k),params$n/params$k)
results$U[999,,] 
#Tracer avec les vraies valeurs 

#Comparaison stat khi-deux avec densité théorique Khi deux

densiteHistKhi2(results$stat,params$d)
#results$outlier_labels
table_contingence(resultsSimul$labelsVrais[1:999999],results$outlier_labels)

#calculErreursSigma(params$Sigma1,results$U,resultsStr$lambdatilde,params$n,params$d)

affiche_erreurs_Sigma(calculErreursSigma(cov_matrix,results$U,results$lambdatilde,params$n/params$k,params$d),params$n/params$k,params$c)
results$U[999999,,] <- results$U[999999,,] %*% diag(1/sqrt(colSums(results$U[999999,,]^2)))
results$U[999999,,] %*% diag(1/sqrt(colSums(results$U[999999,,]^2)))
SigmaEstim <- t(results$U[999999,,])%*%diag(results$lambdaIter[999999,])%*%results$U[999999,,]

plot(cov_matrix,SigmaEstim) 
abline(0,1)

statTheorique <- calculeStatChi2(resultsSimul$Z,params$Sigma1)

densiteHistKhi2(statTheorique,params$d)


