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
  Sigma1 = diag(sqrt(1:10)),
  Sigma2 = 10 * diag(sqrt(1:10)),
  p1 = 1,
  p2 = 0,

  mu1 = rep(0,10 ),
  mu2 = 5 * rep(1, 10),
  r = 1,
  k = 1
)
params$n


#eigen(resultsSimul$Vvrai)$vectors

resultsSimul <- genererEchantillon(params$n,params$d,params$mu1,params$mu2,params$p1,params$p2,params$Sigma1,params$Sigma2)

#resultsSimul$labelsVrais

Z <- resultsSimul$Z
#Z
results <- estimMVOutliers(Z,params$c ,params$n,params$d,params$d,params$r)

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

affiche_erreurs_vectp_mcm(calculErreursVectpMCM(results$U,eigen(params$Sigma1)$vectors,params$n/params$k,params$d),params$n/params$k,params$c)
resultsSimul$VcovVrai
#Affichage de l'erreur d'estimation des valeurs propres de la matrice de covariance

affiche_erreurs_cov(calculErreursValPropreCov(results$lambdaIter,eigen(params$Sigma1)$values,params$n/params$k),params$n/params$k)

#Tracer avec les vraies valeurs 

#Comparaison stat khi-deux avec densité théorique Khi deux

densiteHistKhi2(results$stat,params$d)
#results$outlier_labels
table_contingence(resultsSimul$labelsVrais[1:999999],resultsStr$outlier_labels)

#calculErreursSigma(params$Sigma1,results$U,resultsStr$lambdatilde,params$n,params$d)

affiche_erreurs_Sigma(calculErreursSigma(params$Sigma1,results$U,results$lambdatilde,params$n/params$k,params$d),params$n/params$k,params$c)
results$U[9999,,] <- results$U[9999,,] %*% diag(1/sqrt(colSums(results$U[9999,,]^2)))
plot(params$Sigma1,SigmaEstim) 
abline(0,1)
 