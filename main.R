#install.packages("easystats")
#install.packages("bigutilsr")
#install.packages("xtable")
#install.packages("RGMM")
#install.packages("rrcov")
#install.packages("mvnfast")
#install.packages("reshape")
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
source("~/work/algosto/seuils.R")
for (i in 1:nbruns)
{
  resultats_outliers <- calcule_outliers(contamin = "moyenne")
}


taux_contamination <- c(0,2, 5, 10, 15, 20, 25, 30, 40)

###Test fonction d'estimation online

for (i in seq_along(taux_contamination)) 
{
  contamin = "moyenne"
  
  delta <- taux_contamination[i]
  #delta <- 2
  
  #Génération des échantillons
  p1 <- 1 - delta / 100
  
  p2 <- 1 - p1

  resultsSimul <- genererEchantillon(n, d, mu1, mu2, p1, p2, Sigma1 = Sigma1, Sigma2 = Sigma2,contamin = contamin)
  Z <- resultsSimul$Z
  
  
  #Estimation online
  results <- estimMV(Z, c = sqrt(d),niterRMon = d ,methode = "eigen")
  SigmaEstimOnline <- results$Sigma
  distances <- results$distances
  #Estimation offline
  resultsOffline <- RobVar(Z)
  SigmaEstimOffline <- resultsOffline$variance
  
  erreursonline <- calculErreursNormeFrobenius(SigmaEstimOnline,Sigma1)
  
  #erreursoffline <- calculErreursNormeFrobenius(SigmaEstimOffline,Sigma1)
  #affiche_erreursSigmav2(erreursonline,erreursoffline = rep(0,,nrow(Z)),delta)
  
}

###Calcul du meilleur AUC pour chaque méthode

methodes <- c("Cov Empirique", "OGK", "Online", "Offline", "Comédiane", "Shrinkage")

taux_contamination <- c(2, 5, 10, 15, 20, 25, 30, 40)

#Création d'un dataframe pour contenir les meilleurs AUC

auc_df <- data.frame(matrix(ncol = length(methodes), nrow = length(taux_contamination)))

rownames(auc_df) <- taux_contamination

colnames(auc_df) <- methodes

for (i in seq_along(taux_contamination)) 
{
  contamin = "moyenne"
  
  delta <- taux_contamination[i]
  #delta <- 2
  
  #Génération des échantillons
  p1 <- 1 - delta / 100
  
  p2 <- 1 - p1
  
  resultsSimul <- genererEchantillon(n, d, mu1, mu2, p1, p2, Sigma1 = Sigma1, Sigma2 = Sigma2,contamin = contamin)
  Z <- resultsSimul$Z
  
  for (m in methodes)
  {
    
      distances <- calcule_distances_par_methode(Z,methode = m)
      auc_df[i,m] <- courbeROC(resultsSimul$labelsVrais,distances)
  }
  
  
  
}
