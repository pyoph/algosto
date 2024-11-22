#install.packages("easystats")
#install.packages("bigutilsr")
#install.packages("xtable")
library(xtable)
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
library("easystats")
library("bigutilsr")
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


#Pourcentage de non outliers
p1 <- 0.98

#Pourcentage d'outliers
p2 <- 1-p1

resultsSimul <- genererEchantillon(n,d,mu1,mu2,p1,p2,Sigma1 =Sigma1 ,Sigma2)


Z <- resultsSimul$Z

outliers_list <- check_outliers(Z, method = "mahalanobis_robust")

table_contingence(resultsSimul$labelsVrais[1:9999],as.numeric(outliers_list)[1:9999])

###Fonction pour construire un tableau de résultats

taux_contamination <- c(2, 5, 10, 15, 20, 25, 30, 40)  
faux_positifs_maha <- numeric(length(taux_contamination))
faux_negatifs_maha <- numeric(length(taux_contamination))
faux_positifs_EPP <- numeric(length(taux_contamination))
faux_negatifs_EPP <- numeric(length(taux_contamination))
faux_positifs_online <- numeric(length(taux_contamination))
faux_negatifs_online <- numeric(length(taux_contamination))
i = 1
for (delta in taux_contamination)
{
  p1 <- 1 - delta/100
  p2 <- 1 - p1
  resultsSimul <- genererEchantillon(n,d,mu1,mu2,p1,p2,Sigma1 =Sigma1 ,Sigma2)
  
  
  Z <- resultsSimul$Z
  
  outliers_listMaha <- check_outliers(Z, method = "mahalanobis_robust")
  tc <- table_contingence(resultsSimul$labelsVrais[1:9999],as.numeric(outliers_listMaha)[1:9999])
  faux_positifs_maha[i] <- tc["1","0"]
  faux_negatifs_maha[i] <- tc["0","1"]   
  res.KurtM.Tribe <- EPPlab(Z, PPalg = "PSO", n.simu = 100,maxiter = 1000, sphere = TRUE)
  OUTms <- EPPlabOutlier(res.KurtM.Tribe, k = 1, location = median, scale = sd)
  
  outliersEPP <- OUTms$outlier
  tc <- table_contingence(resultsSimul$labelsVrais[1:9999],as.numeric(outliersEPP)[1:9999])
  
  faux_positifs_EPP[i] <- tc["1","0"]
  faux_negatifs_EPP[i] <- tc["0","1"]
  
  
  results <- estimMVOutliers(Z,params$c ,params$n,params$d,params$d,params$r,aa = 0.75,niter = 1e4)
  
  tc <- table_contingence(resultsSimul$labelsVrais[1:9999],results$outlier_labels[1:9999])
  
  faux_positifs_online[i] <- tc["1","0"]
  faux_negatifs_online[i] <- tc["0","1"]
   
  i = i + 1
}


# Créer un data.frame final
results_outlier_df <- data.frame(
  #"Taux de Contamination (%)" = taux_contamination,
  "FP_Mahalanobis" = faux_positifs_maha,
  "FN_Mahalanobis" = faux_negatifs_maha,
  "FP_EPP" = faux_positifs_EPP,
  "FN_EPP" = faux_negatifs_EPP,
  "FP_Online" = faux_positifs_online,
  "FN_Online" = faux_negatifs_online
)


