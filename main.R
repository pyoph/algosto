#install.packages("easystats")
#install.packages("bigutilsr")
#install.packages("xtable")
#install.packages("RGMM")
#install.packages("rrcov")
#install.packages("mvnfast")
#install.packages("reshape")
#install.packages("matlib")
library(xtable)
library(reshape2)
library(RobRegression)
library("robustbase")
library(Gmedian)
library(ggplot2)
library(far)
library(gridExtra)
library(microbenchmark)
#library("matlib")
library("MASS")
library("corrplot")
library("dplyr")
library("easystats")
library("bigutilsr")
#library("rJava")
#library("REPPlab")
library("RGMM")
library(ROCR)
library(pROC)
source("~/codeThese/algosto/parametres.R")
source("~/codeThese//algosto/simulations.R")
source("~/codeThese/algosto/algorithmes.R")
source("~/codeThese/algosto/resultats.R")
source("~/codeThese/algosto/Outliers.R")
source("~/codeThese/algosto/computeOutliers.R")
source("~/codeThese/algosto/seuils.R")
for (i in 1:nbruns)
{
  resultats_outliers <- calcule_outliers(contamin = "moyenne")
}




taux_contamination <- c(0,2, 5, 10, 15, 20, 25, 30, 40)

###Test fonction d'estimation online affichage des boxplots des erreurs
pdf("Resultats_Erreurs_Sigma Toeplitzvar1sqrtd.pdf", width = 10, height = 7)

erreursSigmaBoxplot <- matrix(0,length(taux_contamination),n)
depart = 100
for (i in seq_along(taux_contamination)) 
{
  contamin = "moyenne"
  
  delta <- taux_contamination[i]
  #delta <- 0
  
  #Génération des échantillons
  p1 <- 1 - delta / 100
  
  p2 <- 1 - p1
  
  resultsSimul <- genererEchantillon(n, d, mu1, mu2, p1, p2, Sigma1 = Sigma1, Sigma2 = Sigma2,contamin = contamin)
  Z <- resultsSimul$Z
  
  
  #Estimation online
  #results <- estimMV(Z,Vinit = Sigma1,methode = "eigen")
  
  resultatsOnline <- detection(Z,depart, methodeEstimation = "online")
  SigmaEstimOnline <- resultatsOnline$SigmaOnline[(nrow(Z) - 1),,]
  distances <- resultatsOnline$distances
  #Estimation offline
  resultsOffline <- detection(Z,depart, methodeEstimation = "offline")
  SigmaEstimOffline <- resultsOffline$SigmaOffline
  
  erreursonline <- calculErreursNormeFrobenius(SigmaEstimOnline,Sigma1)
  erreursSigmaBoxplot[i,] <- erreursonline
  #erreursoffline <- calculErreursNormeFrobenius(SigmaEstimOffline,Sigma1)
 #affiche_erreursSigma(erreurs_online = erreursonline, contamination = delta)

 #Affichage des SigmaEstim Online, Offline et théorique
 
plot_comparaison_sigma(Sigma1, SigmaEstimOffline = SigmaEstimOffline,SigmaEstimOnline = SigmaEstimOnline,delta)

 
 
 #points(Sigma1, SigmaEstimOffline, col=4); abline(0, 1)
 
  
}

creer_boxplot_erreurs(erreursSigmaBoxplot,taux_contamination,methode= "online")
dev.off()

###Calcul du meilleur AUC pour chaque méthode

methodes <- c("Cov Empirique", "OGK", "Online", "Offline", "Comédiane", "Shrinkage")

taux_contamination <- c(2, 5, 10, 15, 20, 25, 30, 40)

#Création d'un dataframe pour contenir les meilleurs AUC

auc_df <- data.frame(matrix(ncol = length(methodes), nrow = length(taux_contamination)))
seuil_df <- data.frame(matrix(ncol = length(methodes), nrow = length(taux_contamination)))

rownames(auc_df) <- taux_contamination


rownames(seuil_df) <- taux_contamination

colnames(auc_df) <- methodes


colnames(seuil_df) <- methodes


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
  
  #Calcul du meilleur seuil selon l'AUC et de l'AUC selon chaque méthode
  for (m in methodes)
  {
      
         
      distances <- calcule_distances_par_methode(Z,methode = m)
      #print(distances)
      resultats <- courbeROC(resultsSimul$labelsVrais,distances)
      auc_df[i,m] <- resultats$auc_max
      seuil_df[i,m] <- resultats$seuil_auc_max
  }
  
}  
  #Enregistrement dans un CSV
  write.csv(auc_df,file = "auc.csv")

  write.csv(seuil_df,file = "seuil.csv")  

