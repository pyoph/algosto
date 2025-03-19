#install.packages("easystats")
#install.packages("bigutilsr")
#install.packages("xtable")
#install.packages("RGMM")
#install.packages("rrcov")
#install.packages("mvnfast")
#install.packages("reshape")
#install.packages("matlib")
#install.packages("RobRegression")
#install.packages("knitr")
#install.packages("Rcpp")
#install.packages("RcppArmadillo")
#install.packages("Rcpp")
#install.packages("RcppEigen")
#library(Rcpp)
#library(RcppEigen)
#install.packages("mvoutlier")
#install.packages("truncdist")  # À installer si nécessaire
#library(truncdist)

library(xtable)
library(reshape2)
library(RobRegression)
library("robustbase")
library(mvtnorm)
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
library("RobRegression")
#library("rJava")
#library("REPPlab")
library("RGMM")
library("mvnfast")
library(ROCR)
library(pROC)
library(tidyr)
library(dplyr)
library(knitr)
library("mvoutlier")
#library(Rcpp)
#library(RcppArmadillo)
#install.packages("Metrics")
#library(Metrics)
source("~/algosto/parametres.R")
source("~/algosto/simulations.R")
source("~/algosto/algorithmes.R")
source("~/algosto/resultats.R")
source("~/algosto/Outliers.R")
source("~/algosto/computeOutliers.R")
source("~/algosto/seuils.R")
source("~/algosto/shrinkageCabana.R")
#sourceCpp("~/algosto/valeursVecteursPropres.cpp")


###Tests CPP
#Tester sur batch = sqrt(d)
data <- genererEchantillon(n,d,mu1,mu2,0.7,0.3,Sigma1,Sigma2,"moyenne")

Z <- data$Z

respcout <- pcout(data$Z,makeplot=FALSE)

outl <- respcout$wfinal01

outl <- (outl + 1) %% 2


table(data$labelsVrais,outl)

temps_execution <- system.time({
  results <- estimMVOnline(data$Z,methode="CPP")
})

print(temps_execution)

norm(results$Sigma[1e4-1,,]-Sigma1,"F")
temps_execution <- system.time({
  results <- estimation(data$Z,methodeOnline = "CPP",methodeEstimation = "online")
})

print(temps_execution)

# Afficher les valeurs propres
print("Valeurs propres :")
print(resultats$valeurs_propres)

# Afficher les vecteurs propres
print("Vecteurs propres :")
print(resultats$vecteurs_propres)

# Afficher les vecteurs propres
print("Vecteurs propres :")
print(resultats$vecteurs_propres)

#######Calcul RMSE AUC et FP######

results_metrics <- RMSEAUCFPdataset(contamin = "UniformeTronquee")

results_metrics <- round(results_metrics,2)

# Définir les valeurs en abscisse
taux_contamination <- c(0, 2, 5, 10, 15, 20, 25, 30, 40)

# Sélectionner uniquement les colonnes RMSE_Sigma_*
rmse_df <- results_metrics[, c("RMSE_Sigma_Cov", "RMSE_Sigma_OGK", "RMSE_Sigma_Comed", 
                               "RMSE_Sigma_Shrink", "RMSE_Sigma_Online", 
                               "RMSE_Sigma_Offline", "RMSE_Sigma_Streaming","RMSE_Sigma_fastmcd")]

# Ajouter la colonne des taux de contamination
rmse_df$taux_contamination <- taux_contamination  

# Conversion en format long pour ggplot
df_long <- melt(rmse_df, id.vars = "taux_contamination", 
                variable.name = "Méthode", value.name = "RMSE")

# Supprimer le préfixe "RMSE_Sigma_" pour un affichage plus clair dans la légende
df_long$Méthode <- gsub("RMSE_Sigma_", "", df_long$Méthode)

# Tracer les courbes RMSE
ggplot(df_long, aes(x = taux_contamination, y = RMSE, color = Méthode)) +
  geom_line(size = 1.2) +  # Courbes
  geom_point(size = 2) +   # Points aux taux spécifiés
  scale_x_continuous(breaks = taux_contamination) +  # Spécifier les valeurs en abscisse
  scale_y_log10() +  # Échelle logarithmique en base 10
    labs(title = "Évolution du RMSE en fonction du taux de contamination",
       x = "Taux de contamination (%)",
       y = "RMSE",
       color = "Méthode") +  # Nom de la légende
  theme_minimal() +
  theme(legend.position = "top")  # Légende en haut


auc_df <- results_metrics[, c("AUC_Cov", "AUC_OGK", "AUC_Comed", 
                               "AUC_Shrink", "AUC_Online", 
                               "AUC_Offline", "AUC_Streaming")]

auc_df$taux_contamination <- taux_contamination  

# Conversion en format long pour ggplot
df_long <- melt(auc_df, id.vars = "taux_contamination", 
                variable.name = "Méthode", value.name = "AUC")

# Supprimer le préfixe "RMSE_Sigma_" pour un affichage plus clair dans la légende
df_long$Méthode <- gsub("AUC_", "", df_long$Méthode)

# Tracer les courbes RMSE
ggplot(df_long, aes(x = taux_contamination, y = AUC, color = Méthode)) +
  geom_line(size = 1.2) +  # Courbes
  geom_point(size = 2) +   # Points aux taux spécifiés
  scale_x_continuous(breaks = taux_contamination) +  # Spécifier les valeurs en abscisse
  labs(title = "Évolution de l'AUC en fonction du taux de contamination",
       x = "Taux de contamination (%)",
       y = "AUC",
       color = "Méthode") +  # Nom de la légende
  theme_minimal() +
  theme(legend.position = "top")  # Légende en haut

#faux positifs

fp_df <- results_metrics[, c("FP_Cov", "FP_OGK", "FP_Comed", 
                              "FP_Shrink", "FP_Online", 
                              "FP_Offline", "FP_Streaming")]

fp_df$taux_contamination <- taux_contamination  

# Conversion en format long pour ggplot
df_long <- melt(fp_df, id.vars = "taux_contamination", 
                variable.name = "Méthode", value.name = "FP")

# Supprimer le préfixe "RMSE_Sigma_" pour un affichage plus clair dans la légende
df_long$Méthode <- gsub("FP_", "", df_long$Méthode)

# Tracer les courbes FP
ggplot(df_long, aes(x = taux_contamination, y = FP, color = Méthode)) +
  geom_line(size = 1.2) +  # Courbes
  geom_point(size = 2) +   # Points aux taux spécifiés
  scale_x_continuous(breaks = taux_contamination) +  # Spécifier les valeurs en abscisse
  labs(title = "Évolution des faux positifs en fonction du taux de contamination",
       x = "Taux de contamination (%)",
       y = "FP",
       color = "Méthode") +  # Nom de la légende
  theme_minimal() +
  theme(legend.position = "top")  # Légende en haut




save(results_metrics,file = "results_metricsMaronnaZamar.Rdata")

# latex_table_results_metrics <- xtable(results_metrics)

results_without_RMSE_Med <- results_metrics %>%
     select(-contains("RMSE_Med")) %>%  # Supprimer toutes les colonnes RMSE_Med
     select(-contains("Cov"))

latex_table <- kable(results_without_RMSE_Med, format = "latex", caption = "Résultats contamination en moyenne 20 runs", label = "tab:results")
writeLines(latex_table, "resultats_contamination_zero.tex")

latex_table_results_metrics_without_RMSE_Med <- xtable(results_without_RMSE_Med)


  
#Création d'une liste vide 
liste_resultats_outliers <- list()

for (i in 1:nbruns)
{
  resultats_outliers <- calcule_outliers(contamin = "uniforme")
  liste_resultats_outliers[[i]] <- resultats_outliers 
}

#resultats_outliers <- round(resultats_outliers,2)


# Calculer la moyenne de chaque colonne sur tous les dataframes
moyenne_resultats <- round(Reduce("+", liste_resultats_outliers) / nbruns,2)

# Afficher la moyenne
print(moyenne_resultats)

save(moyenne_resultats,file = "outliersUnifhypercubeunite.Rdata")

moyenne_resultatsEnreg <- moyenne_resultats  %>% select(-FN_Cov, -FP_Cov, -Tps_Cov)

latex_table <- xtable(moyenne_resultatsEnreg)


taux_contamination <- c(0,2, 5, 10, 15, 20, 25, 30, 40)

###Test fonction d'estimation online affichage des boxplots des erreurs
pdf("Resultats_Erreurs_Sigma Toeplitzvar1sqrtd.pdf", width = 10, height = 7)


#Calcul des erreurs pour tracer les boxplots

erreursSigmaBoxplotOnline <-array(0, dim = c(nbruns,length(taux_contamination),n))

erreursSigmaBoxplotStreaming <- array(0, dim = c(nbruns,length(taux_contamination),n))

erreursVarOracle <- matrix(0,nbruns,length(taux_contamination))

#Somme des erreurs online et streaming pour les moyenner ensuite

somme_erreursOnline <- matrix(0,n,length(taux_contamination))
somme_erreursStreaming <- matrix(0,n,length(taux_contamination))


depart = 100

for (m in (1:nbruns)){
for (i in seq_along(taux_contamination)) 
{
  

  contamin = "moyenne"
  #Initialisation des erreurs online et streaming
  erreursonline <- rep(0,n)
  erreursStr <- rep(0,n)
  delta <- taux_contamination[i]
  #delta <- 0
  
  #Génération des échantillons
  p1 <- 1 - delta / 100
  
  p2 <- 1 - p1
  
  resultsSimul <- genererEchantillon(n, d, mu1, mu2, p1, p2, Sigma1 = Sigma1, Sigma2 = Sigma2,contamin = contamin)
  Z <- resultsSimul$Z
  
  
  #Estimation online
  #results <- estimMV(Z,Vinit = Sigma1,methode = "eigen")
  
  resultatsOnline <- estimation(Z,depart, methodeEstimation = "online")
  SigmaEstimOnline <- resultatsOnline$SigmaOnline
  distances <- resultatsOnline$distances
  miterOn <- resultatsOnline$miter
  
  #Estimation offline
  resultsOffline <- estimation(Z,depart, methodeEstimation = "offline")
  SigmaEstimOffline <- resultsOffline$SigmaOffline
  
  #Estimation en streaming 
  resultsStr <- estimation(Z, depart, methodeEstimation = "streaming")
  miterStr <- estimation(Z, depart, methodeEstimation = "streaming")$miter
  SigmaStr <- resultsStr$SigmaStreaming
  
  
  #calcul des erreurs à chaque itération
  for (k in (1: n)){
  erreursonline[k] <- calculErreursNormeFrobenius(SigmaEstimOnline[k,,],Sigma1)
  
  erreursStr[k] <- calculErreursNormeFrobenius(SigmaStr[k,,],Sigma1)
  }
  erreursSigmaBoxplotOnline[m,i,] <- erreursonline
  erreursSigmaBoxplotStreaming[m,i,] <- erreursStr
  erreursVarOracle[m,i] <- EstimVarMC(nbiter = 1e2, delta = delta, Sigma = Sigma1)
  
  #erreursoffline <- calculErreursNormeFrobenius(SigmaEstimOffline,Sigma1)
  #affiche_erreursSigma(erreurs_online = erreursonline, contamination = delta)
  
  #affiche_erreursSigma(erreurs_online = erreursStr, contamination = delta)
  #Stockage des erreurs :
  print(somme_erreursOnline[1e4,i])
  somme_erreursOnline[,i] <- somme_erreursOnline[,i] + erreursonline
  print(somme_erreursOnline[1e4,i])
  somme_erreursStreaming[,i] <- somme_erreursStreaming[,i] + erreursStr
 #Affichage des SigmaEstim Online, Offline et théorique
 
#plot_comparaison_sigma(Sigma1, SigmaEstimOffline = SigmaEstimOffline,SigmaEstimOnline = SigmaEstimOnline,delta)

 
  }
 #points(Sigma1, SigmaEstimOffline, col=4); abline(0, 1)
 
  
}
#Calcul des moyennes de erreurs

moy_erreursSigmaBoxplotOnline <- round(apply(erreursSigmaBoxplotOnline, c(2,3), mean), 2)
moy_erreursSigmaBoxplotStreaming <- round(apply(erreursSigmaBoxplotStreaming, c(2,3), mean), 2)

moy_erreursVarOracle <- round(apply(erreursVarOracle, 2, mean), 2)

creer_boxplot_erreurs(moy_erreursSigmaBoxplotOnline,taux_contamination,methode= "online",erreursVarOracle[1,])
creer_boxplot_erreurs(moy_erreursSigmaBoxplotStreaming,taux_contamination,methode= "streaming",erreursVarOracle[1,])


dev.off()

#Calcul des moyennes des erreurs 

# Calcul de la moyenne des erreurs par taux de contamination
moy_erreursOnline <- (somme_erreursOnline[100:1e4,])/nbruns
moy_erreursStreaming <-(somme_erreursStreaming[100:1e4,])/nbruns

# Créer une liste pour stocker les graphiques
list_plots <- list()
# Nombre d'itérations à considérer
iterations <- seq(1, 9900)  
for (i in seq_along(taux_contamination)) 
{
  delta <- taux_contamination[i]
# Transformer les données en format long pour ggplot
df_erreurs <- data.frame(
  Iteration = rep(iterations, 2),
  Erreur_Moyenne = c(moy_erreursOnline[iterations,i], moy_erreursStreaming[iterations,i]),
  Méthode = rep(c("Online", "Streaming"), each = length(iterations))
)

# Tracer les courbes
p <- ggplot(df_erreurs, aes(x = Iteration, y = Erreur_Moyenne, color = Méthode)) +
  geom_line(size = 1.2) +  # Ligne plus épaisse
  geom_point(size = 2) +   # Points pour les valeurs précises
  scale_x_log10() +        # Échelle logarithmique pour les itérations
  labs(title = paste("taux de contamination =  ",delta," %"),
       x = "Nombre d'itérations",
       y = "Erreur Moyenne (norme de Frobenius)") +
  theme_minimal() +
  theme(legend.position = "top")
print(p)
list_plots[[i]] <- p
}

grid.arrange(grobs = list_plots[1:4], ncol = 2, nrow = 2)
grid.arrange(grobs = list_plots[5:9], ncol = 3, nrow = 3)




# 
# ###Calcul du meilleur AUC pour chaque méthode
# 
# methodes <- c("Cov Empirique", "OGK", "Online", "Offline", "Comédiane", "Shrinkage")
# 
# taux_contamination <- c(2, 5, 10, 15, 20, 25, 30, 40)
# 
# #Création d'un dataframe pour contenir les meilleurs AUC
# 
# auc_df <- data.frame(matrix(ncol = length(methodes), nrow = length(taux_contamination)))
# seuil_df <- data.frame(matrix(ncol = length(methodes), nrow = length(taux_contamination)))
# 
# rownames(auc_df) <- taux_contamination
# 
# 
# rownames(seuil_df) <- taux_contamination
# 
# colnames(auc_df) <- methodes
# 
# 
# colnames(seuil_df) <- methodes
# 
# 
# for (i in seq_along(taux_contamination)) 
# {
#   contamin = "moyenne"
#   
#   delta <- taux_contamination[i]
#   #delta <- 2
#   
#   #Génération des échantillons
#   p1 <- 1 - delta / 100
#   
#   p2 <- 1 - p1
#   
#   resultsSimul <- genererEchantillon(n, d, mu1, mu2, p1, p2, Sigma1 = Sigma1, Sigma2 = Sigma2,contamin = contamin)
#   Z <- resultsSimul$Z
#   
#   #Calcul du meilleur seuil selon l'AUC et de l'AUC selon chaque méthode
#   for (m in methodes)
#   {
#       
#          
#       distances <- calcule_distances_par_methode(Z,methode = m)
#       #print(distances)
#       resultats <- courbeROC(resultsSimul$labelsVrais,distances)
#       auc_df[i,m] <- resultats$auc_max
#       seuil_df[i,m] <- resultats$seuil_auc_max
#   }
#   
# }  
#   #Enregistrement dans un CSV
#   write.csv(auc_df,file = "auc.csv")
# 
#   write.csv(seuil_df,file = "seuil.csv")  
# 
