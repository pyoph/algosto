#install.packages("easystats")
#install.packages("bigutilsr")
#install.packages("xtable")
#install.packages("RGMM")
#install.packages("rrcov")
#install.packages("mvnfast")
#install.packages("reshape")
#install.packages("matlib")
#install.packages("RobRegression")
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
library("RobRegression")
#library("rJava")
#library("REPPlab")
library("RGMM")
library("mvnfast")
library(ROCR)
library(pROC)
source("~/work/algosto/parametres.R")
source("~/work/algosto/simulations.R")
source("~/work/algosto/algorithmes.R")
source("~/work/algosto/resultats.R")
source("~/work/algosto/Outliers.R")
source("~/work/algosto/computeOutliers.R")
source("~/work/algosto/seuils.R")

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

save(moyenne_resultats,"outliersUnif.Rdata")


taux_contamination <- c(0,2, 5, 10, 15, 20, 25, 30, 40)

###Test fonction d'estimation online affichage des boxplots des erreurs
pdf("Resultats_Erreurs_Sigma Toeplitzvar1sqrtd.pdf", width = 10, height = 7)

erreursSigmaBoxplot <- matrix(0,length(taux_contamination),n)


#Somme des erreurs online et streaming pour les moyenner ensuite

somme_erreursOnline <- matrix(0,n,length(taux_contamination))
somme_erreursStreaming <- matrix(0,n,length(taux_contamination))


depart = 100
for (i in seq_along(taux_contamination)) 
{
  
  for (m in nbruns){
  contamin = "student"
  
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
  SigmaEstimOnline <- resultatsOnline$SigmaOnline
  distances <- resultatsOnline$distances
  #Estimation offline
  resultsOffline <- detection(Z,depart, methodeEstimation = "offline")
  SigmaEstimOffline <- resultsOffline$SigmaOffline
  
  #Estimation en streaming 
  resultsStr <- detection(Z, depart, methodeEstimation = "streaming")
  SigmaStr <- resultsStr$SigmaStreaming
  
  
  #calcul des erreurs à chaque itération
  for (k in (1: n)){
  erreursonline[k] <- calculErreursNormeFrobenius(SigmaEstimOnline[k,,],Sigma1)
  
  erreursStr[k] <- calculErreursNormeFrobenius(SigmaStr[k,,],Sigma1)
  }
  erreursSigmaBoxplot[i,] <- erreursonline
  #erreursoffline <- calculErreursNormeFrobenius(SigmaEstimOffline,Sigma1)
  #affiche_erreursSigma(erreurs_online = erreursonline, contamination = delta)
  
  #affiche_erreursSigma(erreurs_online = erreursStr, contamination = delta)
  #Stockage des erreurs :
  somme_erreursOnline[,i] <- somme_erreursOnline[,i] + erreursonline
  somme_erreursStreaming[,i] <- somme_erreursStreaming[,i] + erreursStr
 #Affichage des SigmaEstim Online, Offline et théorique
 
#plot_comparaison_sigma(Sigma1, SigmaEstimOffline = SigmaEstimOffline,SigmaEstimOnline = SigmaEstimOnline,delta)

 
  }
 #points(Sigma1, SigmaEstimOffline, col=4); abline(0, 1)
 
  
}

creer_boxplot_erreurs(erreursSigmaBoxplot,taux_contamination,methode= "online")
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
