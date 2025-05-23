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
#library(RobRegression)
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
#library("easystats")
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
library(Rcpp)
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
#sourceCpp("~/algosto/RobinsMC2CPP.cpp")
#sourceCpp("~/algosto/valeursVecteursPropres.cpp")


p1 = 0.6
data <- genererEchantillon(n,d,mu1,mu2,p1,1-p1,Sigma1,Sigma2,"moyenne")

Z = data$Z

cumulativeOutlierDetection <- function(resultats, data, pourcentage,titre) {
  total_points <- nrow(data$Z)
  total_outliers_theoriques <- pourcentage / 100 * total_points
  
  nb_outliers_detectes <- 0
  nb_outliers_detectes_vrais <- 0
  nb_outliers_vrais <- 0
  
  taux_outliers_detectes <- numeric(total_points)
  taux_outliers_vrais <- numeric(total_points)
  taux_outliers_detectes_vrais = numeric(total_points)
  
  for (i in 1:total_points) {
    if(data$labelsVrais[i] == 1) nb_outliers_vrais =  nb_outliers_vrais + 1
    
    if (data$labelsVrais[i] == 1 && resultats$outlier_labels[i] == 1) {
      nb_outliers_detectes <- nb_outliers_detectes + 1
      #nb_outliers_vrais <- nb_outliers_vrais + 1
      nb_outliers_detectes_vrais = nb_outliers_detectes_vrais + 1
      #taux_outliers_vrais[i] <- nb_outliers_detectes_vrais
    } else if (data$labelsVrais[i] == 0 && resultats$outlier_labels[i] == 1) {
      nb_outliers_detectes <- nb_outliers_detectes + 1}
      if(nb_outliers_vrais != 0){
        taux_outliers_detectes[i] <- nb_outliers_detectes / nb_outliers_vrais * 100
        #taux_outliers_detectes[i] <- nb_outliers_detectes
      }
      else {taux_outliers_detectes[i] <-  100
        #taux_outliers_detectes[i] <- nb_outliers_detectes
      }
    if(nb_outliers_vrais != 0){
      taux_outliers_detectes_vrais[i] <- nb_outliers_detectes_vrais / nb_outliers_vrais * 100
      #taux_outliers_detectes[i] <- nb_outliers_detectes
    }
    else {taux_outliers_detectes_vrais[i] <-  100
    #taux_outliers_detectes[i] <- nb_outliers_detectes
    }
    
    
      
      taux_outliers_vrais[i] <- 100
      
      print(paste("nb_outliers_detectes_vrais ",nb_outliers_detectes_vrais ))
    print(paste("nb_outliers_vrais ",nb_outliers_vrais ))
    print(paste("taux_outliers_vrais[i] ",taux_outliers_vrais[i] ))
    }
    
   

  
  df <- data.frame(
    index = 1:total_points,
    Detected_rate = taux_outliers_detectes,
    True_outliers = taux_outliers_vrais,
    True_positive_rate = taux_outliers_detectes_vrais
  )
  
  p <- ggplot(df, aes(x = index)) +
    geom_line(aes(y = Detected_rate, color = "True and false positive rate"), size = 1.2) +
    geom_line(aes(y = True_positive_rate, color = "True positives rate"), size = 1.2) +
    geom_line(aes(y = True_outliers, color = "True outliers rate"), size = 1.2) +
    scale_color_manual(values = c(
      "True outliers rate" = "red",
      "True and false positive rate" = "orange",
      "True positives rate" = "purple"
    )) +
    #scale_x_log10() +
    labs(
      title = paste(titre, "-", pourcentage, "% of outliers"),
      x = "Data index",
      y = "Cumulative rate (%)",
      color = "Legend"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  print(taux_outliers_vrais[1:10])
  return(list(p = p,taux_outliers_vrais = taux_outliers_vrais,taux_outliers_detectes = taux_outliers_detectes))
}
plot_obj <- cumulativeOutlierDetection(resultats, data, 20, "Shifted gaussian")$p

plot_objMarZam = cumulativeOutlierDetection(resultats, data, 20,"Maronna Zamar")$p

plot_objTrSt = cumulativeOutlierDetection(resultats, data, 20,"Truncated Student")$p

grid.arrange(plot_obj,plot_objTrSt,plot_objMarZam,nrow = 2,ncol = 2)

par(mfrow = c(2, 2))  # 2 lignes, 2 colonnes
plot_objTrSt
plot_obj
plot_objMarZam
plot.new()

title("Comparaison des scénarios", line = -2, cex.main = 1.4)

# Fermer le fichier pour sauvegarder
dev.off()

temps_calcul = results_metrics$temps_calculBP

methods <- c("Cov", "OGK", "Comed", "Shrink", 
             "Offline", "Online", "Streaming", "FastMCD")

# Aplatir le tableaugra
df <- melt(temps_calcul)
colnames(df) <- c("i", "j", "k", "temps")

# Remplacer la colonne 'k' par le nom de méthode correspondant
df$method <- factor(df$k, labels = methods)

# Enlever les méthodes "Shrink" et "Offline"
df <- subset(df, !method %in% c("Shrink", "Offline"))

# Tracer le boxplot
ggplot(df, aes(x = method, y = temps)) +
  geom_boxplot(fill = "lightblue") +
  scale_y_log10() +
  labs(title = "Boxplot of computation times by method",
       x = "Method",
       y = "Computation time") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


Z <- data$Z

respcout <- pcout(data$Z,makeplot=FALSE)

outl <- respcout$wfinal01

outl <- (outl + 1) %% 2


table(data$labelsVrais,outl)

#######Calcul RMSE AUC et FP######


#results_metrics <- RMSEAUCFPdataset(contamin = "studentTronquee")$result_metrics
#results_metrics <- RMSEAUCFPdataset(contamin = "UniformeTronquee")$result_metrics

results_metrics <- round(results_metrics,2)

# Définir les valeurs en abscisse
taux_contamination <- c(0, 2, 5, 10, 15, 20, 25, 30, 40)

# Sélectionner uniquement les colonnes RMSE_Sigma_*
rmse_df <- rm[, c("RMSE_Sigma_Cov", "RMSE_Sigma_OGK", "RMSE_Sigma_Comed", 
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
plot1 = ggplot(df_long, aes(x = taux_contamination, y = RMSE, color = Méthode)) +
  geom_line(size = 1.2) +  # Courbes
  geom_point(size = 2) +   # Points aux taux spécifiés
  scale_x_continuous(breaks = taux_contamination) +  # Spécifier les valeurs en abscisse
  scale_y_log10() +  # Échelle logarithmique en base 10
    labs(title = "Shifted Gaussian contamination scenario",
       x = "Contamination rate (%)",
       y = "Frobenius norm error",
       color = "Method") +  # Nom de la légende
  theme_minimal() +
  theme(legend.position = "top")  # Légende en haut



# Tracer les courbes RMSE
plot2 = ggplot(df_long, aes(x = taux_contamination, y = RMSE, color = Méthode)) +
  geom_line(size = 1.2) +  # Courbes
  geom_point(size = 2) +   # Points aux taux spécifiés
  scale_x_continuous(breaks = taux_contamination) +  # Spécifier les valeurs en abscisse
  scale_y_log10() +  # Échelle logarithmique en base 10
  labs(title = "Truncated Student contamination scenario",
       x = "Contamination rate (%)",
       y = "Frobenius norm error",
       color = "Method") +  # Nom de la légende
  theme_minimal() +
  theme(legend.position = "top")  # Légende en haut

# Tracer les courbes RMSE
plot3 = ggplot(df_long, aes(x = taux_contamination, y = RMSE, color = Méthode)) +
  geom_line(size = 1.2) +  # Courbes
  geom_point(size = 2) +   # Points aux taux spécifiés
  scale_x_continuous(breaks = taux_contamination) +  # Spécifier les valeurs en abscisse
  scale_y_log10() +  # Échelle logarithmique en base 10
  labs(title = "Maronna Zamar contamination scenario",
       x = "Contamination rate (%)",
       y = "Frobenius norm error",
       color = "Method") +  # Nom de la légende
  theme_minimal() +
  theme(legend.position = "top")  # Légende en haut

grid.arrange(plot2, plot1, plot3, ncol = 2, nrow = 2)

auc_df <- rm[, c("AUC_Cov", "AUC_OGK", "AUC_Comed", 
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
       x = "Contamnination rate (%)",
       y = "AUC",
       color = "Method") +  # Nom de la légende
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


# Sélectionner les colonnes FP et FN
fp_fn_df <- results_metrics[, c("FP_Cov", "FP_OGK", "FP_Comed", 
                                "FP_Shrink", "FP_Online", 
                                "FP_Offline", "FP_Streaming",
                                "FN_Cov", "FN_OGK", "FN_Comed", 
                                "FN_Shrink", "FN_Online", 
                                "FN_Offline", "FN_Streaming")]

fp_fn_df$taux_contamination <- taux_contamination  

# Conversion en format long pour ggplot
df_long <- melt(fp_fn_df, id.vars = "taux_contamination", 
                variable.name = "Méthode", value.name = "Valeur")

# Ajouter une colonne pour distinguer FP et FN
df_long$Type <- ifelse(grepl("FP_", df_long$Méthode), "FP", "FN")

# Nettoyer les noms des méthodes en supprimant le préfixe FP_ ou FN_
df_long$Méthode <- gsub("FP_|FN_", "", df_long$Méthode)

# Tracer les courbes FP et FN
plot2 = ggplot(df_long, aes(x = taux_contamination, y = Valeur, 
                    color = Méthode, linetype = Type)) +
  geom_line(size = 1.2) +  # Courbes
  geom_point(size = 2) +   # Points aux taux spécifiés
  scale_x_continuous(breaks = taux_contamination) +  # Spécifier les valeurs en abscisse
  labs(title = "False positives and false negatives rates",
       x = "Contamination rate (%)",
       y = "Rate",
       color = "Method",
       linetype = "Type") +  # Nom de la légende
  theme_minimal() +
  theme(legend.position = "top")  # Légende en haut

plot1+plot2

# Initialisation

# Dimensions : [100, 9, 8]
rmse_data <- data.frame()

rmseSigmaBP = results_metrics$rmseSigmaBP

# Noms des méthodes
methods <- c("Cov", "OGK", "Comed", "Shrink", 
             "Online", "Offline", "Streaming", "FastMCD")

# Noms des taux de contamination 
#taux_contamination <- paste0("Taux_", 1:9)
taux_contamination <- c(0, 2, 5, 10, 15, 20, 25, 30, 40)

for (m in 1:8) {
  for (c in 1:9) {
    values <- rmseSigmaBP[, c, m]  # 100 valeurs
    temp_df <- data.frame(
      RMSE = values,
      Méthode = methods[m],
      Taux = taux_contamination[c]
    )
    rmse_data <- rbind(rmse_data, temp_df)
  }
}
ggplot(rmse_data, aes(x = factor(Taux), y = RMSE)) +
  geom_boxplot(fill = "lightblue", color = "darkblue") +
  facet_wrap(~ Méthode, scales = "fixed") +  # ou supprime scales complètement
  theme_minimal(base_size = 14) +
  labs(
    title = "Boxplot of Frobenius norm error",
    x = "Contamination rate",
    y = "Frobenius norm"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


####Représentation des contaminations

# Fonction pour générer un échantillon
preparer_donnees <- function(type,p1 = 0.8) {
  set.seed(123)
  p1 = p1
  resultat <- genererEchantillon(n,d,mu1,mu2,p1,1-p1,Sigma1,Sigma2,type)
  df <- as.data.frame(resultat$Z)
  df$label <- factor(resultat$labelsVrais)
  return(df)
}

# Génération des trois jeux de données
donnees_moyenne  <- preparer_donnees("moyenne",p1 = 0.8)
donnees_student  <- preparer_donnees("studentTronquee",p1 = 0.8)
donnees_maronna  <- preparer_donnees("MaronnaZamar",p1 =0.8)

# Fonction de tracé avec échelle log10 sur Y uniquement (données inchangées)
tracer_plot <- function(df, titre) {
  df_inliers <- df[df$label == 0, ]  # Sélection des inliers
  
  ggplot(df, aes(x = V1, y = V2, color = label)) +
    geom_point(alpha = 0.7) +
    #stat_ellipse(data = df_inliers, aes(x = V1, y = V2),
     #            level = 0.95, size = 1, linetype = "dashed", inherit.aes = FALSE, color = "black") +    
    #scale_x_log10() +
    #scale_y_log10() +  # Seulement l’échelle Y en log10 (pas les données)
    theme_minimal() +
    labs(title = titre,
         x = "First composant",
         y =  "Second composant",
         color = "Group") +
    xlim(-10,10) +
    ylim(-10,10)+
    scale_color_manual(values = c("blue", "red"),
                       labels = c("Clean", "Contaminated"))
}

# Graphiques
graph_moyenne  <- tracer_plot(donnees_moyenne,  "Shifted Gaussian")
graph_student  <- tracer_plot(donnees_student,  "Truncated Student")
graph_maronna  <- tracer_plot(donnees_maronna,  "Maronna-Zamar")

# Affichage
grid.arrange(graph_moyenne, graph_student, graph_maronna, ncol = 3)


# Paramètres
alpha <- 0.05
conf_level <- 1 - alpha
p_seq <- seq(0.6, 1, length.out = 200)
N <- nrow(Z)  # N doit être défini au préalable

# Calcul des n0 et des x théoriques
n0 <- floor(p_seq * N)
x <- round(n0 * alpha)

# Librairie pour les intervalles de confiance
library(binom)

# Calcul de l’intervalle de confiance exact (Clopper-Pearson)
conf <- binom.confint(x = x, n = n0, conf.level = conf_level, methods = "exact")

# Tracé des bornes
plot(p_seq, conf$lower, type = "l", col = "blue", lwd = 2,
     ylim = c(0, 0.15), xlab = "p (proportion pour n0 = p * N)",
     ylab = "Intervalle de confiance de la proportion",
     main = "IC à 95% pour Bin(n0, alpha = 0.05)")
lines(p_seq, conf$upper, col = "red", lwd = 2)

# Méthode Cov
rm2 <- rm$FP_Cov
fp_cov_values <- round(rm2 / 100, 4)
fp_cov_positions <- c(1.00, 0.98, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.60)
points(fp_cov_positions, fp_cov_values, col = "black", pch = 19)
text(fp_cov_positions, fp_cov_values + 0.005,
     labels = paste0(round(rm2, 2), "%"), cex = 0.7, pos = 3, col = "black")

# Méthode OGK
rm3 <- rm$FP_OGK
fp_ogk_values <- round(rm3 / 100, 4)
fp_ogk_positions <- fp_cov_positions  # même positions
points(fp_ogk_positions, fp_ogk_values, col = "green", pch = 19)
text(fp_ogk_positions, fp_ogk_values + 0.005,
     labels = paste0(round(rm3, 2), "%"), cex = 0.7, pos = 3, col = "green")

# Méthode Comed
rm4 <- rm$FP_Comed
fp_comed_values <- round(rm4 / 100, 4)
fp_comed_positions <- fp_cov_positions
points(fp_comed_positions, fp_comed_values, col = "gray", pch = 19)
text(fp_comed_positions, fp_comed_values + 0.005,
     labels = paste0(round(rm4, 2), "%"), cex = 0.7, pos = 3, col = "gray")

# Méthode Shrinkage
rm5 <- rm$FP_Shrink
fp_shrink_values <- round(rm5 / 100, 4)
fp_shrink_positions <- fp_cov_positions
points(fp_shrink_positions, fp_shrink_values, col = "orange", pch = 19)
text(fp_shrink_positions, fp_shrink_values + 0.005,
     labels = paste0(round(rm5, 2), "%"), cex = 0.7, pos = 3, col = "orange")


rm9 = rm$FP_fastmcd

fp_fastmcd_values <- round(rm9 / 100, 4)
fp_fastmcd_positions <- fp_cov_positions
points(fp_fastmcd_positions, fp_fastmcd_values, col = "lightblue", pch = 19)
text(fp_fastmcd_positions, fp_fastmcd_values + 0.005,
     labels = paste0(round(rm9, 2), "%"), cex = 0.7, pos = 3, col = "lightblue")


rm6 = round(rm$FP_Online,2)

fp_online_values <- round(rm6 / 100, 4)
fp_online_positions <- fp_cov_positions

points(fp_online_positions, fp_online_values, col = "darkblue", pch = 19)
text(fp_shrink_positions, fp_shrink_values + 0.005,
     labels = paste0(round(rm5, 2), "%"), cex = 0.7, pos = 3, col = "darkblue")



rm7 = round(rm$FP_Offline,2)

fp_offline_values <- round(rm7 / 100, 4)
fp_offline_positions <- fp_cov_positions

points(fp_offline_positions, fp_offline_values, col = "darkgreen", pch = 19)
text(fp_offline_positions, fp_offline_values + 0.005,
     labels = paste0(round(rm7, 2), "%"), cex = 0.7, pos = 3, col = "darkgreen")


rm8 = round(rm$FP_Streaming,2)


fp_streaming_values <- round(rm8 / 100, 4)
fp_streaming_positions <- fp_cov_positions

points(fp_streaming_positions, fp_streaming_values, col = "darkred", pch = 19)
text(fp_offline_positions, fp_offline_values + 0.005,
     labels = paste0(round(rm7, 2), "%"), cex = 0.7, pos = 3, col = "darkred")



# Légende unique
legend("topright",
       legend = c("Borne inférieure", "Borne supérieure",
                  "Cov", "OGK", "Comed", "Shrinkage","Online","Offline","Streaming"),
       col = c("blue", "red", "black", "green", "gray", "orange","lightblue","darkblue","darkgreen","darkred"),
       lwd = c(2, 2, NA, NA, NA, NA,NA,NA,NA),
       pch = c(NA, NA, 19, 19, 19, 19,19,19,19),
       cex = 0.35)  # Taille réduite



taux_contamination = c(0, 2, 5, 10, 15, 20, 25, 30, 40)

cutoff = qchisq(0.95,df = d)

fp_vrai = rep(0,length(taux_contamination))
fp_vrai_corr = rep(0,length(taux_contamination))


fp_offline = rep(0,length(taux_contamination))
fp_offline_corr = rep(0,length(taux_contamination))

fp_online = rep(0,length(taux_contamination))
fp_online_corr = rep(0,length(taux_contamination))

fp_streaming = rep(0,length(taux_contamination))
fp_streaming_corr = rep(0,length(taux_contamination))


 fp_cov = rep(0,length(taux_contamination))
fp_cov_corr = rep(0,length(taux_contamination))


compt = 1




methodes = c("sample_cov_online","offline","online","batch")
rmseSigma = matrix(0, length(taux_contamination),length(methodes))

for (r in taux_contamination){
compt_meth = 1  
data <- genererEchantillon(n,d,mu1,mu2,p1 = 1- r/100,r/100,Sigma1,Sigma2,"moyenne")

Z = data$Z


outliers = rep(0,nrow(Z))
for(i in (1:nrow(Z))){
  distances[i] = mahalanobis_generalizedRcpp(Z[i,],mu1,eigen(Sigma1)$vectors,eigen(Sigma1)$values)
  S = distances[i]
  
  if (S > cutoff) {outliers[i] = 1}
  
  
}

tc <- table(data$labelsVrais[1:(nrow(Z))], as.numeric(outliers)[1:(nrow(Z))])

tc <- safe_access_tc(tc)


if((tc["0","0"] + tc["0","1"]) != 0)
{fp_vrai[compt]   <- round((tc["0", "1"]/(tc["0", "1"] + tc["0", "0"])*100),2)}






resultats = update_mean_Sigma2(Z)

med <- resultats$mean
Sigma <- resultats$Sigma2

#rmseSigma[r,compt_meth] = norm(Sigma - Sigma1,"F")
compt_meth = compt_meth + 1

#distances <- calcule_vecteur_distances(Z, med, Sigma)
#eigen(Sigma1)$values

distances = rep(0, nrow(Z))
outliers = rep(0,nrow(Z))
for (i in (1:nrow(Z))){
distances[i] = mahalanobis_generalizedRcpp(Z[i,],med,eigen(Sigma)$vectors, eigen(Sigma)$values)
S = distances[i]

if (S > cutoff) {outliers[i] = 1}

}

tc <- table(data$labelsVrais[1:(nrow(Z))], as.numeric(outliers)[1:(nrow(Z))])

tc <- safe_access_tc(tc)

resultats = OfflineOutlierDetection(Z)
Sigma = resultats$variance + best_lambda*diag(d)


distances <- resultats$distances
outliers <- resultats$outlier_labels



tc <- table(data$labelsVrais[1:(nrow(Z))], as.numeric(outliers)[1:(nrow(Z))])

tc <- safe_access_tc(tc)
if((tc["0","0"] + tc["0","1"]) != 0)
{fp_offline[compt]   <-(tc["0", "1"]/(tc["0", "1"] + tc["0", "0"])*100)


best_lambda = tune_lambda(Z,mu_hat = resultats$median,Sigma = resultats$variance)$lambda

distances_corr = rep(0,nrow(Z))
outliers = rep(0,nrow(Z))

for(i in (1:nrow(Z))){
  distances_corr[i] = mahalanobis_generalizedRcpp(Z[i,],resultats$median,eigen(Sigma)$vectors,eigen(Sigma)$values)
  S = distances_corr[i]
  
  if (S > cutoff) {outliers[i] = 1}
  
  
}

tc <- table(data$labelsVrais[1:(nrow(Z))], as.numeric(outliers)[1:(nrow(Z))])
tc <- safe_access_tc(tc)
if((tc["0","0"] + tc["0","1"]) != 0)
{fp_offline_corr[compt]   <- (tc["0", "1"]/(tc["0", "1"] + tc["0", "0"])*100)
}
med <- colMeans(Z)
Sigma <- cov(Z)

distances = rep(0, nrow(Z))
outliers = rep(0,nrow(Z))
for (i in (1:nrow(Z))){
  distances[i] = mahalanobis_generalizedRcpp(Z[i,],med,eigen(Sigma)$vectors, eigen(Sigma)$values)
  S = distances[i]
  
  if (S > cutoff) {outliers[i] = 1}
  
}

tc <- table(data$labelsVrais[1:(nrow(Z))], as.numeric(outliers)[1:(nrow(Z))])

tc <- safe_access_tc(tc)

if((tc["0","0"] + tc["0","1"]) != 0)
{fp_cov[compt]   <- round((tc["0", "1"]/(tc["0", "1"] + tc["0", "0"]))*100,2)}
distances_corr = facteur_corr*distances
outliers <- detectionOutliers(distances_corr,cutoff = cutoff)
tc <- table(data$labelsVrais[1:(nrow(Z))], as.numeric(outliers)[1:(nrow(Z))])
tc
tc <- safe_access_tc(tc)
if((tc["0","0"] + tc["0","1"]) != 0)
{fp_cov_corr[compt]   <- round((tc["0", "1"]/(tc["0", "1"] + tc["0", "0"])*100),2)
}


resultats = StreamingOutlierDetection(Z,batch = 1,cutoff = cutoff)
distances <- resultats$distances
outliers = resultats$outlier_labels


tc <- table(data$labelsVrais[1:(nrow(Z))], as.numeric(outliers)[1:(nrow(Z))])
tc <- safe_access_tc(tc)
if((tc["0","0"] + tc["0","1"]) != 0)
{fp_online[compt]   <- round((tc["0", "1"]/(tc["0", "1"] + tc["0", "0"]))*100,2)}


best_lambda = tune_lambda(Z,resultats$moyennem,resultats$Sigma[nrow(Z),,])$lambda
# facteurs = correctionDistanceMahalanobis(distances,Z,methode = "online")
# distances_corr = hadamard.prod(facteurs,distances)
# outliers <- detectionOutliers(distances_corr,cutoff = qchisq(0.95,df = d))
#
Sigma = resultats$Sigma[nrow(Z) ,,] + best_lambda * diag(ncol(Z))
distances = rep(0, nrow(Z))
outliers = rep(0,nrow(Z))
for (i in (1:nrow(Z))){
  distances[i] = mahalanobis_generalizedRcpp(Z[i,],med,eigen(Sigma)$vectors, eigen(Sigma)$values)
  S = distances[i]
  
  if (S > cutoff) {outliers[i] = 1}
  
}

tc <- table(data$labelsVrais[1:(nrow(Z))], as.numeric(outliers)[1:(nrow(Z))])
tc <- safe_access_tc(tc)
if((tc["0","0"] + tc["0","1"]) != 0)
{fp_online_corr[compt]   <- (tc["0", "1"]/(tc["0", "1"] + tc["0", "0"]))*100}


resultats = StreamingOutlierDetection(Z,batch = ncol(Z))
distances <- resultats$distances
outliers = resultats$outlier_labels


tc <- table(data$labelsVrais[1:(nrow(Z))], as.numeric(outliers)[1:(nrow(Z))])
tc <- safe_access_tc(tc)
if((tc["0","0"] + tc["0","1"]) != 0)
{fp_streaming[compt]   <- (tc["0", "1"]/(tc["0", "1"] + tc["0", "0"]))*100}


best_lambda = tune_lambda(Z,resultats$moyennem,resultats$Sigma[nrow(Z),,])$lambda



# facteurs = correctionDistanceMahalanobis(distances,Z,methode = "online")
# distances_corr = hadamard.prod(facteurs,distances)
# outliers <- detectionOutliers(distances_corr,cutoff = qchisq(0.95,df = d))
#
print("ok online")
Sigma = resultats$Sigma[nrow(Z) ,,] + best_lambda * diag(d)
distances = rep(0, nrow(Z))
outliers = rep(0,nrow(Z))
for (i in (1:nrow(Z))){
  distances[i] = mahalanobis_generalizedRcpp(Z[i,],resultats$moyennem,eigen(Sigma)$vectors, eigen(Sigma)$values)
  S = distances[i]
  
  if (S > cutoff) {outliers[i] = 1}
  
}

# facteurs = correctionDistanceMahalanobis(distances,Z,methode = "streaming")
# 
# 
# 
# distances_corr = hadamard.prod(facteurs,distances)
# outliers <- detectionOutliers(distances_corr,cutoff = qchisq(0.95,df = d))

tc <- table(data$labelsVrais[1:(nrow(Z))], as.numeric(outliers)[1:(nrow(Z))])
tc <- safe_access_tc(tc)
if((tc["0","0"] + tc["0","1"]) != 0)
{fp_streaming_corr[compt]   <- (tc["0", "1"]/(tc["0", "1"] + tc["0", "0"]))*100}


print("ok streaming")

compt = compt + 1
}}

# Paramètres
alpha <- 0.05
conf_level <- 1 - alpha
p_seq <- seq(0.6, 1, length.out = 200)
N <- nrow(Z)  # N doit être défini au préalable

# Calcul des n0 et des x théoriques
n0 <- floor(p_seq * N)
x <- round(n0 * alpha)




# Calcul de l’intervalle de confiance exact (Clopper-Pearson)
conf <- binom.confint(x = x, n = n0, conf.level = conf_level, methods = "exact")

#Plot offline

# Tracé des bornes
plot(p_seq, conf$lower, type = "l", col = "blue", lwd = 2,
     ylim = c(0, 0.1), xlab = "p (proportion pour n0 = p * N)",
     ylab = "Intervalle de confiance de la proportion",
     main = "IC à 95% pour Bin(n0, alpha = 0.05) pour offline")
lines(p_seq, conf$upper, col = "red", lwd = 2)

fp_offline_positions = c(1.00, 0.98, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.60)

fp_offline_values = fp_offline/100


fp_offline_values_corr = round(fp_offline_corr/100,2)


points(fp_offline_positions, fp_offline_values_corr, col = "darkgreen", pch = 19)

points(fp_offline_positions, fp_offline_values, col = "green", pch = 19)
#text(fp_offline_positions, fp_offline_values + 0.005,
     #labels = paste0(round(, 2), "%"), cex = 0.7, pos = 3, col = "darkgreen")
legend("topright",
       legend = c("Borne inférieure (IC 95%)", 
                  "Borne supérieure (IC 95%)", 
                  "Taux brut hors-ligne", 
                  "Taux corrigé offline"),
       col = c("blue", "red", "green", "darkgreen"),
       lwd = c(2, 2, NA, NA),
       pch = c(NA, NA, 19, 19),
       bty = "n")

#Plot online

# Tracé des bornes
plot(p_seq, conf$lower, type = "l", col = "blue", lwd = 2,
     ylim = c(0, 0.1), xlab = "p (proportion pour n0 = p * N)",
     ylab = "Intervalle de confiance de la proportion",
     main = "IC à 95% pour Bin(n0, alpha = 0.05) pour online")
lines(p_seq, conf$upper, col = "red", lwd = 2)

fp_online_positions = c(1.00, 0.98, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.60)

fp_online_values = round(fp_online/100,2)


fp_online_values_corr = round(fp_online_corr/100,2)


points(fp_online_positions, fp_online_values_corr, col = "darkgreen", pch = 19)

points(fp_online_positions, fp_online_values, col = "green", pch = 19)
#text(fp_offline_positions, fp_offline_values + 0.005,
#labels = paste0(round(, 2), "%"), cex = 0.7, pos = 3, col = "darkgreen")
legend("topright",
       legend = c("Borne inférieure (IC 95%)", 
                  "Borne supérieure (IC 95%)", 
                  "Taux brut onlinr", 
                  "Taux corrigé onligne"),
       col = c("blue", "red", "green", "darkgreen"),
       lwd = c(2, 2, NA, NA),
       pch = c(NA, NA, 19, 19),
       bty = "n")

#Plot streaming

# Tracé des bornes
plot(p_seq, conf$lower, type = "l", col = "blue", lwd = 2,
     ylim = c(0, 0.1), xlab = "p (proportion pour n0 = p * N)",
     ylab = "Intervalle de confiance de la proportion",
     main = "IC à 95% pour Bin(n0, alpha = 0.05) pour streaming")
lines(p_seq, conf$upper, col = "red", lwd = 2)

fp_streaming_positions = c(1.00, 0.98, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.60)

fp_streaming_values = round(fp_streaming/100,2)

fp_streaming_values_corr = round(fp_streaming_corr/100,2)

points(fp_streaming_positions, fp_streaming_values_corr, col = "darkgreen", pch = 19)

points(fp_streaming_positions, fp_streaming_values, col = "green", pch = 19)

fp_vrais_positions = c(1.00, 0.98, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.60)

fp_vrais = round(fp_vrai/100,2)

fp_streaming_values_corr = round(fp_streaming_corr/100,2)

points(fp_vrais_positions, fp_vrais, col = "darkgreen", pch = 19)


points(fp_streaming_positions, fp_streaming_values_corr, col = "darkgreen", pch = 19)



fp_cov_positions = c(1.00, 0.98, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.60)

fp_cov_values = round(fp_cov/100,2)


points(fp_cov_positions, 1.1*fp_cov_values, col = "darkgreen", pch = 19)


points(fp_streaming_positions, fp_streaming_values_corr, col = "darkgreen", pch = 19)


#text(fp_offline_positions, fp_offline_values + 0.005,
#labels = paste0(round(, 2), "%"), cex = 0.7, pos = 3, col = "darkgreen")
legend("topright",
       legend = c("Borne inférieure (IC 95%)", 
                  "Borne supérieure (IC 95%)", 
                  "Taux brut streaming", 
                  "Taux corrigé streaming"),
       col = c("blue", "red", "green", "darkgreen"),
       lwd = c(2, 2, NA, NA),
       pch = c(NA, NA, 19, 19),
       bty = "n")


set.seed(123)
tune_lambda <- function(Z,mu_hat,Sigma, n_grid = 1000, alpha = 0.05,epsilon = 1e-5) {
  
  
  
  eig <- eigen(Sigma)
  eigvecs <- eig$vectors
  eigvals <- eig$values
  
  lambda_min_allowed <- -min(eigvals) + epsilon
  lambda_max <- 50
  
  grid_lambda <- seq(lambda_min_allowed, lambda_max, length.out = n_grid)
  results <- data.frame(lambda = grid_lambda, false_pos_rate = NA)
  
  
  for (i in seq_along(grid_lambda)) {
    lambda <- grid_lambda[i]
    eigvals_reg <- eigvals + lambda*rep(1,ncol(Z))
    #print(lambda)
    dists <- apply(Z, 1, function(x) {
mahalanobis_generalizedRcpp(x,mu_hat,eigvecs,eigvals_reg)    })
    
    threshold <- qchisq(1 - alpha, df = ncol(Z))
    #print(threshold)
    false_pos_rate <- mean(dists > threshold)
    
    results$false_pos_rate[i] <- false_pos_rate
  }
  
  # Choix de lambda le plus proche de 5 % de faux positifs
  best <- results[which.min(abs(results$false_pos_rate - alpha)), ]
  return(best)
}


best_lambda <- tune_lambda(Z,resultats$moyennem,resultats$Sigma[nrow(Z),,])
print(best_lambda)


fauxpositifs <- function(Z, data,alpha = 0.05,methode = "offline") {
  d <- ncol(Z)
  
  
  
  
  if (methode == "offline"){
  resultats <- OfflineOutlierDetection(Z)
  mu_hat = resultats$median
  Sigma = resultats$variance
  }
  # Étape 2 : Détection initiale sans régularisation
  cutoff <- qchisq(1 - alpha, df = d)
  
  outliers = resultats$outlier_labels
  
  tc <- table(data$labelsVrais[1:nrow(Z)], as.numeric(outliers))
  tc <- safe_access_tc(tc)
  
  if ((tc["0", "0"] + tc["0", "1"]) != 0) {
    fp_offline <- (tc["0", "1"] / (tc["0", "1"] + tc["0", "0"])) * 100
  }
  
  # Étape 3 : Recherche de lambda optimal
  best_lambda <- tune_lambda(Z, mu_hat , Sigma, alpha = alpha)$lambda
  Sigma_reg <- Sigma + best_lambda * diag(d)
  
  # Étape 4 : Détection régularisée
  distances_corr <- rep(0, nrow(Z))
  outliers_corr <- rep(0, nrow(Z))
  
  eig <- eigen(Sigma_reg)
  
  for (i in 1:nrow(Z)) {
    distances_corr[i] <- mahalanobis_generalizedRcpp(
      Z[i, ],
      mu_hat,
      eig$vectors,
      eig$values
    )
    if (distances_corr[i] > cutoff) {
      outliers_corr[i] <- 1
    }
  }
  
  tc <- table(data$labelsVrais[1:nrow(Z)], as.numeric(outliers_corr))
  tc <- safe_access_tc(tc)
    if ((tc["0", "0"] + tc["0", "1"]) != 0) {
    fp_offline_corr <- (tc["0", "1"] / (tc["0", "1"] + tc["0", "0"])) * 100
  }
  
  return(list(
    variance = resultats$Sigma2,
    median = resultats$mean,
    distances = distances,
    outlier_labels = outliers,
    distances_reg = distances_corr,
    outlier_labels_reg = outliers_corr,
    lambda_opt = best_lambda,
    fp_offline = fp_offline,
    fp_offline_corr = fp_offline_corr
  ))
}

fp_offline = rep(0,length(taux_contamination))

fp_offline_corr = rep(0,length(taux_contamination))
compt = 1
for (r in taux_contamination){
  
  data <- genererEchantillon(n,d,mu1,mu2,p1 = 1- r/100,r/100,Sigma1,Sigma2,"moyenne")
  
  Z = data$Z
  
  res = fauxpositifs(Z,data)
  fp_offline[compt] = res$fp_offline
  fp_offline_corr[compt] = res$fp_offline_corr
  #fp_offline_corr[compt]
  compt = compt + 1
}

fp_offline
