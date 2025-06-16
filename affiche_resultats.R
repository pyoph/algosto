
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(patchwork)  # pour agencer les 4 graphiques



rmseSigma = res1run$rmseSigmaRec

dim(resMoyenne$rmseSigmaRec)

dim(rmseSigma)
rmseSigma_moy = rmseSigma[,,,1]

rmseSigma_moy <- apply(rmseSigma, c(1, 2, 3), mean)

#rmseSigma_moy = res1run$rmseSigmaRec[,,,1]


#rmseSigma_moy = res1run$rmseSigmaRec[,,,1]
# On prépare les données en long format pour ggplot2
# Extraire les méthodes et taux souhaités

methodes = c(1,9,10)

taux_indices = c(1,3,5,9)

# Noms lisibles pour les méthodes
method_labels <- c(
  "1" = "Sample covariance (online)",
  "9" = "Online",
  "10" = "Streaming"
)

# Noms lisibles pour les taux
rate_labels <- c(
  "1" = "0%",
  "3" = "5%",
  "5" = "15%",
  "9" = "40%"
)

df_long = list()
# Création des données en format long
data_list <- list()
for (j in taux_indices) {
  for (k in methodes) {
    df <- data.frame(
      index = 1:n,
      RMSE = rmseSigma_moy[, j, k],
      Method = method_labels[as.character(k)],
      Rate = rate_labels[as.character(j)]
    )
    data_list[[length(data_list) + 1]] <- df
  }
}
df_long <- do.call(rbind, data_list)

# Forcer l’ordre d’affichage des taux
df_long$Rate <- factor(df_long$Rate, levels = c("0%", "5%", "15%", "40%"))

# Graphique
gg <- ggplot(df_long, aes(x = index, y = RMSE, color = Method)) +
  geom_line(size = 0.8) +
  facet_wrap(~ Rate, ncol = 2) +
  labs(
    title = "Frobenius Norm Error Across Different Contamination Rates",
    subtitle = "Comparison of Covariance Estimation Methods (shifted Gaussian with reduction of variance)",
    x = "Observation Index",
    y = "Frobenius Norm Error"
  ) +
  theme_minimal(base_size = 12) +
  scale_y_log10(
    breaks = 10^seq(0, 5, by = 1),
    labels = c("1", "10", "100", "1000", "10000", "100000")
  ) +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "bottom") +
  annotation_logticks(sides = "l")

# Affichage
print(gg)


methodes <- c(1, 9, 10)

taux_contamination <- c(0, 2, 5, 10, 15, 20, 25, 30, 40)


method_labels <- c(
  "1" = "Sample covariance (online)",
  "9" = "Online",
  "10" = "Streaming"
)

# Taux de contamination correspondants à j = 1:9
#taux_contamination <- seq(0, 40, length.out = 9)

# Extraction du dernier RMSE (ligne 10000) pour chaque taux et méthode
df_summary <- data.frame()
for (j in 1:9) {
  for (k in methodes) {
    rmse_value <- rmseSigma_moy[1e4, j, k]
    df_summary <- rbind(df_summary, data.frame(
      Taux = taux_contamination[j],
      Methode = method_labels[as.character(k)],
      RMSE = rmse_value
    ))
  }
}

# Graphique
gg <- ggplot(df_summary, aes(x = Taux, y = RMSE, color = Methode)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  labs(
    title = "Final Frobenius norm error vs. Contamination rate (shifted Gaussian with variance reduction)",
    x = "Contamination rate (%)",
    y = "Final Frobenius norm error"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(gg)



# #########################################
#             AUC
# ########################################


aucTout = resMoyenne$aucRec

aucTout = res$aucRec

#aucTout = res1run$aucRec

auc_moy <- apply(aucTout, c(1, 2), mean)

auc_moy = res1run$aucRec[,,1]
#auc_moy = res$aucRec[,,1]
#auc_moy = res$aucRec[,,1]

methodes_sel <- c(1, 9, 10)

# Taux de contamination et méthodes à utiliser
taux_indices <- 2:9
methodes <- c(1, 9, 10)
taux_valeurs <- c(2, 5, 10, 20, 25, 40)

# Construction du data frame long
data_list <- list()
for (i in seq_along(taux_indices)) {
  j <- taux_indices[i]
  for (k in methodes) {
    df <- data.frame(
      ContaminationRate = taux_valeurs[i],
      AUC = auc_moy[j, k],
      Method = factor(k, levels = c(1, 2, 9, 10),
                      labels = c("Cov Online", "Cov Online Trimmed", "Online", "Streaming"))
    )
    data_list[[length(data_list) + 1]] <- df
  }
}
df_long <- do.call(rbind, data_list)

ggplot(df_long, aes(x = ContaminationRate, y = AUC, color = Method)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "AUC vs. Contamination Rate",
    x = "Contamination Rate (%)",
    y = "AUC",
    color = "Method"
  ) +
  scale_x_continuous(
    #limits = c(2, max(df_long$ContaminationRate)),
    breaks = c(2, 5,10, 15, 20, 25, 30,40)
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")


# ##############################################
# Temps calculs
# ##############################################

temps_calculTout = res$temps_calcul
  

# Méthodes et indices souhaités
methodes <- c("sampleCovOnline", "samplecovTrimmed", "sampleCovOffline", "comedianeOffline",
              "comedianeOfflineShrinkage", "OGK", "FASTMCD", "offline", "online", "streaming")

# Indices à garder : 1, 2, 6 à 10
indices_gardes <- c(1, 6, 7, 8, 9, 10)

# Supposons que taux_index est défini
taux_index <- 3  # par exemple

# Extraction des données pour ce taux et méthodes sélectionnées
temps_sel <- temps_calculTout[taux_index, indices_gardes, ]  # dims : méthodes sélectionnées x runs

# Transformation en data frame long
df_temps <- melt(temps_sel, varnames = c("MethodeIndex", "Run"), value.name = "Temps")

# Remplacement par les noms des méthodes sélectionnées
df_temps$Methode <- factor(df_temps$MethodeIndex, 
                           levels = 1:length(indices_gardes), 
                           labels = methodes[indices_gardes])

# Plot boxplot
ggplot(df_temps, aes(x = Methode, y = Temps)) +
  geom_boxplot(fill = "lightblue", outlier.color = "red", outlier.shape = 1) +
  labs(
    title = "Boxplot of computation times",
    x = "Method",
    y = "Time (seconds)"
  ) +
  scale_y_log10() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

######Outliers labels######################

outliers_labelsTout = resMoyenne$outliersLabelsRec

outliers_labelsTout = res$outliersLabelsRec
outliers_labelsTout = res1run$outliersLabelsRec[,,,1]

# Fonction pour vote majoritaire
majority_vote <- function(x, tie = c("first", "min", "all")) {
  tie <- match.arg(tie)
  
  ux <- unique(x)
  counts <- tabulate(match(x, ux))
  max_count <- max(counts)
  major_values <- ux[counts == max_count]
  
  if (length(major_values) == 1 || tie == "all") {
    return(major_values)
  } else if (tie == "first") {
    return(major_values[1])
  } else if (tie == "min") {
    return(min(major_values))
  }
}


# Appliquer vote majoritaire sur la 4e dimension
outliers_majority <- apply(outliers_labelsTout, c(1, 2, 3), majority_vote)

#outliers_majority = res$outliersLabelsRec[,,,1]

cumulativeOutlierDetection <- function(labelsVrais,outlier_labels , pourcentage,titre) {
  total_points <- length(labelsVrais)
  total_outliers_theoriques <- pourcentage / 100 * total_points
  
  nb_outliers_detectes <- 0
  nb_outliers_detectes_vrais <- 0
  nb_outliers_vrais <- 0
  
  taux_outliers_detectes <- numeric(total_points)
  taux_outliers_vrais <- numeric(total_points)
  taux_outliers_detectes_vrais = numeric(total_points)
  
  for (i in 1:total_points) {
    if(labelsVrais[i] == 1) nb_outliers_vrais =  nb_outliers_vrais + 1
    
    if (labelsVrais[i] == 1 & outlier_labels[i] == 1) {
      nb_outliers_detectes <- nb_outliers_detectes + 1
      #nb_outliers_vrais <- nb_outliers_vrais + 1
      nb_outliers_detectes_vrais = nb_outliers_detectes_vrais + 1
      #taux_outliers_vrais[i] <- nb_outliers_detectes_vrais
    } else if (labelsVrais[i] == 0 && outlier_labels[i] == 1) {
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
    
    # print(paste("nb_outliers_detectes_vrais ",nb_outliers_detectes_vrais ))
    # print(paste("nb_outliers_vrais ",nb_outliers_vrais ))
    # print(paste("taux_outliers_vrais[i] ",taux_outliers_vrais[i] ))
    # print(paste("taux_outliers_vrais[i] ",taux_outliers_vrais[i] ))
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
  # print(taux_outliers_vrais[1:10])
  return(list(p = p,taux_outliers_vrais = taux_outliers_vrais,taux_outliers_detectes = taux_outliers_detectes,taux_outliers_detectes_vrais = taux_outliers_detectes_vrais))
}

outliers_majority = outliers_labelsTout[,,1]
res1run$faux_negatifsRec[9900:1e4]
table(res1run$labelsVraisRec[,5],res1run$outliersLabelsRec[,5,9,1])

cumulativeOutlierDetection(res1run$labelsVraisRec[,5],res1run$outliersLabelsRec[,5,9,1],15,"Shifted Gaussian contamination scenario")