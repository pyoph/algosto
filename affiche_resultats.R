
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(patchwork)  # pour agencer les 4 graphiques



rmseSigma = resMoyenne$rmseSigmaRec

dim(resMoyenne$rmseSigmaRec)

dim(rmseSigma)

rmseSigma = resMoyenne$rmseSigmaRec[,,,]


rmseSigma_moy <- apply(rmseSigma, c(1, 2, 3), mean)







# On prépare les données en long format pour ggplot2
# Extraire les méthodes et taux souhaités
taux_indices <- c(1, 3, 5, 9)
methodes <- c(9, 10)

# Create long format data.frame
data_list <- list()
for (j in taux_indices) {
  for (k in methodes) {
    df <- data.frame(
      index = 1:n,
      RMSE = rmseSigma_moy[, j, k],
      Method = paste0("Method ", k),
      Rate = paste0("Rate ", j)
    )
    data_list[[length(data_list) + 1]] <- df
  }
}
df_long <- do.call(rbind, data_list)

# Create plot with ggplot2
gg <- ggplot(df_long, aes(x = index, y = RMSE, color = Method)) +
  geom_line(size = 0.8) +
  scale_color_manual(values = rainbow(length(unique(df_long$Methode)))) +
  facet_wrap(~ Rate, ncol = 2) +
  labs(title = "Frobenius Norm Error Across Different Contamination Rates",
       subtitle = "Comparison of Robust Covariance Estimation Methods",
       x = "Observation Index", 
       y = "Frobenius Norm Error (log scale)") +
  theme_minimal(base_size = 12) +
  scale_y_log10() +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold")) +
  scale_color_brewer(palette = "Set1")

# Display the plot
print(gg)



# #########################################
# AUC
# ########################################


aucTout = resMoyenne$aucRec


auc_moy <- apply(aucTout, c(1, 2), mean)

methodes_sel <- c(1, 2, 9, 10)

# Taux de contamination et méthodes à utiliser
taux_indices <- 2:9
methodes <- c(1, 2, 9, 10)
taux_valeurs <- c(2, 5, 10, 15, 20, 25, 30, 40)

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
    breaks = c(2, 5, 10, 15, 20, 25, 30, 40)
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")


# ##############################################
# Temps calculs
# ##############################################

temps_calculTout = resMoyenne$temps_calcul
  

# Méthodes et indices souhaités
methodes <- c("sampleCovOnline", "samplecovTrimmed", "sampleCovOffline", "comedianeOffline",
              "comedianeOfflineShrinkage", "OGK", "FASTMCD", "offline", "online", "streaming")

# Indices à garder : 1, 2, 6 à 10
indices_gardes <- c(1, 2, 6, 7, 8, 9, 10)

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







