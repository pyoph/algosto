rmseSigma = resMoyenne$rmseSigmaRec

dim(resMoyenne$rmseSigmaRec)

dim(rmseSigma)

rmseSigma = resMoyenne$rmseSigmaRec[,,,]


rmseSigma_moy <- apply(rmseSigma, c(1, 2, 3), mean)






library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(patchwork)  # pour agencer les 4 graphiques

# Supposons que rmseSigma_moy est un tableau [10000 x 9 x 10]

# On prépare les données en long format pour ggplot2
# Extraire les méthodes et taux souhaités
taux_indices <- c(1, 4, 6, 9)
methodes <- c(1, 2, 9, 10)

# Création d'un data.frame long
data_list <- list()
for (j in taux_indices) {
  for (k in methodes) {
    df <- data.frame(
      index = 1:10000,
      RMSE = rmseSigma_moy[, j, k],
      Methode = paste0("M", k),
      Taux = paste0("Taux ", j)
    )
    data_list[[length(data_list) + 1]] <- df
  }
}
df_long <- do.call(rbind, data_list)

# Plot avec ggplot2
gg <- ggplot(df_long, aes(x = index, y = RMSE, color = Methode)) +
  geom_line(size = 0.8) +
  facet_wrap(~ Taux, ncol = 2) +
  labs(title = "RMSE Sigma pour différents taux de contamination",
       x = "Index des données", y = "RMSE Sigma") +
  theme_minimal() +
  scale_y_log10() +
  theme(legend.position = "bottom")

# Afficher le graphique
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
    title = paste("Boxplot of computation times", taux_index),
    x = "Method",
    y = "Time (seconds)"
  ) +
  scale_y_log10() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))