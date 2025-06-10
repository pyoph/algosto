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


#########################################
AUC
########################################


aucTout = resMoyenne$aucRec


auc_moy <- apply(aucTout, c(1, 2), mean)

methodes_sel <- c(1, 2, 9, 10)

# Sous-tableau avec uniquement les méthodes sélectionnées
auc_sel <- auc_moy[, methodes_sel]

# Convertir en data.frame avec noms clairs
df_auc <- data.frame(
  Taux = 2:9,
  M1 = auc_sel[2:9, 1],
  M2 = auc_sel[2:9, 2],
  M9 = auc_sel[2:9, 3],
  M10 = auc_sel[2:9, 4]
)

# Transformer en long format pour ggplot2
df_long <- melt(df_auc, id.vars = "Taux", variable.name = "Méthode", value.name = "AUC")

# Tracer avec ggplot2
ggplot(df_long, aes(x = Taux, y = AUC, color = Méthode)) +
  geom_line(size = 1) +
  geom_point() +
  labs(title = "AUC selon les taux de contamination",
       x = "Taux de contamination", y = "AUC") +
  theme_minimal() +
  theme(legend.position = "bottom")