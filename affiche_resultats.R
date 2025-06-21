
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(patchwork)  # pour agencer les 4 graphiques

###################################################
#Affiche contamination scenarios
##################################################



afficherContaminationScenarios = function(k,l,rho,contamination,rate){
  
  sigmaSq0 <- (1:d); sigmaSq0 <- sigmaSq0 / mean(sigmaSq0)
  SigmaContamin = diag(sqrt(sigmaSq0)) %*% toeplitz(rho^(0:(d-1))) %*% diag(sqrt(sigmaSq0))

  data <- genererEchantillon(n,d,mu1,mu2 = k*rep(1/sqrt(d), d),p1 = 1- rate/100,rate/100,Sigma1,Sigma2 = l*Sigma2,contamin = contamination,cluster = FALSE)
  
 Z = data$Z
 # Création d'un dataframe pour ggplot
 df <- data.frame(X1 = data$Z[,1], X2 = data$Z[,2], 
                  Label = factor(data$labelsVrais, levels = c(0, 1)))
 #Création du plot
 # Create the plot (English version)
 p <- ggplot(df, aes(x = X1, y = X2, color = Label)) +
   geom_point() +
   scale_color_manual(
     values = c("0" = "blue", "1" = "red"),
     labels = c("Inlier", "Outlier"),
     name = "Status"  # Removed the empty second name
   ) +
   labs(
     title = "Contamination Scenario: Mean and Variance",
     subtitle = paste("Rate:", rate, "%, k =", k, ", l =", l, ", rho =", rho),
     x = "X1",  # Added axis labels
     y = "X2"
   ) +
   theme_minimal() +
   theme(
     legend.position = "bottom",  # Optional: moves legend to bottom
     plot.title = element_text(face = "bold")  # Optional: makes title bold
   )
  return(p)
}

p0 = afficherContaminationScenarios(0,1,0.3,"moyenne_variance",20)

p1 = afficherContaminationScenarios(1,0.1,0.6,"moyenne_variance",20)

p2 = afficherContaminationScenarios(30,0.01,0.995,"moyenne_variance",20)


# Créer une disposition 2x2 avec les 3 plots et un espace vide

library(cowplot)

# Créer une grille 2x2 avec p0, p1, p2 et un espace vide
plot_grid(
  p0, p1,  # Première ligne (p0 et p1 côte à côte)
  p2, NULL, # Deuxième ligne (p2 + case vide)
  ncol = 2,
  align = "hv",  # Alignement horizontal et vertical
  rel_widths = c(1, 1),  # Largeurs égales
  rel_heights = c(1, 1)   # Hauteurs égales
)



######################################################################
#######RMSE Sigma Iter
######################################################################

rmseSigma = res100runNearesScenario$rmseSigmaRec

rmseSigma_moy = res1run$rmseSigmaRec[,,,1]

dim(rmseSigma)
#rmseSigma_moy = rmseSigma[,,,1]

rmseSigma_moy <- apply(rmseSigma, c(1, 2, 3), mean)

#rmseSigma_moy = res1run$rmseSigmaRec[,,,1]


#rmseSigma_moy = res1run$rmseSigmaRec[,,,1]
# On prépare les données en long format pour ggplot2
# Extraire les méthodes et taux souhaités
methodes <- c(1, 9, 10)
taux_indices <- c(3, 4, 5, 8)

method_labels <- c(
  "1" = "Sample covariance (online)",
  "9" = "Online",
  "10" = "Streaming"
)

rate_labels <- c(
  "3" = "5%",
  "4" = "10%",
  "5" = "20%",
  "8" = "30%"
)

# Mappage des identifiants méthode à leur position dans le tableau (axe 3)
method_position_map <- setNames(1:length(methodes), methodes)

data_list <- list()
n <- dim(rmseSigma_moy)[1]

for (j in taux_indices) {
  for (k in methodes) {
    k_pos <- method_position_map[as.character(k)]
    df <- data.frame(
      index = 1:n,
      RMSE = rmseSigma_moy[, j, k_pos],
      Method = method_labels[as.character(k)],
      Rate = rate_labels[as.character(j)]
    )
    data_list[[length(data_list) + 1]] <- df
  }
}

########################################
#Frobenius norm error iterations
#######################################

df_long <- do.call(rbind, data_list)
df_long$Rate <- factor(df_long$Rate, levels = c( "5%", "10%", "20%","30%"))

gg <- ggplot(df_long, aes(x = index, y = RMSE, color = Method)) +
  geom_line(size = 0.8) +
  facet_wrap(~ Rate, ncol = 2) +
  labs(
    title = "Frobenius Norm Error Across Different Contamination Rates",
    subtitle = "Comparison of Online Covariance Estimation Methods (nearest scenario)",
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

print(gg)


# #########################################
#             AUC
# ########################################


aucTout = res100runNearesScenario$aucRec


#aucTout = res1run$aucRec

auc_moy <- apply(aucTout, c(1, 2), mean)

auc_moy = res10run$aucRec[,,1]
#auc_moy = res$aucRec[,,1]
#auc_moy = res$aucRec[,,1]


methodes <- c(1, 9, 10)
method_labels <- c("1" = "Cov Online", "9" = "Online", "10" = "Streaming")
method_pos_map <- setNames(1:3, methodes)  # map méthode → colonne

taux_indices <- 1:9
taux_valeurs <- c(2, 5, 10, 15, 20, 25, 30, 35, 40)

data_list <- list()
for (i in seq_along(taux_indices)) {
  j <- taux_indices[i]
  for (k in methodes) {
    k_pos <- method_pos_map[as.character(k)]
    df <- data.frame(
      ContaminationRate = taux_valeurs[i],
      AUC = auc_moy[j, k_pos],
      Method = method_labels[as.character(k)]
    )
    data_list[[length(data_list) + 1]] <- df
  }
}
df_long <- do.call(rbind, data_list)
df_long$Method <- factor(df_long$Method, levels = method_labels)

# Graphique
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
    breaks = taux_valeurs
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")


dim(res10run$faux_positifsRec)

#########################################
#Faux négatifs moyenne
########################################


res1run$faux_negatifsRec[1e4,9,3,1]


faux_neg_moy = apply(res1run$faux_negatifsRec, c(1, 2, 3), mean)

faux_neg_10000 <- faux_neg_moy[10000, , ]


methodes <- c(1, 9, 10)
method_labels <- c("1" = "Cov Online", "9" = "Online", "10" = "Streaming")
method_pos_map <- setNames(1:3, methodes)

taux_valeurs <- c(2, 5, 10, 15, 20, 25, 30, 35, 40)
data_list <- list()
for (i in seq_along(taux_valeurs)) {
  total_outliers = taux_valeurs[i]/100*n
    
    for (k in methodes) {
    k_pos <- method_pos_map[as.character(k)]
    data_list[[length(data_list) + 1]] <- data.frame(
      ContaminationRate = taux_valeurs[i],
      FalseNegatives = faux_neg_10000[i, k_pos]/(total_outliers)*100,
      Method = method_labels[as.character(k)]
    )
  }
}
df_long <- do.call(rbind, data_list)
df_long$Method <- factor(df_long$Method, levels = method_labels)

ggplot(df_long, aes(x = ContaminationRate, y = FalseNegatives, color = Method)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "Faux Négatifs vs. Taux de Contamination",
    x = "Taux de Contamination (%)",
    y = "Faux Négatifs",
    color = "Méthode"
  ) +
  scale_x_continuous(breaks = taux_valeurs) +
  theme_minimal() +
  theme(legend.position = "bottom")


#########################################
#Faux positifs moyenne
########################################


fprec = res100runNearesScenario$faux_positifsRec
fprec = res1run$faux_positifsRec
fprec_moy <- apply(fprec, c(1, 2, 3), mean)

fprec100run = fprec_moy[1e4,,]

#dim(fprec1run)

methodes <- c(1, 9, 10)
method_labels <- c("1" = "Cov Online", "9" = "Online", "10" = "Streaming")
method_pos_map <- setNames(1:3, methodes)

taux_valeurs <- c(0,2, 5, 10, 15, 20, 25, 30,  40)
data_list <- list()
for (i in seq_along(taux_valeurs)) {
  total_inliers = (1 - taux_valeurs[i]/100)*n
  
  for (k in methodes) {
    k_pos <- method_pos_map[as.character(k)]
    data_list[[length(data_list) + 1]] <- data.frame(
      ContaminationRate = taux_valeurs[i],
      FalsePositives = fprec100run[i, k_pos]/(total_inliers)*100,
      Method = method_labels[as.character(k)]
    )
  }
}
df_long <- do.call(rbind, data_list)
df_long$Method <- factor(df_long$Method, levels = method_labels)

pnearDist = ggplot(df_long, aes(x = ContaminationRate, y = FalsePositives, color = Method)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "False Positives (near distribution scenario) (k,l,rho) = (2,1,0.6)",
    x = "Contamination rates (%)",
    y = "False positives",
    color = "Methods"
  ) +
  scale_x_continuous(breaks = taux_valeurs) +
  theme_minimal() +
  theme(legend.position = "bottom")

####Sauvegarde du graphique pour les faux positifs

saveRDS(pnearDist, file = "pnearDistFP.rds")

# ##############################################
# Temps calculs
# ################# #############################

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

outliers_labelsTout = res100runNearesScenario$outliersLabelsRec
#outliers_labelsTout = res1run$outliersLabelsRec[,,,1]

# Fonction pour vote majoritaire
majority_vote <- function(x) {
  ux <- unique(x)
  counts <- tabulate(match(x, ux))
  max_count <- max(counts)
  major_values <- ux[counts == max_count]
  
  if (length(major_values) == 1) {
    return(major_values)
  } else {
    return(1)
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
    geom_line(aes(y = Detected_rate, color = "True and false positive rate"), size = 0.5) +
    geom_line(aes(y = True_positive_rate, color = "True positives rate"), size = 0.5) +
    geom_line(aes(y = True_outliers, color = "True outliers rate"), size = 0.5) +
    scale_color_manual(values = c(
      "True outliers rate" = "red",
      "True and false positive rate" = "orange",
      "True positives rate" = "purple"
    )) +
    #scale_x_log10() +
    labs(
      title = paste(pourcentage, "% of outliers"),
      x = "Data index",
      y = "Cumulative rate (%)",
      color = "Legend"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  # print(taux_outliers_vrais[1:10])
  return(list(p = p,taux_outliers_vrais = taux_outliers_vrais,taux_outliers_detectes = taux_outliers_detectes,taux_outliers_detectes_vrais = taux_outliers_detectes_vrais))
}

#outliers_majority = outliers_labelsTout[,,1]
# res1run$faux_negatifsRec[9900:1e4]
# table(res1run$labelsVraisRec[,5],res1run$outliersLabelsRec[,5,9,1])
dim(res10run$outliersLabelsRec)

########Graphs

pCumOutDetRateNearScOnl5 = cumulativeOutlierDetection(res1run$labelsVraisRec[,3],res1run$outliersLabelsRec[,3,2,1],5,"")

pCumOutDetRateNearScOnl10 = cumulativeOutlierDetection(res1run$labelsVraisRec[,4],res1run$outliersLabelsRec[,4,2,1],10,"")

pCumOutDetRateNearScOnl20 = cumulativeOutlierDetection(res1run$labelsVraisRec[,6]$labelsVraisRec[,6],res1run$outliersLabelsRec[,6,2,1],20,"")

pCumOutDetRateNearScOnl30 = cumulativeOutlierDetection(res1run$labelsVraisRec[,8],outliers_majority[,8,2],30,"")

pnearDist

install.packages("cowplot")  # Une seule fois
library(cowplot)

# Crée des cases vides avec draw_plot(NULL)
plot_grid(
  pnearDist,
  pCumOutDetRateNearScOnl5[[1]],
  pCumOutDetRateNearScOnl10[[1]],
  pCumOutDetRateNearScOnl20[[1]],
  pCumOutDetRateNearScOnl30[[1]],
  ncol = 2
  #labels = "AUTO"  # Optionnel : ajoute des lettres A, B, C...
)


pCumStrmDetRateNearSc5 = cumulativeOutlierDetection(res1run$labelsVraisRec[,3],res1run$outliersLabelsRec[,3,3,1],5,"")

pCumStrmDetRateNearSc10 = cumulativeOutlierDetection(res1run$labelsVraisRec[,4],res1run$outliersLabelsRec[,4,3,1],10,"")

pCumStrmDetRateNearSc20 = cumulativeOutlierDetection(res1run$labelsVraisRec[,6],res1run$outliersLabelsRec[,6,3,1],20,"")

pCumOutDetRateNearSc30 = cumulativeOutlierDetection(res1run$labelsVraisRec[,8],res1run$outliersLabelsRec[,6,3,1],30,"")

# Crée des cases vides avec draw_plot(NULL)
plot_grid(
  pnearDist,
  pCumStrmDetRateNearSc5 [[1]],
  pCumStrmDetRateNearSc10 [[1]],
  pCumStrmDetRateNearSc20 [[1]],
  pCumStrmDetRateNearSc30[[1]],
  ncol = 2
  #labels = "AUTO"  # Optionnel : ajoute des lettres A, B, C...
)

