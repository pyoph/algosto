#library(reshape2)
#library(RobRegression)
#library(Gmedian)
#library(ggplot2)
#library(far)
#library(gridExtra)
#library(microbenchmark)
#library("matlib")
#library("MASS")
#library("corrplot")
#library("dplyr")
#library(ROCR)
#install.packages("ROCR")

#Calcule les erreurs d'estimation de la médiane géométrique critère RMSE

calculErreursM <- function(miter,mvrai)
{
  nb_iterations <- dim(miter)[2]
  erreursM <- rep(0,nb_iterations)
  
  for (i in 1:(nb_iterations-1))
  {
    erreursM[i] <- sqrt(sum((mvrai - miter[i,])^2))
    
  }
  return (erreursM)
  
}

#Calcul des erreurs en norme de Frobenius

calculErreursNormeFrobenius <- function (Estim,Vrai)
{
  nbiterations <- nrow(Estim)
  
  erreurs <- rep(0,nbiterations)
  
  for (i in (1:nbiterations))
  {
    erreurs[i] <-  norm(Vrai - Estim,type = "F")
  }
  return (erreurs)
}



# Fonction pour afficher les erreurs d'estimation de m online et offline pour différents taux de contamination. 

#Les erreurs offline seront par la suite ajoutées
affiche_erreursM <- function(erreurs_online, erreurs_offline = NULL, contamination = 2) {
  
  nberreurs_online <- length(erreurs_online)
  
  if (!is.null(erreurs_offline)) {
    nberreurs_offline <- length(erreurs_offline)
    
    if (nberreurs_online != nberreurs_offline) {
      stop("Les vecteurs 'erreurs_online' et 'erreurs_offline' doivent avoir la même longueur.")
    }
    
    # Créer le data frame combiné
    data <- data.frame(
      x = rep(1:(nberreurs_online - 1), 2),
      erreurs = c(erreurs_online[1:(nberreurs_online - 1)], erreurs_offline[1:(nberreurs_offline - 1)]),
      type = c(rep("Online", nberreurs_online - 1), rep("Offline", nberreurs_offline - 1))
    )
    
    # Définir les couleurs et labels
    colors <- c("Online" = "blue", "Offline" = "green")
    labels <- c(
      paste("Online (contamination = ", contamination , "%)", sep = ""),
      paste("Offline (contamination = ", contamination , "%)", sep = "")
    )
  } else {
    # Créer le data frame uniquement pour les erreurs online
    data <- data.frame(
      x = 1:(nberreurs_online - 1),
      erreurs = erreurs_online[1:(nberreurs_online - 1)],
      type = rep("Online", nberreurs_online - 1)
    )
    
    # Définir les couleurs et labels pour online uniquement
    colors <- c("Online" = "blue")
    labels <- c(paste("Online (contamination = ", contamination , "%)", sep = ""))
  }
  
  # Créer le graphique
  p <- ggplot(data, aes(x = x, y = erreurs, color = type, linetype = type)) +
    geom_line(size = 1) +
    scale_x_log10() +   # Échelle logarithmique pour l'axe x
    scale_y_log10() +   # Échelle logarithmique pour l'axe y
    scale_color_manual(values = colors, labels = labels) +
    labs(x = "Taille de l'échantillon", 
         y = "Erreur RMSE",
         title = "Comparaison des erreurs d'estimation de m",
         color = "Méthode", 
         linetype = "Méthode") +
    theme_minimal()
  
  print(p)
}

#Boucle de construction du tableau de résultats

construction_tableau_resultats <- function(nbrunsParam = nbruns,contamin = "moyenne")
{
  methodes = c("Comédiane","Shrinkage","OGK","Cov Empirique","offline","online","streaming")
  taux_contamination <- c(0,2, 5, 10, 15, 20, 25, 30, 40)
  rmseMed <- matrix(0,length(taux_contamination),length(methodes))
  rmseSigma <- matrix(0,length(taux_contamination),length(methodes))
  faux_positifs <- matrix(0,length(taux_contamination),length(methodes))
  faux_negatifs <- matrix(0,length(taux_contamination),length(methodes))
  temps_calcul  <- matrix(0,length(taux_contamination),length(methodes))
  auc <- matrix(0,length(taux_contamination),length(methodes))
  for (i in seq_along(taux_contamination)){
  #Génération des données 
    
    delta <- taux_contamination[i]
    p1 <- 1 - delta/100
    p2 <- 1 - p1
    data <- genererEchantillon(n,d,mu1,mu2,p1,p2 ,Sigma1,Sigma2,contamin = contamin)
    
    for (k in (1:nbrunsParam)){
    for (j in seq_along(methodes))
  {
      
    #m = "Shrinkage"
    methode <- methodes[j]
    resultats <- calcule_RMSE_FP_AUC_par_methode(data,methode = methode)
    rmseMed[i,j] <- resultats$rmseMed + rmseMed[i,j]
    rmseSigma[i,j] <- resultats$rmseSigma +  rmseSigma[i,j]
    faux_positifs[i,j] <- resultats$faux_positifs + faux_positifs[i,j]
    
    faux_negatifs[i,j] <- resultats$faux_negatifs + faux_negatifs[i,j]
    temps_calcul[i,j] <- resultats$temps_calcul + temps_calcul[i,j]
    auc[i,j] <- resultats$auc + auc[i,j]
    }
     
    }
  }
  
  return(list(resultats = resultats, rmseMed = rmseMed/nbruns,rmseSigma = rmseSigma/nbruns,faux_positifs = faux_positifs/nbruns,faux_negatifs = faux_negatifs/nbruns,auc =auc/nbruns,temps_calcul= temps_calcul/nbruns))
  }


#Construction du dataset

RMSEAUCFPdataset<- function(nbrunsParam = nbruns,contamin = "moyenne")

{
  
  #Initialisation des taux de contamination
  taux_contamination <- c(0,2, 5, 10, 15, 20, 25, 30, 40)
  
  
  #Initialisation faux positifs  
  
  faux_positifs_offline <- rep(0,(length(taux_contamination)))
  faux_positifs_online <- rep(0,(length(taux_contamination)))
  faux_positifs_streaming <- rep(0,(length(taux_contamination)))
  faux_positifs_covEmp <- rep(0,(length(taux_contamination)))
  faux_positifs_ogk <- rep(0,(length(taux_contamination)))
  faux_positifs_comed <- rep(0,(length(taux_contamination)))
  faux_positifs_shrink <- rep(0,(length(taux_contamination)))
  
  #Intialisation faux négatifs
  
  faux_negatifs_offline <- rep(0,(length(taux_contamination)))
  faux_negatifs_online <- rep(0,(length(taux_contamination)))
  faux_negatifs_streaming <- rep(0,(length(taux_contamination)))
  faux_negatifs_covEmp <- rep(0,(length(taux_contamination)))
  faux_negatifs_ogk <- rep(0,(length(taux_contamination)))
  faux_negatifs_comed <- rep(0,(length(taux_contamination)))
  faux_negatifs_shrink <- rep(0,(length(taux_contamination)))
  
  
  # Initialisation des vecteurs pour les temps de calcul
  temps_calcul_offline <- rep(0, length(taux_contamination))
  temps_calcul_online <- rep(0, length(taux_contamination))
  temps_calcul_streaming <- rep(0, length(taux_contamination))
  temps_calcul_covEmp <- rep(0, length(taux_contamination))
  temps_calcul_ogk <- rep(0, length(taux_contamination))
  temps_calcul_comed <- rep(0, length(taux_contamination))
  temps_calcul_shrink <- rep(0, length(taux_contamination))
  
  #   
  #   #initialisation RMSE médiane
  #   
  #   
  rmse_med_offline <- rep(0,(length(taux_contamination)))
  rmse_med_online <- rep(0,(length(taux_contamination)))
  rmse_med_streaming <- rep(0,(length(taux_contamination)))
  rmse_med_covEmp <- rep(0,(length(taux_contamination)))
  rmse_med_ogk <- rep(0,(length(taux_contamination)))
  rmse_med_comed <- rep(0,(length(taux_contamination)))
  rmse_med_shrink <- rep(0,(length(taux_contamination)))
  
  
  
  
  
  #   #initialisation RMSE Sigma
  
  
  rmse_Sigma_offline <- rep(0,(length(taux_contamination)))
  rmse_Sigma_online <- rep(0,(length(taux_contamination)))
  rmse_Sigma_streaming <- rep(0,(length(taux_contamination)))
  rmse_Sigma_covEmp <- rep(0,(length(taux_contamination)))
  rmse_Sigma_ogk <- rep(0,(length(taux_contamination)))
  rmse_Sigma_comed <- rep(0,(length(taux_contamination)))
  rmse_Sigma_shrink <- rep(0,(length(taux_contamination)))
  
  
  
  
  #   #initialisation AUC
  
  auc_offline <- rep(0,(length(taux_contamination)))
  auc_online <- rep(0,(length(taux_contamination)))
  auc_streaming <- rep(0,(length(taux_contamination)))
  auc_covEmp <- rep(0,(length(taux_contamination)))
  auc_ogk <- rep(0,(length(taux_contamination)))
  auc_comed <- rep(0,(length(taux_contamination)))
  auc_shrink <- rep(0,(length(taux_contamination)))
  
  # Appeler la fonction pour obtenir les résultats
  resultats <- construction_tableau_resultats(nbrunsParam = nbruns,contamin = contamin)
  
  # Extraire les matrices des résultats
  rmseMed <- resultats$rmseMed
  rmseSigma <- resultats$rmseSigma
  faux_positifs <- resultats$faux_positifs
  faux_negatifs <- resultats$faux_negatifs
  temps_calcul <- resultats$temps_calcul
  auc <- resultats$auc
  
  # Remplir les vecteurs pour chaque méthode
  for (i in seq_along(taux_contamination)) {
    faux_positifs_comed[i] <- faux_positifs[i, 1]  # Comédiane
    faux_positifs_shrink[i] <- faux_positifs[i, 2] # Shrinkage
    faux_positifs_ogk[i] <- faux_positifs[i, 3]    # OGK
    faux_positifs_covEmp[i] <- faux_positifs[i, 4] # Cov Empirique
    faux_positifs_offline[i] <- faux_positifs[i, 5] # Offline
    faux_positifs_online[i] <- faux_positifs[i, 6]  # Online
    faux_positifs_streaming[i] <- faux_positifs[i, 7] # Streaming
    
    faux_negatifs_comed[i] <- faux_negatifs[i, 1]  # Comédiane
    faux_negatifs_shrink[i] <- faux_negatifs[i, 2] # Shrinkage
    faux_negatifs_ogk[i] <- faux_negatifs[i, 3]    # OGK
    faux_negatifs_covEmp[i] <- faux_negatifs[i, 4] # Cov Empirique
    faux_negatifs_offline[i] <- faux_negatifs[i, 5] # Offline
    faux_negatifs_online[i] <- faux_negatifs[i, 6]  # Online
    faux_negatifs_streaming[i] <- faux_negatifs[i, 7] # Streaming
    
    
    rmse_med_comed[i] <- rmseMed[i, 1]  # Comédiane
    rmse_med_shrink[i] <- rmseMed[i, 2] # Shrinkage
    rmse_med_ogk[i] <- rmseMed[i, 3]    # OGK
    rmse_med_covEmp[i] <- rmseMed[i, 4] # Cov Empirique
    rmse_med_offline[i] <- rmseMed[i, 5] # Offline
    rmse_med_online[i] <- rmseMed[i, 6]  # Online
    rmse_med_streaming[i] <- rmseMed[i, 7] # Streaming
    
    rmse_Sigma_comed[i] <- rmseSigma[i, 1]  # Comédiane
    rmse_Sigma_shrink[i] <- rmseSigma[i, 2] # Shrinkage
    rmse_Sigma_ogk[i] <- rmseSigma[i, 3]    # OGK
    rmse_Sigma_covEmp[i] <- rmseSigma[i, 4] # Cov Empirique
    rmse_Sigma_offline[i] <- rmseSigma[i, 5] # Offline
    rmse_Sigma_online[i] <- rmseSigma[i, 6]  # Online
    rmse_Sigma_streaming[i] <- rmseSigma[i, 7] # Streaming
    
    auc_comed[i] <- auc[i, 1]  # Comédiane
    auc_shrink[i] <- auc[i, 2] # Shrinkage
    auc_ogk[i] <- auc[i, 3]    # OGK
    auc_covEmp[i] <- auc[i, 4] # Cov Empirique
    auc_offline[i] <- auc[i, 5] # Offline
    auc_online[i] <- auc[i, 6]  # Online
    auc_streaming[i] <- auc[i, 7] # Streaming

    
    temps_calcul_comed[i] <- temps_calcul[i, 1]  # Comédiane
    temps_calcul_shrink[i] <- temps_calcul[i, 2] # Shrinkage
    temps_calcul_ogk[i] <- temps_calcul[i, 3]    # OGK
    temps_calcul_covEmp[i] <- temps_calcul[i, 4] # Cov Empirique
    temps_calcul_offline[i] <- temps_calcul[i, 5] # Offline
    temps_calcul_online[i] <- temps_calcul[i, 6]  # Online
    temps_calcul_streaming[i] <- temps_calcul[i, 7] # Streaming
      }
  
  # Créer le dataframe avec les champs RMSE_Sigma, RMSE_Med, AUC et FP pour chaque méthode
  results_metrics <- data.frame(
    RMSE_Sigma_Cov = rmse_Sigma_covEmp,
    RMSE_Med_Cov = rmse_med_covEmp,
    AUC_Cov = auc_covEmp,
    FP_Cov = faux_positifs_covEmp,
    FN_Cov = faux_negatifs_covEmp,
    TC_Cov = temps_calcul_covEmp,
    RMSE_Sigma_OGK = rmse_Sigma_ogk,
    RMSE_Med_OGK = rmse_med_ogk,
    AUC_OGK = auc_ogk,
    FP_OGK = faux_positifs_ogk,
    FN_OGK = faux_negatifs_ogk,
    TC_OGK = temps_calcul_ogk,
    RMSE_Sigma_Comed = rmse_Sigma_comed,
    RMSE_Med_Comed = rmse_med_comed,
    AUC_Comed = auc_comed,
    FP_Comed = faux_positifs_comed,
    FN_Comed = faux_negatifs_comed,
    TC_Comed = temps_calcul_comed,
    RMSE_Sigma_Shrink = rmse_Sigma_shrink,
    RMSE_Med_Shrink = rmse_med_shrink,
    AUC_Shrink = auc_shrink,
    FP_Shrink = faux_positifs_shrink,
    FN_Shrink = faux_negatifs_shrink,
    TC_Shrink = temps_calcul_shrink,
    RMSE_Sigma_Online = rmse_Sigma_online,
    RMSE_Med_Online = rmse_med_online,
    AUC_Online = auc_online,
    FP_Online = faux_positifs_online,
    FN_Online = faux_negatifs_online,
    TC_Online = temps_calcul_online,
    RMSE_Sigma_Offline = rmse_Sigma_offline,
    RMSE_Med_Offline = rmse_med_offline,
    AUC_Offline = auc_offline,
    FP_Offline = faux_positifs_offline,
    FN_Offline = faux_negatifs_offline,
    TC_offline = temps_calcul_offline,
    RMSE_Sigma_Streaming = rmse_Sigma_streaming,
    RMSE_Med_Streaming = rmse_med_streaming,
    AUC_Streaming = auc_streaming,
    FP_Streaming = faux_positifs_streaming,
    FN_Streaming = faux_negatifs_streaming,
    TC_Streaming = temps_calcul_streaming
  )
  
  row.names(results_metrics) <- taux_contamination
  
  # Afficher le dataframe
  print(results_metrics)
  # 
  # results_without_RMSE_Med <- results_metrics %>%
  #   select(-contains("RMSE_Med")) %>%  # Supprimer toutes les colonnes RMSE_Med
  #   select(-contains("Cov"))
  # 
  # # 2. Dataframe sans RMSE_Sigma
  # results_without_RMSE_Sigma <- results_metrics %>%
  #   select(-contains("RMSE_Sigma"))
  # 
  # # Afficher les dataframes
  # print("Dataframe sans RMSE_Med :")
  # print(results_without_RMSE_Med)
  # 
  # print("Dataframe sans RMSE_Sigma :")
  # print(results_without_RMSE_Sigma)
  # results_without_RMSE_Med <- round(results_without_RMSE_Med,2)
  # 
  # save(results_without_RMSE_Med,file = "results_without_RMSE_Med.RData")
  # results_metrics <- round(results_metrics,2)
  # latex_table_results_metrics <- xtable(results_metrics)
  # 
  # latex_table_results_metrics_sansRMSESigma <- xtable( results_without_RMSE_Sigma)
  # latex_table_results_metrics_sansRMSEMed <- xtable( results_without_RMSE_Med)
  # 
  # 
  save(results_metrics,file = "results_metrics.RData")
  
  return(results_metrics)
}


# Fonction pour afficher les erreurs d'estimation de Sigma online et offline pour différents taux de contamination. 

#Les erreurs offline seront par la suite ajoutées
affiche_erreursSigma <- function(erreurs_online, erreurs_offline = NULL, erreurs_Str = NULL,contamination = 0) {
  
  nberreurs_online <- length(erreurs_online)
  
  if (!is.null(erreurs_offline)) {
    nberreurs_offline <- length(erreurs_offline)
    
    if (nberreurs_online != nberreurs_offline) {
      stop("Les vecteurs 'erreurs_online' et 'erreurs_offline' doivent avoir la même longueur.")
    }
    
    # Créer le data frame combiné
    data <- data.frame(
      x = rep(1:(nberreurs_online - 1), 2),
      erreurs = c(erreurs_online[1:(nberreurs_online - 1)], erreurs_offline[1:(nberreurs_offline - 1)]),
      type = c(rep("Online", nberreurs_online - 1), rep("Offline", nberreurs_offline - 1))
    )
    
    # Définir les couleurs et labels
    colors <- c("Online" = "blue", "Offline" = "green")
    labels <- c(
      paste("Online (contamination = ", contamination , "%)", sep = ""),
      paste("Offline (contamination = ", contamination , "%)", sep = "")
    )
  } else {
    # Créer le data frame uniquement pour les erreurs online
    data <- data.frame(
      x = 1:(nberreurs_online - 1),
      erreurs = erreurs_online[1:(nberreurs_online - 1)],
      type = rep("Online", nberreurs_online - 1)
    )
    
    # Définir les couleurs et labels pour online uniquement
    colors <- c("Online" = "blue")
    labels <- c(paste("Online (contamination = ", contamination , "%)", sep = ""))
  }
  
  # Créer le graphique
  p <- ggplot(data, aes(x = x, y = erreurs, color = type)) +
    geom_line(size = 1) +
    scale_x_log10() +   # Échelle logarithmique pour l'axe x
    scale_y_log10() +   # Échelle logarithmique pour l'axe y
    scale_color_manual(values = colors, labels = labels) +
    labs(x = "Taille de l'échantillon", 
         y = "Erreur norme de Frobenius",
         title = "Comparaison des erreurs d'estimation de Sigma",
         color = "Méthode") +   # Une seule rubrique pour "Méthode"
    theme_minimal()
  
  print(p)
}


#Affichage comparaison Sigma

plot_comparaison_sigma <- function(Sigma1, SigmaEstimOnline, SigmaEstimOffline, delta) {
  par(mfrow = c(1, 2))  # Afficher deux graphiques côte à côte
  
  # Comparaison Online
  plot(Sigma1, SigmaEstimOnline, col = 2, 
       main = paste("Comparaison Online (Contamination:", delta, "%)"), 
       xlab = "Sigma réelle", ylab = "Estimations Online")
  abline(0, 1)
  points(Sigma1, SigmaEstimOnline, col=2)
  
  # Comparaison Offline
  plot(Sigma1, SigmaEstimOffline, col = 4, 
       main = paste("Comparaison Offline (Contamination:", delta, "%)"), 
       xlab = "Sigma réelle", ylab = "Estimation Offline")
  abline(0, 1)
  
  points(Sigma1, SigmaEstimOffline, col=2)
  par(mfrow = c(1, 1))  # Réinitialisation du paramètre graphique
}


#Affichage boxplot des erreurs

creer_boxplot_erreurs <- function(erreursSigmaBoxplot, taux_contamination, methode,erreursMonteCarlo ) {
  # Vérifier que la matrice d'erreurs a bien des lignes correspondant aux taux de contamination
  if (length(taux_contamination) != nrow(erreursSigmaBoxplot)) {
    stop("Le nombre de lignes de erreursSigmaBoxplot doit correspondre à la longueur de taux_contamination")
  }
  
  
  # Conversion de la matrice en data frame long
  df_erreurs <- data.frame(erreursSigmaBoxplot)
  
  df_monte_carlo <- data.frame(
    Taux_Contamination = as.factor(taux_contamination),  
    MonteCarlo = erreursMonteCarlo,
    Type = "Erreur Monte Carlo" 
  )
  
  df_erreurs$taux_contamination <- factor(taux_contamination)  # Facteur pour l'axe X
  df_long <- melt(df_erreurs, id.vars = "taux_contamination")
  
  # Création du boxplot avec ggplot2
  p <- ggplot(df_long, aes(x = taux_contamination, y = value)) +
    geom_boxplot(fill = "lightblue", color = "blue", outlier.color = "red") +
    #geom_jitter(width = 0.2, alpha = 0.2, color = "black") +  # Ajout de dispersion pour la visibilité
    #geom_crossbar(data = df_monte_carlo, 
    #  aes(y = MonteCarlo, ymin = MonteCarlo, ymax = MonteCarlo, color = "Monte Carlo"), 
    # width = 0.7, fatten = 2, size = 1.2, linetype = "dashed") +  # Une seule ligne par taux
    geom_line(data = df_monte_carlo, aes(x = Taux_Contamination, y = MonteCarlo, group = 1), 
              color = "red", size = 1.2) +   
    scale_y_log10() + 
    scale_fill_manual(values = c("Boxplot Erreurs" = "lightblue")) +  # 
    scale_color_manual(values = c("Boxplot Erreurs" = "blue", "Erreur Monte Carlo" = "red")) +  
    labs(title = paste("Erreur de l'estimation de Sigma (", methode, ")"),
         x = "Taux de contamination (%)",
         y = "Erreur norme de Frobenius",
         color = "Méthode",  # Légende des couleurs
         fill = "Méthode") +  # Légende des boxplots
    theme_minimal() +
    annotate("text", x = 5, y = max(erreursMonteCarlo, na.rm = TRUE) * 1.5, 
             label = "Trait rouge =  E(||S_{n_0} - V||_F^2) estimation par Monte Carlo",
             color = "black", size = 5, hjust = 0)  + # 
    theme(legend.position = "top")  
  print(p)
}

#Fonction pour tracer les courbes ROC

Performance <- function(score, class, thresh, add=FALSE, col=1){
  perf <- performance(prediction(score, class), "tpr", "fpr")
  tab <- cbind(unlist(perf@alpha.values), unlist(perf@x.values), unlist(perf@y.values))
  numThresh <- min(which(unlist(perf@alpha.values) < thresh))
  perfThresh <- tab[numThresh, -1]; 
  names(perfThresh) <- c("tpr", "fpr")
  plot(perf, add=add, col=col); abline(0, 1, col=8)
  points(perfThresh[1], perfThresh[2], col=col, cex=2)
  return(list(perf=perf, tab=tab, perfThresh=perfThresh))
}


#Estimation de E(||1/n_0 sum Xi Xi^T - V||_F^2) par Monte Carlo

EstimVarMC <- function(nbiter,delta,Sigma)
  
{
  erreurs_MC <- rep(0,nbiter)
  
  for (i in (1 : nbiter))
  {
    p1 <- 1 - delta / 100
    
    p2 <- 1 - p1
    #mu1
    
    resultsSimul <- genererEchantillon(n, d, mu1, mu2, p1, p2, Sigma1 = Sigma1, Sigma2 = Sigma2,contamin = contamin)
    
    Z <- resultsSimul$Z
    
    n1 <- floor(p1*n)
    
    #Recherche des inliers
    inliers <- Z[resultsSimul$labelsVrais == 0, , drop = FALSE]
    
   M <- matrix(0,d,d)
   
   for (k in (1:nrow(inliers)))
   {
     M <- M + (inliers[k,])%*%t(inliers[k,])
   }
    
   
   erreurs_MC[i] <- norm(M/n1 - Sigma,"F")^2
   
   
   
  }
  return (mean(erreurs_MC))
  
}
