#install.packages("microbenchmark")
#install.packages("matlib")
#install.packages("MASS")
#install.packages("reshape2")
#install.packages("corrplot")
#install.packages("dlpyr")
# library(reshape2)
# library(RobRegression)
# library(Gmedian)
# library(ggplot2)
# library(far)
# library(gridExtra)
# library(microbenchmark)
# library("matlib")
# library("MASS")
# library("corrplot")
# library("dplyr")

#Facteur de correction de la distance de Mahalanobis

correctionDistanceMahalanobis <- function(distances,Z,methode = "offline")
{
  if(methode == "offline"){
 return(sqrt(qchisq(0.5,df = ncol(Z)))/(median(distances))) }
  else {
    facteur = rep(0,nrow(Z))
    for (i in (1:nrow(Z)))
    {
      facteur[i] = qchisq(0.5,df = ncol(Z))/(median(sqrt(distances[1:i])))^2
    }
    return (facteur)
  }
}

#Calcule le vecteur distances selon la méthode choisie

calcule_RMSE_FP_AUC_par_methode <- function(data, methode, cutoff = qchisq(0.95, df = ncol(data$Z)),SigmaVrai = Sigma1,muVrai = mu1) {
  distances <- NULL
  Z <- data$Z
  #faux_negatifs <- 0
  start_time <- Sys.time()  # Démarrage du chronomère
  
  if (methode == "Shrinkage") {
    # Méthode Shrinkage
    #med <- covComed(Z)$center
    #Sigma <- covCor(Z)
    
    med <- shrinkage_med(Z)$muShrink
    Sigma <- shrinkage_SCCM(Z,k = 1)$SCCM_shrinked
    
    
    #print(Sigma)
    rmseSigma <- norm(Sigma - SigmaVrai,"F")
    #rmseSigma <- 0
    print("OK Shrinkage")
    rmseMed <- sqrt(sum((muVrai - med)^2))
    distances <- calcule_vecteur_distances(Z, med, Sigma)
    outliers <- detectionOutliers(distances,cutoff)
    
    tc <- table(data$labelsVrais[1:(nrow(Z))], as.numeric(outliers)[1:(nrow(Z))])
    tc
    tc <- safe_access_tc(tc)
    if((tc["0","0"] + tc["0","1"]) != 0)
    {faux_positifs   <- round((tc["0", "1"]/(tc["0", "1"] + tc["0", "0"]))*100,2)}
    if((tc["1","0"] + tc["1","1"]) != 0)
    {faux_negatifs  <- round((tc["1", "0"]/(tc["1", "1"] + tc["1", "0"]))*100,2)}
    else faux_negatifs <- 0
    
    #auc <- round(auc(outliers,data$labelsVrais),2)*100
    if (length(unique(outliers)) == 2) {
      auc <- round(auc(outliers, data$labelsVrais), 2) * 100
    } else {
      auc <- 50  # Valeur par défaut pour un cas non exploitable
    }
  }
    else if (methode == "FASTMCD")
    {
      res <- covMcd(Z)
      med <- res$center
      Sigma <- res$cov
      
      #print(Sigma)
      rmseSigma <- norm(Sigma - SigmaVrai,"F")
      #rmseSigma <- 0
      print("OK FastMCD")
      rmseMed <- sqrt(sum((muVrai - med)^2))
      distances <- calcule_vecteur_distances(Z, med, Sigma)
      outliers <- detectionOutliers(distances,cutoff)
      
      tc <- table(data$labelsVrais[1:(nrow(Z))], as.numeric(outliers)[1:(nrow(Z))])
      tc
      tc <- safe_access_tc(tc)
      if((tc["0","0"] + tc["0","1"]) != 0)
      {faux_positifs   <- round((tc["0", "1"]/(tc["0", "1"] + tc["0", "0"]))*100,2)}
      else faux_positifs <- 0
      if((tc["1","0"] + tc["1","1"]) != 0)
      {faux_negatifs   <- round((tc["1", "0"]/(tc["1", "0"] + tc["1", "1"]))*100,2)}
      else faux_negatifs <- 0
      if (length(unique(outliers)) == 2) {
        auc <- round(auc(outliers, data$labelsVrais), 2) * 100
      } else {
        auc <- 50  # Valeur par défaut pour un cas non exploitable
      }
      
   
  } else if (methode == "Comédiane") {
    # Méthode Comédiane
    res <- covComed(Z,n.iter = 0)
    med <- res$center
    Sigma <- res$cov
    #print(Sigma)
    
    rmseSigma <- norm(Sigma - SigmaVrai,"F")
    print("OK Comed")
    rmseMed <- sqrt(sum((muVrai - med)^2))
    distances <- calcule_vecteur_distances(Z, med, Sigma)
    #rmseSigma <- norm(Sigma - Sigma1,"F")
    #rmseMed <- sqrt(sum((mu1 - med)^2))
    distances <- calcule_vecteur_distances(Z, med, Sigma)
    outliers <- detectionOutliers(distances,cutoff)
    
    tc <- table(data$labelsVrais[1:(nrow(Z))], as.numeric(outliers)[1:(nrow(Z))])
    tc
    tc <- safe_access_tc(tc)
    if((tc["0","0"] + tc["0","1"]) != 0)
    {faux_positifs   <- round((tc["0", "1"]/(tc["0", "1"] + tc["0", "0"]))*100,2)}
    if((tc["1","0"] + tc["1","1"]) != 0)
    {faux_negatifs  <- round((tc["1", "0"]/(tc["1", "1"] + tc["1", "0"]))*100,2)}
    else faux_negatifs <- 0
    
    auc <- round(auc(outliers,data$labelsVrais),2)*100
    
  } else if (methode == "OGK") {
    # Méthode OGK
    OGK_result <- covOGK(Z,sigmamu = s_mad)
    med <- OGK_result$center
    Sigma <- OGK_result$cov
    rmseSigma <- norm(Sigma - SigmaVrai,"F")
    print("OK OGK")
    rmseMed <- sqrt(sum((muVrai - med)^2))
    distances <- calcule_vecteur_distances(Z, med, Sigma)
    outliers <- detectionOutliers(distances,cutoff)
    
    tc <- table(data$labelsVrais[1:(nrow(Z))], as.numeric(outliers)[1:(nrow(Z))])
    tc
    tc <- safe_access_tc(tc)
    if((tc["0","0"] + tc["0","1"]) != 0)
    {faux_positifs   <- round((tc["0", "1"]/(tc["0", "1"] + tc["0", "0"]))*100,2)}
    if((tc["1","0"] + tc["1","1"]) != 0)
    {faux_negatifs  <- round((tc["1", "0"]/(tc["1", "1"] + tc["1", "0"]))*100,2)}
    else faux_negatifs <- 0
    
    auc <- round(auc(outliers,data$labelsVrais),2)*100
    
    
  } else if (methode == "Cov Empirique") {
    med <- colMeans(Z)
    Sigma <- cov(Z)
    distances <- calcule_vecteur_distances(Z, med, Sigma)
    rmseSigma <- norm(Sigma - SigmaVrai,"F")
    print("OK Cov Emp")
    rmseMed <- sqrt(sum((muVrai - med)^2))
    outliers <- detectionOutliers(distances,cutoff)
    
    tc <- table(data$labelsVrais[1:(nrow(Z))], as.numeric(outliers)[1:(nrow(Z))])
    tc
    tc <- safe_access_tc(tc)
    if((tc["0","0"] + tc["0","1"]) != 0)
    {faux_positifs   <- round((tc["0", "1"]/(tc["0", "1"] + tc["0", "0"]))*100,2)}
    if((tc["1","0"] + tc["1","1"]) != 0)
    {faux_negatifs  <- round((tc["1", "0"]/(tc["1", "1"] + tc["1", "0"]))*100,2)}
    else faux_negatifs <- 0
    
    
    auc <- round(auc(outliers,data$labelsVrais),2)*100
    
  } else if (methode == "offline") {
    # Méthode Offline
    #med <- RobVar(Z)$median
    #resultats <- estimation(Z,methodeEstimation = "offline")
    resultats = OfflineOutlierDetection(Z)
    med <- resultats$median
    Sigma <- resultats$variance
    #rmseSigma <- ifelse(norm(resultats$SigmaOfflinePoids,"F") < 0.001,norm(resultats$SigmaOffline - SigmaVrai,"F"),min(norm(resultats$SigmaOfflinePoids - SigmaVrai,"F"),norm(resultats$SigmaOffline- SigmaVrai,"F")))
    rmseSigma = norm(Sigma - SigmaVrai,"F")
    print("OK offline")
    rmseMed <- sqrt(sum((muVrai - med)^2))
    distances <- resultats$distances
    #distances = calcule_vecteur_distances(Z,mu1,SigmaVrai)
    #distancescorr = distances*correctionDistanceMahalanobis(distances,Z)
    #distances <- calcule_vecteur_distances(Z,med,Sigma)
    #outliers <- detectionOutliers(distances,cutoff = qchisq(0.95,df = ncol(Z))) 
    outliers <- resultats$outlier_labels
    tc <- table(data$labelsVrais[1:(nrow(Z))], as.numeric(outliers)[1:(nrow(Z))])
    tc
    tc <- safe_access_tc(tc)
    if((tc["0","0"] + tc["0","1"]) != 0)
    {faux_positifs   <- round((tc["0", "1"]/(tc["0", "1"] + tc["0", "0"]))*100,2)}
    if((tc["1","0"] + tc["1","1"]) != 0)
    {faux_negatifs  <- round((tc["1", "0"]/(tc["1", "1"] + tc["1", "0"]))*100,2)}
    else faux_negatifs <- 0
    
    #auc <- round(auc(outliers,data$labelsVrais),2)*100
    if (length(unique(outliers)) == 2) {
      auc <- round(auc(outliers, data$labelsVrais), 2) * 100
    } else {
      auc <- 50  # Valeur par défaut pour un cas non exploitable
    }
    #Vérification que la distance de Mahalanobis suit la loi du Khi-deux
#     hist(distances, probability = TRUE, breaks = 30,    ylim = c(0,.1),
# 
#          main = "Histogramme des distances de Mahalanobis",
#          xlab = "Distance", col = "lightblue", border = "white")
# curve(dchisq(x, df = ncol(Z)), col = "red", lwd = 2, add = TRUE)
#     legend("topright", legend = c("Empirique", "Théorique Chi²"),
#            fill = c("lightblue", NA), border = c("white", NA), lty = c(NA, 1), col = c(NA, "red"))
  } else if (methode == "online") {
    # Méthode Online
    #resultats <- estimation(Z,methodeEstimation = "online")
    
    resultats = StreamingOutlierDetection(Z,batch = 1,mc_sample_size = 1)
    med <- resultats$moyennem
    Sigma <- resultats$Sigma[nrow(Z),,]
    #rmseSigma <- ifelse(norm(resultats$SigmaOnlinePoids,"F") < 0.001,norm(resultats$SigmaOnline - SigmaVrai,"F"),min(norm(resultats$SigmaOnlinePoids - SigmaVrai,"F"),norm(resultats$SigmaOnline - SigmaVrai,"F")))
    
    rmseSigma <- norm(Sigma - SigmaVrai,"F")
    print("OK online")
    rmseMed <- sqrt(sum((muVrai - med)^2))
    #miter <- results$miter
    #U <- results$U
    #lambda <- results$lambdaIter
    distances <- resultats$distances
    outliers = resultats$outlier_labels
    #outliers <- detectionOutliers(distances,cutoff)
    
    tc <- table(data$labelsVrais[1:(nrow(Z))], as.numeric(outliers)[1:(nrow(Z))])
    tc
    tc <- safe_access_tc(tc)
    if((tc["0","0"] + tc["0","1"]) != 0)
    {faux_positifs   <- round((tc["0", "1"]/(tc["0", "1"] + tc["0", "0"]))*100,2)}
    if((tc["1","0"] + tc["1","1"]) != 0)
    {faux_negatifs  <- round((tc["1", "0"]/(tc["1", "1"] + tc["1", "0"]))*100,2)}
    else faux_negatifs <- 0
    
    if (length(unique(outliers)) == 2) {
      auc <- round(auc(outliers, data$labelsVrais), 2) * 100
    } else {
      auc <- 50  
    }
    # #Vérification qu'on a une Khi-Deux
    # hist(distances, probability = TRUE, breaks = 30,
    #      main = "Histogramme des distances de Mahalanobis",
    #      xlab = "Distance", col = "lightblue", border = "white")
    # 
    # curve(dchisq(x, df = ncol(Z)), col = "red", lwd = 2, add = TRUE)
    # legend("topright", legend = c("Empirique", "Théorique Chi²"),
    #        fill = c("lightblue", NA), border = c("white", NA), lty = c(NA, 1), col = c(NA, "red"))
    } 
  else if (methode == "streaming") {
    # Méthode streaming
    #resultats <- estimation(Z,methodeEstimation = "streaming")
    resultats = StreamingOutlierDetection(Z,batch = ncol(Z))
    
    med <- resultats$moyennem
    Sigma <- resultats$Sigma[nrow(Z),,]
    rmseSigma <- norm(Sigma - SigmaVrai,"F")
    
    #rmseSigma <- ifelse(norm(resultats$SigmaStreamingPoids,"F") < 0.01,norm(resultats$SigmaStreaming - SigmaVrai,"F"),min(norm(resultats$SigmaStreamingPoids - SigmaVrai,"F"),norm(resultats$SigmaStreaming - SigmaVrai,"F")))
       
    print("OK Streaming")
    rmseMed <- sqrt(sum((muVrai - med)^2))
    #miter <- results$miter
    #U <- results$U
    #lambda <- results$lambdaIter
    distances <- resultats$distances
    #distances = calcule_vecteur_distances(Z,resultats$SigmaStreamingPoids,cutoff)
    outliers = resultats$outlier_labels
    #outliers <- detectionOutliers(distances,cutoff)
    
    tc <- table(data$labelsVrais[1:(nrow(Z))], as.numeric(outliers)[1:(nrow(Z))])
    tc
    tc <- safe_access_tc(tc)
    if((tc["0","0"] + tc["0","1"]) != 0)
    {faux_positifs   <- round((tc["0", "1"]/(tc["0", "1"] + tc["0", "0"]))*100,2)}
    if((tc["1","0"] + tc["1","1"]) != 0)
    {faux_negatifs   <- round((tc["1", "0"]/(tc["1", "1"] + tc["1", "0"]))*100,2)}
    else faux_negatifs <- 0
    if (length(unique(outliers)) == 2) {
      auc <- round(auc(outliers, data$labelsVrais),2) * 100
    } else {
      auc <- 50  
    }
    # #Vérification distance Khi-deux
    # hist(distances, probability = TRUE, breaks = 30,
    #      main = "Histogramme des distances de Mahalanobis",
    #      xlab = "Distance", col = "lightblue", border = "white")
    # 
    # curve(dchisq(x, df = ncol(Z)), col = "red", lwd = 2, add = TRUE)
    # legend("topright", legend = c("Empirique", "Théorique Chi²"),
    #        fill = c("lightblue", NA), border = c("white", NA), lty = c(NA, 1), col = c(NA, "red"))
  } 
  # Fin du chrono
  end_time <- Sys.time()
  temps_calcul <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  return(list(distances = distances,med = med,Sigma = Sigma,rmseMed = rmseMed,rmseSigma = rmseSigma,faux_positifs = faux_positifs,faux_negatifs = faux_negatifs,temps_calcul = temps_calcul,auc = auc))
}


#Calcul des distances de Mahalanobis à partir de la médiane géométrique et de la covariance estimées

calcule_vecteur_distances <- function(Z,m,Sigma)

{

  distances <- rep(0,n)

  for  (i in 1:nrow(Z))
  {
    distances[i] = as.numeric(Z[i,] - m) %*% solve(Sigma) %*% (as.numeric(t(Z[i,] - m)))

  }

  return (distances)
}

# Fonction pour calculer le cutoff (deuxième seuil)

calcule_cutoff <- function(distances, d) {

  # Constantes et quantiles
  quantile_95 <- qchisq(0.95, df = d)
  quantile_50 <- qchisq(0.5, df = d)
  median_d <- median(distances)

  # Calcul du cutoff
  cv <- (quantile_95 * median_d) / quantile_50

  return(cv)
}


#Détection des outliers à partir d'un vecteur de distances de Mahalanobis et d'un cutoff


detectionOutliers <- function(distances,cutoff)
  
{
  
  outliers_labels <- rep(0,length(distances))

  
 
  
  for (i in 1:length(distances))
  {
    

    if (distances[i] > cutoff) {outliers_labels[i] <- 1}
    else {outliers_labels[i] <- 0}
  }
  
  return (outliers_labels)
  
}


