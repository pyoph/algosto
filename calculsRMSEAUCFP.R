#install.packages("easystats")
#install.packages("bigutilsr")
#install.packages("xtable")
#install.packages("RGMM")
#install.packages("rrcov")
#install.packages("mvnfast")
#library(xtable)
#library(reshape2)
#library(RobRegression)
#library("robustbase")
# library(Gmedian)
# library(ggplot2)
# library(far)
# library(gridExtra)
# library(microbenchmark)
# #library("matlib")
# library("MASS")
# library("corrplot")
# library("dplyr")
# library("easystats")
# library("bigutilsr")
# #library("rJava")
# #library("REPPlab")
# library("RGMM")

#Rajout de 0 dans les tables de contingence si les champs sont vides
safe_access_tc <- function(tc, default = 0) {
  # Identifier toutes les combinaisons possibles des noms de lignes et colonnes
  all_rows <- unique(c(rownames(tc), "0", "1"))  # Ajoutez vos besoins ici
  all_cols <- unique(c(colnames(tc), "0", "1"))
  
  # Initialiser une matrice complète
  full_tc <- matrix(default, nrow = length(all_rows), ncol = length(all_cols),
                    dimnames = list(all_rows, all_cols))
  
  # Remplir avec les valeurs existantes
  for (row in rownames(tc)) {
    for (col in colnames(tc)) {
      full_tc[row, col] <- tc[row, col]
    }
  }
  return(full_tc)
}


#########Calcul rmse, faux positifs et auc pour un dataset et toutes les méthodes



################Boucle construction tableau#############
# 
# calculeRMSEAUCFP <- function(nbruns = 20,cutoff = qchisq(p = 0.95,df = ncol(data)),contamin ="moyenne")
# {
# 
#   
#   taux_contamination <- c(0,2, 5, 10, 15, 20, 25, 30, 40)
#   
#   #initialisation faux positifs 
# 
#   faux_positifs_offline <- matrix(0,nbruns,(length(taux_contamination)))
#   faux_positifs_online <- matrix(0,nbruns,(length(taux_contamination)))
#   faux_positifs_streaming <- matrix(0,nbruns,(length(taux_contamination)))
#   faux_positifs_covEmp <- matrix(0,nbruns,(length(taux_contamination)))
#   faux_positifs_ogk <- matrix(0,nbruns,(length(taux_contamination)))
#   faux_positifs_comed <- matrix(0,nbruns,(length(taux_contamination)))
#   faux_positifs_shrink <- matrix(0,nbruns,(length(taux_contamination)))
#   
#   
#   
#   #initialisation RMSE médiane
#   
#   
#   rmse_med_offline <- matrix(0,nbruns,(length(taux_contamination)))
#   rmse_med_online <- matrix(0,nbruns,(length(taux_contamination)))
#   rmse_med_streaming <- matrix(0,nbruns,(length(taux_contamination)))
#   rmse_med_covEmp <- matrix(0,nbruns,(length(taux_contamination)))
#   rmse_med_ogk <- matrix(0,nbruns,(length(taux_contamination)))
#   rmse_med_comed <- matrix(0,nbruns,(length(taux_contamination)))
#   rmse_med_shrink <- matrix(0,nbruns,(length(taux_contamination)))
#   
#   
#   
#   
#   
#   #initialisation RMSE Sigma
#       
#   
#   rmse_Sigma_offline <- matrix(0,nbruns,(length(taux_contamination)))
#   rmse_Sigma_online <- matrix(0,nbruns,(length(taux_contamination)))
#   rmse_Sigma_streaming <- matrix(0,nbruns,(length(taux_contamination)))
#   rmse_Sigma_covEmp <- matrix(0,nbruns,(length(taux_contamination)))
#   rmse_Sigma_ogk <- matrix(0,nbruns,(length(taux_contamination)))
#   rmse_Sigma_comed <- matrix(0,nbruns,(length(taux_contamination)))
#   rmse_Sigma_shrink <- matrix(0,nbruns,(length(taux_contamination)))
#   
#   
#   
#   
#   #initialisation AUC
#   
#   auc_offline <- matrix(0,nbruns,(length(taux_contamination)))
#   auc_online <- matrix(0,nbruns,(length(taux_contamination)))
#   auc_streaming <- matrix(0,nbruns,(length(taux_contamination)))
#   auc_covEmp <- matrix(0,nbruns,(length(taux_contamination)))
#   auc_ogk <- matrix(0,nbruns,(length(taux_contamination)))
#   auc_comed <- matrix(0,nbruns,(length(taux_contamination)))
#   auc_shrink <- matrix(0,nbruns,(length(taux_contamination)))
#   
#   
#   for (i in seq_along(taux_contamination))
#   {
#     delta <- taux_contamination[i]
#     #delta <- 10
#     #contamin = "moyenne"
#     print(contamin)
#     p1 <- 1 - delta / 100
#     
#     p2 <- 1 - p1
#     #mu1
#     resultsSimul <- genererEchantillon(n, d, mu1, mu2, p1, p2, Sigma1 = Sigma1, Sigma2 = Sigma2,contamin = contamin)
#     Z <- resultsSimul$Z
#     
#     
#     for (j in (1:nbruns)){
#     
#       
#       #Méthode Offline
#    
#       resOffline <- estimation(Z,methodeEstimation = "offline")
#       SigmaOffline<- resOffline$Sigma
#       medOffline <- resOffline$med
#       rmse_Sigma_offline[j,i] <- norm(SigmaOffline - Sigma1,"F")
#       rmse_med_offline[j,i] <-  sqrt(sum((mu1 - medOffline)^2))
#       distances <- calcule_vecteur_distances(Z,medOffline,SigmaOffline)
#       #cutoff = calcule_cutoff(distances,d)
#       #cutoff =  qchisq(p = 0.95, df = ncol(Z))
#       #outliers_listMaha <- check_outliers(Z, method = "mahalanobis")
#       outliers_Offline <- detectionOutliers(distances,cutoff)
#       
#       
#       tc <- table(resultsSimul$labelsVrais[1:(nrow(Z))], as.numeric(outliers_Offline)[1:(nrow(Z))])
#       tc
#       tc <- safe_access_tc(tc)
#       if((tc["0","0"] + tc["0","1"]) != 0)
#       {faux_positifs_offline[i]   <- round((tc["0", "1"]/(tc["0", "1"] + tc["0", "0"]))*100,2)}
#       auc_offline <- round(auc(outliers_OGK,resultsSimul$labelsVrais),2)*100
#    
#       
#       #Méthode Online
#       
#       resOnline <- estimation(Z,methodeEstimation = "online")
#       SigmaOnline<- resOnline$Sigma
#       medOnline <- resOnline$med
#       rmse_Sigma_online[j,i] <- norm(SigmaOnline[nrow(Z) -1,,] - Sigma1,"F")
#       rmse_med_online[j,i] <-  sqrt(sum((mu1 - medOnline)^2))
#       distances <- calcule_vecteur_distances(Z,medOnline,SigmaOffline)
#       #cutoff = calcule_cutoff(distances,d)
#       #cutoff =  qchisq(p = 0.95, df = ncol(Z))
#       #outliers_listMaha <- check_outliers(Z, method = "mahalanobis")
#       outliers_Online <- detectionOutliers(distances,cutoff)
#       
#       
#       tc <- table(resultsSimul$labelsVrais[1:(nrow(Z))], as.numeric(outliers_Online)[1:(nrow(Z))])
#       tc
#       tc <- safe_access_tc(tc)
#       if((tc["0","0"] + tc["0","1"]) != 0)
#       {faux_positifs_online[i]   <- round((tc["0", "1"]/(tc["0", "1"] + tc["0", "0"]))*100,2)}
#       auc_online <- round(auc(outliers_Online,resultsSimul$labelsVrais),2)*100
#       
#       #Méthode Streaming
#       
#       resStreaming <- estimation(Z,methodeEstimation = "streaming")
#       SigmaStreaming<- resStreaming$Sigma
#       medStreaming <- resStreaming$med
#       rmse_Sigma_streaming[j,i] <- norm(SigmaStreaming[nrow(Z) -1,,] - Sigma1,"F")
#       rmse_med_streaming[j,i] <-  sqrt(sum((mu1 - medStreaming)^2))
#       distances <- calcule_vecteur_distances(Z,medStreaming,SigmaStreaming)
#       #cutoff = calcule_cutoff(distances,d)
#       #cutoff =  qchisq(p = 0.95, df = ncol(Z))
#       #outliers_listMaha <- check_outliers(Z, method = "mahalanobis")
#       outliers_Online <- detectionOutliers(distances,cutoff)
#       
#       
#       tc <- table(resultsSimul$labelsVrais[1:(nrow(Z))], as.numeric(outliers_Online)[1:(nrow(Z))])
#       tc
#       tc <- safe_access_tc(tc)
#       if((tc["0","0"] + tc["0","1"]) != 0)
#       {faux_positifs_online[i]   <- round((tc["0", "1"]/(tc["0", "1"] + tc["0", "0"]))*100,2)}
#       auc_online <- round(auc(outliers_Online,resultsSimul$labelsVrais),2)*100
#       
#          
#       #Méthode OGK
#       
#       resOGK <- covOGK(Z, sigmamu = s_mad)
#       SigmaOGK <- res$cov
#       medOGK <- res$center
#       rmse_Sigma_ogk <- norm(SigmaOGK - Sigma1,"F")
#       
#       distances <- calcule_vecteur_distances(Z,med,SigmaOGK)
#       #cutoff = calcule_cutoff(distances,d)
#       #cutoff =  qchisq(p = 0.95, df = ncol(Z))
#       #outliers_listMaha <- check_outliers(Z, method = "mahalanobis")
#       outliers_OGK <- detectionOutliers(distances,cutoff)
#       
#       
#       tc <- table(resultsSimul$labelsVrais[1:(nrow(Z) - 1)], as.numeric(outliers_OGK)[1:(nrow(Z) - 1)])
#       tc
#       tc <- safe_access_tc(tc)
#       if((tc["0","0"] + tc["0","1"]) != 0)
#       {faux_positifs_ogk[i]   <- round((tc["0", "1"]/(tc["0", "1"] + tc["0", "0"]))*100,2)}
#       auc_ogk <- auc(outliers_OGK,resultsSimul$labelsVrais)
#       
#       #Comédiane
#       
#       resComed <-  covComed(Z)
#       SigmaComed <- resComed$cov
#       medComed <- resComed$center
#       distances <- calcule_vecteur_distances(Z, medComed, SigmaComed)
#       
#       rmse_Sigma_comed[j,i] <- norm(SigmaComed - Sigma1,"F")
#       rmse_med_comed[j,i] <-  sqrt(sum((mu1 - medComed)^2))
#       
#       
#       outliers_Comed <- detectionOutliers(distances,cutoff)
#       
#       
#       tc <- table(resultsSimul$labelsVrais[1:(nrow(Z) - 1)], as.numeric(outliers_Comed)[1:(nrow(Z) - 1)])
#       tc
#       tc <- safe_access_tc(tc)
#       if((tc["0","0"] + tc["0","1"]) != 0)
#       {faux_positifs_comed[j,i]   <- round((tc["0", "1"]/(tc["0", "1"] + tc["0", "0"]))*100,2)}
#       auc_comed[j,i] <- auc(outliers_Comed,resultsSimul$labelsVrais)
#       
#       
#       #Covariance empirique
#       
#       moyEmp <- as.numeric(colMeans(Z))
#       covEmp <- cov(Z)
#       
#       distances <- calcule_vecteur_distances(Z,moyEmp, covEmp)
#       
#       outliers_Cov <- detectionOutliers(distances,cutoff = cutoff)
#       
#       tc <- table(resultsSimul$labelsVrais[1:(nrow(Z) - 1)], as.numeric(outliers_Cov)[1:(nrow(Z) - 1)])
#       tc
#       tc <- safe_access_tc(tc)
#       if((tc["0","0"] + tc["0","1"]) != 0)
#       {faux_positifs_covEmp[j,i]   <- round((tc["0", "1"]/(tc["0", "1"] + tc["0", "0"]))*100,2)}
#       auc_covEmp[j,i] <- auc(outliers_Cov,resultsSimul$labelsVrais)
#       
#       
#     }
#     
#       
#       
#     }
#   
#   results_data <- data.frame(
#     #Taux_Contamination = taux_contamination,
#     rmse_Sigma_ogk = rmse_Sigma_ogk[1,],
#     rmse_med_ogk = rmse_med_ogk[1,],
#     FP_OGK = faux_positifs_ogk[1,],
#     auc_ogk = auc_ogk[1,],
#     rmse_Sigma_comed = rmse_Sigma_comed[1,],
#     rmse_med_comed = rmse_med_comed[1,],
#     FP_Comed = faux_positifs_comed[1,],
#     auc_comed = auc_comed[1,],
#     rmse_Sigma_shrink = rmse_Sigma_shrink[1,],
#     rmse_med_shrink = rmse_med_shrink[1,],
#     FP_Shrink = faux_positifs_shrink[1,],
#     auc_shrink = auc_shrink[1,],
#     rmse_Sigma_Cov = rmse_Sigma_covEmp[1,],
#     rmse_med_covEmp = rmse_med_covEmp[1,],
#     FP_covEmp = faux_positifs_covEmp[1,],
#     auc_covEmp = auc_covEmp[1,],
#     rmse_Sigma_Offline =  rmse_Sigma_offline[1,] ,
#     rmse_med_offline = rmse_med_offline[1,],
#     auc_offline = auc_offline[1,],
#     FP_Offline = faux_positifs_offline[1,],
#     RMSE_Sigma_Online =  rmse_Sigma_online[1,] ,
#     rmse_med_online = rmse_med_online[1,],
#     auc_online = auc_online[1,],
#     FP_Online = faux_positifs_online[1,],
#     RMSE_Sigma_Streaming =  rmse_Sigma_streaming[1,] ,
#     rmse_med_streaming = rmse_med_streaming[1,],
#     auc_streaming = auc_streaming[1,],
#     FP_Streaming = faux_positifs_streaming[1,]
#     
#     
#     
#   )
#   
#   
#   row.names(results_data) <- taux_contamination
#   
#   return(results_data)
#     
#   }
#   




######Boucle calcul outliers#####


calcule_outliers <- function(n = 1e4, d = 10, c = sqrt(d),rho = 0.8, mu1 = rep(0,d),mu2 = 5*rep(1,d),Sigma1 = creerMatriceToeplitz(rho,d) ,Sigma2 = permuterLignesColonnes(Sigma1,lignes_a_permuter = c(1,2)),colonnes_a_permuter = c(1,2),cutoff = qchisq(p = 0.95,df = ncol(Z)),contamin = "moyenne",depart = 100)  
{

taux_contamination <- c(0,2, 5, 10, 15, 20, 25, 30, 40)

faux_positifs_covEmp <- numeric(length(taux_contamination))
faux_negatifs_covEmp <- numeric(length(taux_contamination))
faux_positifs_maha <- numeric(length(taux_contamination))
faux_negatifs_maha <- numeric(length(taux_contamination))
faux_positifs_EPP <- numeric(length(taux_contamination))
faux_negatifs_EPP <- numeric(length(taux_contamination))
faux_positifs_online <- numeric(length(taux_contamination))
faux_negatifs_online <- numeric(length(taux_contamination))
faux_positifs_offline <- numeric(length(taux_contamination))
faux_negatifs_offline <- numeric(length(taux_contamination))
faux_positifs_streaming <- numeric(length(taux_contamination))
faux_negatifs_streaming <- numeric(length(taux_contamination))
faux_positifs_comed <- numeric(length(taux_contamination))
faux_negatifs_comed <- numeric(length(taux_contamination))
faux_positifs_shrink <- numeric(length(taux_contamination))
faux_negatifs_shrink <- numeric(length(taux_contamination))



# Vecteurs pour stocker les temps de calcul
temps_maha <- numeric(length(taux_contamination))
temps_EPP <- numeric(length(taux_contamination))
temps_online <- numeric(length(taux_contamination))
temps_offline <- numeric(length(taux_contamination))
temps_streaming <- numeric(length(taux_contamination))
temps_comed <- numeric(length(taux_contamination))
temps_shrink <- numeric(length(taux_contamination))
temps_covEmp <- numeric(length(taux_contamination))

#Vecteur pour stocker les distances

distances <- rep(0,n)

#erreursSigma <-  array(0, dim = c(n, length(taux_contamination), 10))

for (i in seq_along(taux_contamination)) {
  delta <- taux_contamination[i]
  #delta <- 0
  #contamin = "moyenne"
  print(contamin)
  p1 <- 1 - delta / 100
  
  p2 <- 1 - p1
  #mu1
  resultsSimul <- genererEchantillon(n, d, mu1, mu2, p1, p2, Sigma1 = Sigma1, Sigma2 = Sigma2,contamin = contamin)
  Z <- resultsSimul$Z
  # Temps pour la cov empirique
  
  temps_covEmp[i] <- system.time({
    #distances <- calcule_vecteur_distancesEmpirique(Z)
    #distances <- calcule_vecteur_distancesEmpiriqueVrai(Z,Sigma1)
    #distances
    distances <- calcule_vecteur_distances(Z,as.numeric(colMeans(Z)), cov(Z))
    
    outliers <- detectionOutliers(distances,cutoff = cutoff)
    print("ok Cov")
    print(delta)
    
    tc <- table(resultsSimul$labelsVrais[1:(nrow(Z) -1)], as.numeric(outliers)[1:(nrow(Z) - 1)])
    tc <- safe_access_tc(tc)
    
    
    if((tc["0","0"] + tc["0","1"]) != 0)
    {faux_positifs_covEmp[i]   <- round((tc["0", "1"]/(tc["0", "1"] + tc["0", "0"]))*100,2)}
    if((tc["1","0"] + tc["1","1"]) != 0)
    {faux_negatifs_covEmp[i] <- round((tc["1", "0"]/(tc["1","0"] + tc["1","1"]))*100,2)}
    #faux_negatifs_maha[i]
    #tc
  })["elapsed"]
  
  
  
  
  # Temps pour la méthode OGK
  temps_maha[i] <- system.time({
    
    
    res <- covOGK(Z, sigmamu = s_mad)
    Sigma <- res$cov
    med <- res$center
    distances <- calcule_vecteur_distances(Z,med,Sigma)
    #cutoff = calcule_cutoff(distances,d)
    #cutoff =  qchisq(p = 0.95, df = ncol(Z))
    #outliers_listMaha <- check_outliers(Z, method = "mahalanobis")
    outliers_listMaha <- detectionOutliers(distances,cutoff)
    print("ok OGK")
    print(delta)
    
    tc <- table(resultsSimul$labelsVrais[1:(nrow(Z) - 1)], as.numeric(outliers_listMaha)[1:(nrow(Z) - 1)])
    tc
    tc <- safe_access_tc(tc)
    if((tc["0","0"] + tc["0","1"]) != 0)
    {faux_positifs_maha[i]   <- round((tc["0", "1"]/(tc["0", "1"] + tc["0", "0"]))*100,2)}
    #faux_positifs_maha[i]
    #else {faux_positifs_maha[i]   <- 0}
    if((tc["1","0"] + tc["1","1"]) != 0)
    {faux_negatifs_maha[i] <- round((tc["1", "0"]/(tc["1","0"] + tc["1","1"]))*100,2)}
    #else {faux_negatifs_maha[i] <- tc["1", "0"]}
    #faux_negatifs_maha[i]
    #faux_negatifs_maha[i]
    #tc
  })["elapsed"]
  
  # Temps pour la méthode EPP
  #temps_EPP[i] <- system.time({
  #res.KurtM.Tribe <- EPPlab(Z, PPalg = "GA", n.simu = 100, maxiter = 1000, sphere = TRUE)
  #OUTms <- EPPlabOutlier(res.KurtM.Tribe, k = qchisq(p = 0.975, df = ncol(Z)), location = median, scale = sd)
  #outliersEPP <- OUTms$outlier
  #tc <- table_contingence(resultsSimul$labelsVrais[1:9999], as.numeric(outliersEPP)[1:9999])
  #tc
  #if(i == 1){faux_negatifs_EPP[i] <- 0 }else{faux_negatifs_EPP[i] <- tc["1", "0"]}
  #faux_positifs_EPP[i] <- tc["0", "1"]
  #if ("0" %in% rownames(tc) && "1" %in% colnames(tc)) {
  #  faux_positifs_EPP[i]   <- tc["0", "1"]
  #}
  #else {faux_positifs_maha[i]   <- 0}
  #if ("1" %in% rownames(tc) && "0" %in% colnames(tc)) {faux_negatifs_EPP[i] <- tc["1","0"]}
  
  #  })["elapsed"]

#  Temps pour la méthode Online
  temps_online[i] <- system.time({
    #params <- initialiser_parametres(
      #Y = Z,
     # c = sqrt(10),
      #Sigma = Sigma1,
      #r = 1.5,
      #k = 1
    #)
    #depart = 100
   results <- estimation(Z,depart = depart,methodeEstimation = "online")
    #miter <- results$miter
    #U <- results$U
    #lambda <- results$lambdaIter
    distances <- results$distances
    #c <- calcule_cutoff(distances,d)
    outliers_labels <- detectionOutliers(distances,cutoff =  cutoff)
    print("ok online")
    print(delta)
    
    tc <- table(resultsSimul$labelsVrais[1:(nrow(Z) - 1)], outliers_labels[1:(nrow(Z) - 1)])
    tc <- safe_access_tc(tc)

    if((tc["0","0"] + tc["0","1"]) != 0)
    {faux_positifs_online[i]   <- round((tc["0", "1"]/(tc["0","1"] + tc["0","0"]))*100,2)}
    else {faux_positifs_maha[i]   <- 0}
    if((tc["1","0"] + tc["1","1"]) != 0){
      faux_negatifs_online[i] <- round((tc["1", "0"]/(tc["1","0"] + tc["1","1"]))*100,2)
    }
  })["elapsed"]
  # 
  # Temps pour la méthode Offline
  temps_offline[i] <- system.time({
    
    resultsOffline <- estimation(Z,methodeEstimation = "offline")
    
    m <- resultsOffline$med
    SigmaEstim <- resultsOffline$SigmaOffline
    #distances <- calcule_vecteur_distances(Z,m,SigmaEstim)
    
    #cutoff <- calcule_cutoff(distances,d)
    outliers_labels <- detectionOutliers(resultsOffline$distances,cutoff =  cutoff)
    print("ok offline")
    print(delta)
    
    #auc <- courbeROCAUC (resultsSimul$labelsVrais,outliers_labels)
    tc <- table(resultsSimul$labelsVrais[1:(nrow(Z) - 1)], as.numeric(outliers_labels[1:(nrow(Z) - 1)]))
    tc <- safe_access_tc(tc)
    #tc
    #print(i)
    if((tc["0","1"] + tc["0","0"])!= 0)
    {faux_positifs_offline[i]   <- round((tc["0", "1"]/(tc["0","1"] + tc["0","0"]))*100,2)}
    #else {faux_positifs_maha[i]   <- 0}
    if((tc["1","0"] + tc["1","1"]) != 0){
      faux_negatifs_offline[i] <-  round((tc["1","0"]/(tc["1","0"] + tc["1","1"]))*100,2)
    }
  })["elapsed"]
  
  
  # Temps pour la méthode streaming
 temps_streaming[i] <- system.time({
    
    resultsStr <- estimation(Z,methodeEstimation = "streaming")
    
    m <- resultsStr$med
    SigmaEstim <- resultsStr$SigmaOffline
    #distances <- calcule_vecteur_distances(Z,m,SigmaEstim)
    distances <- resultsStr$distances
    
    #cutoff <- calcule_cutoff(distances,d)
    outliers_labels <- detectionOutliers(resultsStr$distances,cutoff =  cutoff)
    print("ok streaming")
    print(delta)
    
    tc <- table(resultsSimul$labelsVrais[1:(nrow(Z) - 1)], as.numeric(outliers_labels[1:(nrow(Z) - 1)]))
    tc <- safe_access_tc(tc)
    #tc
    #print(i)
    if((tc["0","1"] + tc["0","0"])!= 0)
{faux_positifs_streaming[i]   <- round((tc["0", "1"]/(tc["0","1"] + tc["0","0"]))*100,2)}
else {faux_positifs_maha[i]   <- 0}
if((tc["1","0"] + tc["1","1"]) != 0){
faux_negatifs_streaming[i] <-  round((tc["1","0"]/(tc["1","0"] + tc["1","1"]))*100,2)
}
 })["elapsed"]

  
  # Temps pour la comédiane
  temps_comed[i] <- system.time({
    med <- covComed(Z)$center
    comed <- covComed(Z)$cov
    length(med)
    dim(comed)
    #(Z[1,] - med) %*% solve(comed) %*% (t(Z[1,] - med))
    distances <- calcule_vecteur_distances(Z,med,comed)
    #cutoff <- calcule_cutoff(distances,d)
    cutoff <- qchisq(p = 0.95, df = ncol(Z))
    outliers_labels <- detectionOutliers(distances,cutoff)
    print("ok comed")
    print(delta)
    
    tc <- table(resultsSimul$labelsVrais[1:(nrow(Z) - 1)], as.numeric(outliers_labels[1:(nrow(Z)- 1)]))
    tc <- safe_access_tc(tc)
    
    #print(i)
    if((tc["0","1"] + tc["0","0"])!= 0)
    {faux_positifs_comed[i]   <- round((tc["0", "1"]/(tc["0","1"] + tc["0","0"]))*100,2)}
    #else {faux_positifs_maha[i]   <- 0}
    if((tc["1","1"] + tc["1","0"])!= 0)
    {faux_negatifs_comed[i] <- round((tc["1","0"]/(tc["1","0"] + tc["1","1"]))*100,2)}
    
    
  })["elapsed"]
  
  
  # Temps pour la comédiane
  temps_shrink[i] <- system.time({
    med <- covComed(Z)$center
    SigmaShrink <- covCor(Z)
    #length(med)
    #dim(comed)
    #(Z[1,] - med) %*% solve(comed) %*% (t(Z[1,] - med))
    distances <- calcule_vecteur_distances(Z,med,SigmaShrink)
    #cutoff <- calcule_cutoff(distances,d)
    cutoff <- qchisq(p = 0.95, df = ncol(Z))
    outliers_labels <- detectionOutliers(distances,cutoff)
    print("ok shrink")
    print(delta)
    
    tc <- table(resultsSimul$labelsVrais[1:(nrow(Z) - 1)], as.numeric(outliers_labels[1:(nrow(Z) - 1)]))
    #tc
    tc <- safe_access_tc(tc)
    
    if((tc["0","1"] + tc["0","0"])!= 0)                   
    {faux_positifs_shrink[i]   <- round((tc["0", "1"]/(tc["0","1"] + tc["0","0"]))*100,2)}
    #else {faux_positifs_maha[i] if   <- 0}
    if((tc["1","0"] + tc["1","1"])!= 0){
      faux_negatifs_shrink[i] <- round((tc["1","0"]/(tc["1","0"] + tc["1","1"]))*100,2)
    }  })["elapsed"]
  
  
}


results_outliers <- data.frame(
  #Taux_Contamination = taux_contamination,
  FN_Cov = faux_negatifs_covEmp,
  FP_Cov= faux_positifs_covEmp,
  Tps_Cov = temps_covEmp,
  FN_OGK = faux_negatifs_maha,
  FP_OGK = faux_positifs_maha,
  Tps_OGK = temps_maha,
  #FP_EPP = faux_positifs_EPP,
  #FN_EPP = faux_negatifs_EPP,
  #Temps_EPP = temps_EPP,
  FN_Online = faux_negatifs_online,
  FP_Online = faux_positifs_online,
  Tps_Online = temps_online,
  FN_Offline = faux_negatifs_offline,
  FP_Offline = faux_positifs_offline,
  Tps_Offline = temps_offline,
  FN_Streaming = faux_negatifs_streaming,
  FP_Streaming = faux_positifs_streaming,
  Tps_Streaming = temps_streaming,
  FN_Comed = faux_negatifs_comed,
  FP_Comed = faux_positifs_comed,
  Tps_Comed = temps_comed,
  FN_Shrink = faux_negatifs_shrink,
  FP_Shrink = faux_positifs_shrink,
  Tps_Shrink = temps_shrink
)


row.names(results_outliers) <- taux_contamination

return(results_outliers)
}

