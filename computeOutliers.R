#install.packages("easystats")
#install.packages("bigutilsr")
#install.packages("xtable")
#install.packages("RGMM")
#install.packages("rrcov")
#install.packages("mvnfast")
library(xtable)
library(reshape2)
library(RobRegression)
library("robustbase")
library(Gmedian)
library(ggplot2)
library(far)
library(gridExtra)
library(microbenchmark)
library("matlib")
library("MASS")
library("corrplot")
library("dplyr")
library("easystats")
library("bigutilsr")
#library("rJava")
#library("REPPlab")
library("RGMM")
#source("~/algosto/parametres.R")
source("~/algosto/simulations.R")
source("~/algosto/algorithmes.R")
source("~/algosto/resultats.R")
source("~/algosto/Outliers.R")

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


######Boucle calcul outliers#####


calcule_outliers <- function(n = 1e4, d = 10, c = sqrt(d),rho = 0.8, distances = rep(0,n), mu1 = rep(0,d),mu2 = 5*rep(1,d),Sigma1 = creerMatriceToeplitz(rho,d) ,Sigma2 = permuterLignesColonnes(Sigma1,lignes_a_permuter = c(1,2),colonnes_a_permuter = c(1,2))
, contamin = "student") 
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
faux_positifs_comed <- numeric(length(taux_contamination))
faux_negatifs_comed <- numeric(length(taux_contamination))
faux_positifs_shrink <- numeric(length(taux_contamination))
faux_negatifs_shrink <- numeric(length(taux_contamination))



# Vecteurs pour stocker les temps de calcul
temps_maha <- numeric(length(taux_contamination))
temps_EPP <- numeric(length(taux_contamination))
temps_online <- numeric(length(taux_contamination))
temps_offline <- numeric(length(taux_contamination))
temps_comed <- numeric(length(taux_contamination))
temps_shrink <- numeric(length(taux_contamination))
temps_covEmp <- numeric(length(taux_contamination))

#Matrice pour stocker les erreurs de calcul de Sigma

#erreursSigma <-  array(0, dim = c(n, length(taux_contamination), 10))

for (i in seq_along(taux_contamination)) {
  delta <- taux_contamination[i]
  #delta <- 2
  p1 <- 1 - delta / 100
  p2 <- 1 - p1
  #mu1
  resultsSimul <- genererEchantillon(n, d, mu1, mu2, p1, p2, Sigma1 = Sigma1, Sigma2 = Sigma2,contamin = "student")
  Z <- resultsSimul$Z
  # Temps pour la cov empirique
  
  temps_covEmp[i] <- system.time({
    #distances <- calcule_vecteur_distancesEmpirique(Z)
    distances <- calcule_vecteur_distancesEmpiriqueVrai(Z,Sigma1)
    #distances
    outliers <- detectionOutliers(distances, nrow(Z),ncol(Z),cutoff = qchisq(p = 0.95,df = ncol(Z)))
    tc <- table_contingence(resultsSimul$labelsVrais[1:9999], as.numeric(outliers)[1:9999])
    tc <- safe_access_tc(tc)
    
    
    if((tc["0","0"] + tc["0","1"]) != 0)
    {faux_positifs_covEmp[i]   <- tc["0", "1"]/(tc["0","0"] + tc["1","0"])}
    if((tc["1","0"] + tc["1","1"]) != 0)
    {faux_negatifs_covEmp[i] <- tc["1","0"]/(tc["1","0"] + tc["1","1"])}
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
    cutoff =  qchisq(p = 0.95, df = ncol(Z))
    #outliers_listMaha <- check_outliers(Z, method = "mahalanobis")
    outliers_listMaha <- detectionOutliers(distances,nrow(Z),ncol(Z), cutoff)
    tc <- table_contingence(resultsSimul$labelsVrais[1:9999], as.numeric(outliers_listMaha)[1:9999])
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
  
  # Temps pour la méthode Online
  temps_online[i] <- system.time({
    #params <- initialiser_parametres(
      #Y = Z,
     # c = sqrt(10),
      #Sigma = Sigma1,
      #r = 1.5,
      #k = 1
    #)
    
    results <- estimMV(Z, c = sqrt(d), nrow(Z), ncol(Z), ncol(Z), r = 1.5, aa = 0.75, niter = 1e4,niterRMon = d ,methode = "eigen")
    miter <- results$miter
    U <- results$U
    lambda <- results$lambdaIter
    distances <- calcule_vecteur_distancesOnline(Z,miter,U,lambda)
    #c <- calcule_cutoff(distances,d)
    outliers_labels <- detectionOutliers(distances, nrow(Z),ncol(Z),cutoff =  qchisq(p = 0.95, df = ncol(Z)))
    tc <- table_contingence(resultsSimul$labelsVrais[1:9999], outliers_labels[1:9999])
    tc <- safe_access_tc(tc)
    
    if((tc["0","0"] + tc["0","1"]) != 0)
    {faux_positifs_online[i]   <- round((tc["0", "1"]/(tc["0","1"] + tc["0","0"]))*100,2)}
    #else {faux_positifs_maha[i]   <- 0}
    if((tc["1","0"] + tc["1","1"]) != 0){
      faux_negatifs_online[i] <- round((tc["1", "0"]/(tc["1","0"] + tc["1","1"]))*100,2)
    }
  })["elapsed"]
  
  # Temps pour la méthode Offline
  temps_offline[i] <- system.time({
    Rvar <- RobVar(Z)
    
    
    m <- Rvar$median
    SigmaEstim <- Rvar$variance
    distances <- calcule_vecteur_distances(Z,m,SigmaEstim)
    #cutoff <- calcule_cutoff(distances,d)
    outliers_labels <- detectionOutliers(distances, nrow(Z),ncol(Z),cutoff =  qchisq(p = 0.95, df = ncol(Z)))
    #outliers_labels <- detectionOffline(Z, SigmaEstim, m, 0.025,cutoff)
    
    tc <- table_contingence(resultsSimul$labelsVrais[1:9999], as.numeric(outliers_labels[1:9999]))
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
    outliers_labels <- detectionOutliers(distances,nrow(Z),ncol(Z),cutoff)
    tc <- table_contingence(resultsSimul$labelsVrais[1:9999], as.numeric(outliers_labels[1:9999]))
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
    outliers_labels <- detectionOutliers(distances,nrow(Z),ncol(Z),cutoff)
    tc <- table_contingence(resultsSimul$labelsVrais[1:9999], as.numeric(outliers_labels[1:9999]))
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
  #FN_Cov = faux_negatifs_covEmp,
  #FP_Cov= faux_positifs_covEmp,
  #Tps_Cov = temps_covEmp,
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
  FN_Comed = faux_negatifs_comed,
  FP_Comed = faux_positifs_comed,
  Tps_Comed = temps_comed,
  FN_Shrink = faux_negatifs_shrink,
  FP_Shrink = faux_positifs_shrink,
  Tps_Shrink = temps_shrink
)

#Tester mahalanobis threshold = 0,05

row.names(results_outliers) <- taux_contamination

return(results_outliers)
}

