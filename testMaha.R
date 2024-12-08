#install.packages("easystats")
#install.packages("bigutilsr")
#install.packages("xtable")
#install.packages("RGMM")
library(xtable)
library(reshape2)
library(RobRegression)
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
library("rJava")
library("REPPlab")
library("RGMM")
source("~/algosto/parametres.R")
source("~/algosto/simulations.R")
source("~/algosto/algorithmes.R")
source("~/algosto/resultats.R")
#usethis::create_package("C:/Users/Paul/Documents/codeRTheseLPSM")
#devtools::load_all("C:/Users/Paul/Documents/codeRTheseLPSM")

n <- 1e4

d <- 10

rho <- 0.3

seuil_p_value <- 0.05

Sigma1 <- creerMatriceToeplitz(rho,d)

Sigma2 <- creerMatriceToeplitz(0.7,d)
Sigma1 <- diag(sqrt(1:d))
#Sigma1 <- diag(d)
mu1 <- rep(0,d)

mu2 <- 5*rep(1,d)




p1 <- 0.98
p2 <- 1 - p1


resultsSimul <- genererEchantillon(n, d, mu1, mu2, p1, p2, Sigma1 = Sigma1, Sigma2 = Sigma2)

#cov(Z)

Z <- resultsSimul$Z






Rvar <- RobVar(Z)

m <- Rvar$median
#Z[1,]

Sigma <- Rvar$variance


outliers_labels <- detectionOffline(Z, SigmaEstim = Sigma,m, 0.025)
tc <- table_contingence(resultsSimul$labelsVrais[1:9999], as.numeric(outliers_labels[1:9999]))

tc
resultsSimul$labelsVrais[1:9999]

outliers_labels <- rep(0,n)
#m = rep(0,d)
for (i in 1:nrow(Z)) 
  {
  S <- 0
  
  S <- (Z[i,] - m)%*%solve(Sigma1)%*%t(Z[i,] - m) 
  #print(S)
  # Calcul de la p-value basée sur la statistique du Chi2
  phat <- pchisq(S, df = ncol(Z),lower.tail =FALSE)
  #print(phat)
  # Détection de l'outlier
  if (phat < 0.05) {
    outliers_labels[i] <- 1  # Indiquer qu'il s'agit d'un outlier
    #print("OK")
  } else {
    outliers_labels[i] <- 0  # Indiquer qu'il ne s'agit pas d'un outlier
  }
  
  #if (S > cutoff) {outliers_labels[i] <- 1}
  #else {outliers_labels[i] <- 0}
  #}
}

#dim(t(Z[1,] -m))



indices_labels_1 <- which(outliers_labels[1:9999] == 1)

# Affichage des indices
print(indices_labels_1)

taux_contamination <- c(2, 5, 10, 15, 20, 25, 30, 40)  
faux_positifs_maha <- numeric(length(taux_contamination))
faux_negatifs_maha <- numeric(length(taux_contamination))
faux_positifs_EPP <- numeric(length(taux_contamination))
faux_negatifs_EPP <- numeric(length(taux_contamination))
faux_positifs_online <- numeric(length(taux_contamination))
faux_negatifs_online <- numeric(length(taux_contamination))
faux_positifs_offline <- numeric(length(taux_contamination))
faux_negatifs_offline <- numeric(length(taux_contamination))

# Vecteurs pour stocker les temps de calcul
temps_maha <- numeric(length(taux_contamination))
temps_EPP <- numeric(length(taux_contamination))
temps_online <- numeric(length(taux_contamination))
temps_offline <- numeric(length(taux_contamination))

# Vecteurs pour stocker les temps de calcul
temps_maha <- numeric(length(taux_contamination))
temps_EPP <- numeric(length(taux_contamination))
temps_online <- numeric(length(taux_contamination))
temps_offline <- numeric(length(taux_contamination))

for (i in seq_along(taux_contamination)) {
  delta <- taux_contamination[i]
  
  p1 <- 1 - delta / 100
  p2 <- 1 - p1
  
  resultsSimul <- genererEchantillon(n, d, mu1, mu2, p1, p2, Sigma1 = Sigma1, Sigma2 = Sigma2)
  Z <- resultsSimul$Z
  
  # Temps pour la méthode Mahalanobis
  temps_maha[i] <- system.time({
    outliers_listMaha <- check_outliers(Z, method = "mahalanobis_robust", qchisq(p = 0.975, df = ncol(Z)))
    tc <- table_contingence(resultsSimul$labelsVrais[1:9999], as.numeric(outliers_listMaha)[1:9999])
    tc
    faux_negatifs_maha[i] <- tc["1", "0"]
    faux_positifs_maha[i] <- tc["0", "1"]
  })["elapsed"]
  
  # Temps pour la méthode EPP
  temps_EPP[i] <- system.time({
    res.KurtM.Tribe <- EPPlab(Z, PPalg = "GA", n.simu = 100, maxiter = 1000, sphere = TRUE)
    OUTms <- EPPlabOutlier(res.KurtM.Tribe, k = 1, location = median, scale = sd)
    outliersEPP <- OUTms$outlier
    tc <- table_contingence(resultsSimul$labelsVrais[1:9999], as.numeric(outliersEPP)[1:9999])
    tc
    faux_negatifs_EPP[i] <- tc["1", "0"]
    faux_positifs_EPP[i] <- tc["0", "1"]
  })["elapsed"] 
  
  # Temps pour la méthode Online
  temps_online[i] <- system.time({
    params <- initialiser_parametres(
      Y = Z,
      c = sqrt(10),
      Sigma = Sigma1,
      r = 1.5,
      k = 1
    )
    results <- estimMVOutliers(Z, params$c, params$n, params$d, params$d, params$r, aa = 0.75, niter = 1e4)
    tc <- table_contingence(resultsSimul$labelsVrais[1:9999], results$outlier_labels[1:9999])
    tc
    faux_negatifs_online[i] <- tc["1", "0"]
    faux_positifs_online[i] <- tc["0", "1"]
  })["elapsed"]
  
  # Temps pour la méthode Offline
  temps_offline[i] <- system.time({
    Rvar <- RobVar(Z)
    SigmaEstim <- Rvar$variance
    m <- Rvar$median
    outliers_labels <- detectionOffline(Z, SigmaEstim, m, 0.025)
    tc <- table_contingence(resultsSimul$labelsVrais[1:9999], as.numeric(outliers_labels[1:9999]))
    tc

        faux_negatifs_offline[i] <- tc["1", "0"]
      
        faux_positifs_offline[i] <- tc["0","1"]
      
      
      
  })["elapsed"]
}




results_outliers <- data.frame(
  #Taux_Contamination = taux_contamination,
  FN_Mahalanobis = faux_negatifs_maha,
  FP_Mahalanobis = faux_positifs_maha,
  Temps_Mahalanobis = temps_maha,
  FP_EPP = faux_positifs_EPP,
  FN_EPP = faux_negatifs_EPP,
  Temps_EPP = temps_EPP,
  FN_Online = faux_negatifs_online,
  FP_Online = faux_positifs_online,
  Temps_Online = temps_online,
  FN_Offline = faux_negatifs_offline,
  FP_Offline = faux_positifs_offline,
  Temps_Offline = temps_offline
)

#Tester mahalanobis threshold = 0,05

row.names(results_outliers) <- taux_contamination

results_outliers

table_latex <- xtable(results_outliers, caption = "Faux positifs et faux négatifs pour chaque méthode", label = "tab:results_outliers")

print(table_latex, file = "results_outliersAvecTempsContamStudentVar.tex",digits = 0)

