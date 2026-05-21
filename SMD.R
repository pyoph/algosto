remove_high_corr <- function(Z, threshold = 0.9999)
{
  C <- cor(Z)
  
  C[lower.tri(C, diag = TRUE)] <- 0
  
  remove <- which(apply(abs(C) > threshold, 2, any))
  
  if(length(remove) > 0)
  {
    cat("Removing cols:", remove, "\n")
  }
  
  Z[, -remove, drop = FALSE]
return(Z)
  }


url <- "https://github.com/NetManAIOps/OmniAnomaly/archive/refs/heads/master.zip"
dest <- "smd.zip"
setwd("~")

download.file(url, destfile = dest, mode = "wb")

unzip(dest)
# chercher tous les fichiers SMD
files <- list.files(".", recursive = TRUE, full.names = TRUE)

test_files  <- grep("test/machine", files, value = TRUE)
label_files <- grep("test_label/machine", files, value = TRUE)

length(test_files)
length(label_files)


sequential_clean_columns <- function(Z, begin = 1, cut = NULL) {
  
  Z <- as.matrix(Z)
  p <- ncol(Z)
  
  if (is.null(cut)) cut <- 0.95
  
  start <- NULL
  
  # =========================
  # 1) trouver prefix stable
  # =========================
  
  for (k in (begin + 1):p) {
    
    ok <- tryCatch({
      
      onlineRobustVariance(
        Z[, begin:k, drop = FALSE],
        computeOutliers = TRUE,
        cutoff = cut,
        cutinit = 0.6,
        nDataInit = 200,
        c_m = 2,
        batch = 3
      )
      
      TRUE
      
    }, error = function(e) {
      FALSE
    })
    
    if (ok) {
      start <- k
      cat("✔ first stable prefix found at:", k, "\n")
      break
    }
  }
  
  if (is.null(start)) {
    stop("No stable prefix found")
  }
  
  keep <- c(begin, start)
  
  for (j in (start + 1):p) {
    
    cat("j =", j, "keep =", keep, "\n")
    
    ok <- tryCatch({
      
      onlineRobustVariance(
        Z[, c(keep, j), drop = FALSE],
        computeOutliers = TRUE,
        cutoff = cut,
        cutinit = 0.6,
        nDataInit = 200,
        c_m = 2,
        batch = 3
      )
      
      TRUE
      
    }, error = function(e) FALSE)
    
    if (ok) {
      keep <- c(keep, j)
      cat("✔ col:", j, "OK\n")
    } else {
      cat("✘ col:", j, "KO\n")
    }
  }
  
  return(list(
    Z = Z[, keep, drop = FALSE],
    keep = keep
  ))
}


load_smd_machine <- function(test_path, label_path) {
  
  Z <- as.matrix(
    read.table(test_path, sep = ",")
  )
  
  labels <- scan(label_path)
  
  list(
    Z = Z,
    labels = labels,
    name = basename(test_path)
  )
}

smd_data <- list()

for (i in seq_along(test_files)) {
  
  smd_data[[i]] <- load_smd_machine(
    test_files[i],
    label_files[i]
  )
}

nbmachines = length(smd_data)

nb_anomalies <- sapply(
  smd_data,
  function(x) sum(x$labels == 1)
)

nb_total <- sapply(
  smd_data,
  function(x) length(x$labels)
)

ratio_anomalies <- 100*round(nb_anomalies / nb_total,2)

erreursSigmaSMD = array(0,dim = c(6,nbmachines))
faux_positifsSMD= array(0,dim = c(7,nbmachines))
faux_negatifsSMD = array(0,dim = c(7,nbmachines))
outliersLabelsSMD <- vector("list", nbmachines)
temps = array(0,dim=c(7))
#keep_columns_ok = list()

for(j in 1:(nbmachines)){
  
  cat("\n====================\n")
  cat("Machine", j, "\n")
  cat("====================\n")
  
  
Z <- smd_data[[j]]$Z
labels <- smd_data[[j]]$labels
Z <- as.matrix(Z)

outliersLabelsSMD[[j]] <- matrix(
  0,
  nrow = nrow(Z),
  ncol = 7
)

#Enlèvement des colonnes qui posent problème d'inversibilité
#Z <- Z[, !(colnames(Z) %in% c("V5","V8","V9","V10","V11", "V13","V17","V18","V23","V25","V27", "V29", "V30","V32","V33","V35","V36","V37","V38")), drop = FALSE]

# 
# #Filtre 1 : enlever les colonnes qui font planter onlineRobustVariance
 res <- tryCatch({
#   
   clean_cols <- sequential_clean_columns(Z)
#   
  list(
     Z = clean_cols$Z,
     keep = clean_cols$keep
   )
#   
 }, error = function(e) {
#   
   cat("Erreur détectée -> retry sans colonne 1\n")
#   
   Z2 <- Z[, -1, drop = FALSE]
#   
   clean_cols <- sequential_clean_columns(Z2)
#   
   list(
     Z = clean_cols$Z,
     keep = clean_cols$keep + 1
   )
 })
# 
 Z <- res$Z
 keep_columns_ok[[j]] <- res$keep
# 
# 
# 
# #Filtre 2 : enlever les colonnes avec MAD ~ 0
# 
# 
tol = 1e-3
# 
mad_cols <- apply(Z, 2, mad, na.rm = TRUE)
# 
keep <- which(mad_cols > tol)
# 
 removed <- setdiff(1:ncol(Z), keep)
# 
 if (length(removed) > 0) {
   cat("Colonnes retirées (MAD ~ 0):", removed, "\n")
 }

keep_columns_ok[[j]] = keep



Z = Z[,keep_columns_ok[[j]],drop = FALSE]

Z = remove_high_corr(Z,threshold = 0.99)
Z_clean = Z[labels == 0,]

###########Extraction moyenne et covariance sans outliers
meanTrueCov = colMeans(Z_clean)

Sigma_trueCov = cov(Z_clean)

Sigma_trueMCD = covMcd(Z_clean)$cov

Sigma_trueOGK = covOGK(Z_clean)$cov

Sigma_trueOffl = offlineRobustVariance(Z_clean)$variance

Sigma_trueOnl = onlineRobustVariance(Z_clean,batch = 1)$variance

Sigma_trueStrm = onlineRobustVariance(Z_clean)$variance


resSamplecov= SampleCovOnline(Z,quantcutoff = TRUE,nDataInit = 1e3 )

erreursSigmaSMD[1,j] = norm(resSamplecov$Sigma - Sigma_trueCov,"F")

table(resSamplecov$outliers_labels,labels)

# resStrm = onlineRobustVariance(Z,batch = ncol(Z),computeOutliers = TRUE)
# resOnline = onlineRobustVariance(Z,batch = 1,computeOutliers = TRUE)
# resOffline = offlineRobustVariance(Z,computeOutliers = TRUE)
resmcd = covMcd(Z)

erreursSigmaSMD[2,j] = norm(resmcd$cov - Sigma_trueMCD,"F")

#resOfflineClean = offlineRobustVariance(Z_clean,computeOutliers = TRUE)

resogk = covOGK(Z, sigmamu = scaleTau2)

erreursSigmaSMD[3,j] = norm(resogk$cov - Sigma_trueOGK,"F")

distogk = rep(0,nrow(Z))
distoffl = rep(0,nrow(Z))
distmcd = rep(0,nrow(Z))
invSigmaOGK = inv_safe(resogk$cov)
invSigmaMCD = inv_safe(resmcd$cov)
invSigmaTrue = inv_safe(Sigma_true)
#med = t(WeiszfeldMedian(Z))
#quantemp = 0
outliersLabelsSMD[[j]][, 1] <- resSamplecov$outliers_labels
distrue = rep(0,nrow(Z_clean))

for(i in 1:nrow(Z))
{
  distmcd[i] = t(Z[i,] - (resmcd$center))%*%invSigmaMCD%*%(Z[i,] - (resmcd$center))
  distogk[i] = t(Z[i,] - resogk$center)%*%invSigmaOGK%*%(Z[i,] - resogk$center)
  distrue[i] = t(Z[i,] - meanTrue)%*%invSigmaTrue%*%(Z[i,] - meanTrue)
  
}

#############################################Quantile oracle###############################
quantoracle = quantile(distrue,.95)
###########################################################################################

#oracle = rep(0,nrow(Z))

#quantemp = median(distoffl)
#pred = rep(0,nrow(Z))

cutoff_true = calcule_cutoff(distrue,c_m = 2,cutinit=0.6,cutoff = 0.9,n = nrow(Z))
cut_off_ogk = calcule_cutoff(distogk,c_m = 2,cutinit=0.6,cutoff = 0.9,n = nrow(Z))
cut_off_mcd = calcule_cutoff(distmcd,c_m = 2,cutinit=0.6,cutoff = 0.9,n = nrow(Z))


for (i in 1:nrow(Z)) {
  
  if (distrue[i] > cutoff_true) {
    outliersLabelsSMD[[j]][i, 7] <- 1
  }
  
  if (distogk[i] > cut_off_ogk) {
    outliersLabelsSMD[[j]][i, 3] <- 1
  }
  
  if (distmcd[i] > cut_off_mcd) {
    outliersLabelsSMD[[j]][i, 2] <- 1
  }
  
}
#pred = distoffl > quantemp

#table(pred,labels)

#table(oracle,labels)
cut=0.9
res0 <- offlineRobustVariance(Z,computeOutliers = TRUE,cutoff=cut,cutinit=0.6,c_m=1.5)

erreursSigmaSMD[4,j] = norm(res0$variance - Sigma_trueOffl,"F")
outliersLabelsSMD[[j]][, 4] = res0$outlier_labels


table(res0$outlier_labels,labels)
resStrm <- onlineRobustVariance(Z,computeOutliers = TRUE,cutoff=cut,cutinit=0.6,nDataInit = 200,c_m=2,batch = floor(nrow(Z)/ncol(Z)))
erreursSigmaSMD[5,j] = norm(resStrm$variance - Sigma_true,"F")
outliersLabelsSMD[[j]][, 5] = resStrm$outlier_labels
table(resStrm$outlier_labels,labels)

resOnline <- onlineRobustVariance(Z,computeOutliers = TRUE,cutoff=cut,cutinit=0.6,nDataInit = 200,c_m=2,batch = 1)
outliersLabelsSMD[[j]][, 6] = resOnline$outlier_labels

erreursSigmaSMD[6,j] = norm(resOnline$variance - Sigma_trueOnl,"F")

table(resOnline$outlier_labels,labels)


}
setwd("~")


file <- paste0("boxploterreursSigma_SMD", ".pdf")

pdf(file, width = 18, height = 6) 

boxplot(
  erreursSigmaSMD[1,],erreursSigmaSMD[2,],erreursSigmaSMD[3,],erreursSigmaSMD[4,],erreursSigmaSMD[5,],erreursSigmaSMD[6,],
  log = "y",                     # échelle logarithmique sur Y
  names = c("sample covariance online","mcd","ogk","offline","online","streaming"),
  main = "",
  ylab = "",
  col = rainbow(ncol(erreursSigmaSMD))
)
dev.off()

#############################ROC#######################################


compute_roc <- function(score, labels) {
  
  thresholds <- sort(unique(score))
  
  TPR <- numeric(length(thresholds))
  FPR <- numeric(length(thresholds))
  
  for (i in seq_along(thresholds)) {
    
    th <- thresholds[i]
    pred <- ifelse(score >= th, 1, 0)
    
    TP <- sum(pred == 1 & labels == 1)
    TN <- sum(pred == 0 & labels == 0)
    FP <- sum(pred == 1 & labels == 0)
    FN <- sum(pred == 0 & labels == 1)
    
    TPR[i] <- TP / (TP + FN + 1e-12)
    FPR[i] <- FP / (FP + TN + 1e-12)
  }
  
  ord <- order(FPR, TPR)
  
  list(
    FPR = FPR[ord],
    TPR = TPR[ord],
    auc = sum(diff(c(0, FPR[ord])) * (head(TPR[ord], -1) + tail(TPR[ord], -1))/2)
  )
  
  
}

roc_samplecov <- compute_roc(resSamplecov$distances,labels)
roc_strm    <- compute_roc(resStrm$distances, labels)
roc_online  <- compute_roc(resOnline$distances, labels)
roc_offline <- compute_roc(resOffline$distances, labels)
roc_mcd     <- compute_roc(distmcd, labels)
roc_ogk     <- compute_roc(distogk, labels)


file <- paste0("roc", ".pdf")

pdf(file, width = 18, height = 6) 


plot(roc_strm$FPR, roc_strm$TPR,
     type = "l", col = "red", lwd = 2,
     xlim = c(0,1), ylim = c(0,1),
     xlab = "FPR", ylab = "TPR",
     main = "",cex = 2.5)
lines(roc_samplecov$FPR, roc_samplecov$TPR, col = "darkgreen", lwd = 2)
lines(roc_online$FPR, roc_online$TPR, col = "blue", lwd = 2)
lines(roc_offline$FPR, roc_offline$TPR, col = "orange", lwd = 2)
lines(roc_mcd$FPR, roc_mcd$TPR, col = "black", lwd = 2)
lines(roc_ogk$FPR, roc_ogk$TPR, col = "brown", lwd = 2)

abline(0,1,lty=2)

dev.off()
save(roc_samplecov,roc_strm,roc_online,roc_offline,roc_mcd,roc_ogk,file = "roc.RData")


#######################TRUE LABELS FALSE NEGATIVES AND FALSE POSITVES AND THRESHOLDS TRAJECTORIES ####################

compute_rates <- function(pred, labels) {

  n <- length(pred)
  
  FP_rate <- rep(0, n)
  FN_rate <- rep(0, n)
  fp <- 0
  tn <- 0
  tp <- 0
  fn <- 0
  
  for (t in 1:n) {
    
    if (pred[t] == 1 && labels[t] == 0) fp <- fp + 1
    if (pred[t] == 0 && labels[t] == 1) fn <- fn + 1
    if (pred[t] == 1 && labels[t] == 1) tp <- tp + 1
    if (pred[t] == 0 && labels[t] == 0) tn <- tn + 1
    
    # vrai taux streaming
    FP_rate[t] <- fp / (fp + tn + 1e-10)
    FN_rate[t] <- fn / (fn + tp + 1e-10)
  }
  
  list(FP_rate = FP_rate,
       FN_rate = FN_rate)
}

####Calcul des trajectoires pour les 3 méthodes online et l'oracle######
for(j in 1:nbmachines){
outlmach = outliersLabelsSMD[[j]]
labels = smd_data[[j]]$labels

rates_samplecov = compute_rates(outlmach[,1], labels)
rates_online  <- compute_rates(outlmach[,5], labels)
rates_strm <- compute_rates(outlmach[,6], labels)
rates_oracle <- compute_rates(outlmach[,7], labels)

##################################Trajectoire des seuils################################"
# 
# distonl = resOnline$distances
# 
# distStrm = resStrm$distances
# 
# scal_factor_onl = rep(0,nrow(Z))
# scal_factor_str = rep(0,nrow(Z))
# 
# for (i in 1:nrow(Z)){
#   scal_factor_onl[i] = qchisq(.5,df = ncol(Z))/median(distonl[1:i])
#   scal_factor_str[i] = qchisq(.5,df = ncol(Z))/median(distStrm[1:i])
#   
# }
setwd("~/figures")
nom_fichier = paste0("trajectories_SMD_mach-",j,".pdf")
pdf(nom_fichier, width = 14, height = 10)

par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

x_vals = 1:length(rates_strm$FN_rate)


# =====================================================
# 1. LABELS
# =====================================================

plot(x_vals, labels,
     col = "purple",
     xlab = "", ylab = "",
     yaxt = "n", xaxt = "n",
     ylim = c(0, 1),
     main = "Ground truth"
)

axis(2, at = c(0, 1), las = 1, cex.axis = 1.2)
axis(1, at = seq(1000, max(x_vals), by = 1000),
     las = 1, cex.axis = 1.2)
box()

# =====================================================
# 2. SCALE FACTORS
# =====================================================
# 
# plot(x_vals, scal_factor_onl,
#      type = "l", lwd = 3, col = "blue",
#      xlab = "", ylab = "",
#      yaxt = "n", xaxt = "n",
#      main = "Scale factors",
#      ylim = range(c(scal_factor_onl, scal_factor_str))
# )
# 
# lines(x_vals, scal_factor_str,
#       col = "red", lwd = 3)
# 
# axis(2, las = 1, cex.axis = 1.2)
# axis(1, at = seq(1000, max(x_vals), by = 1000),
#      las = 1, cex.axis = 1.2)
# box()

# =====================================================
# 2. FALSE NEGATIVE RATE
# =====================================================

plot(x_vals, rates_strm$FN_rate*100,
     type = "l", lwd = 3, col = "red",
     xlab = "", ylab = "",
     yaxt = "n", xaxt = "n",
     main = "False Negative Rate",
     ylim = c(0,100)
)

lines(x_vals, rates_samplecov$FN_rate*100,
      lty = "dotted", col = "darkgreen", lwd = 3)

lines(x_vals, rates_online$FN_rate*100,
      lty = "dashed", col = "blue", lwd = 3)
lines(x_vals, rates_oracle$FN_rate*100,
      lty = "dashed", col = "purple", lwd = 3)

axis(2, las = 1, cex.axis = 1.2)
axis(1, at = seq(1000, max(x_vals), by = 1000),
     las = 1, cex.axis = 1.2)
box()

# =====================================================
# 3. FALSE POSITIVE RATE
# =====================================================

plot(x_vals, rates_strm$FP_rate*100,
     type = "l", lwd = 3, col = "red",
     xlab = "", ylab = "",
     yaxt = "n", xaxt = "n",
     main = "False Positive Rate",
     ylim = range(c(0,20))
)

lines(x_vals, rates_samplecov$FP_rate*100,
      lty = "dotted", col = "darkgreen", lwd = 3)

lines(x_vals, rates_online$FP_rate*100,
      lty = "dashed", col = "blue", lwd = 3)

lines(x_vals, rates_oracle$FP_rate*100,
      lty = "dashed", col = "purple", lwd = 3)

axis(2, las = 1, cex.axis = 1.2)
axis(1, at = seq(1000, max(x_vals), by = 1000),
     las = 1, cex.axis = 1.2)
box()

dev.off()
}
