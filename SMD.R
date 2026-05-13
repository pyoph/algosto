url <- "https://github.com/NetManAIOps/OmniAnomaly/archive/refs/heads/master.zip"
dest <- "smd.zip"

download.file(url, destfile = dest, mode = "wb")

unzip(dest)

# chercher tous les fichiers SMD
files <- list.files(".", recursive = TRUE, full.names = TRUE)

test_files  <- grep("test/machine", files, value = TRUE)
label_files <- grep("test_label/machine", files, value = TRUE)

length(test_files)
length(label_files)


load_smd_machine <- function(test_path, label_path) {
  
  Z <- as.matrix(read.table(test_path, sep = ","))
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

Z <- smd_data[[1]]$Z
labels <- smd_data[[1]]$labels
Z <- as.matrix(Z)


# =========================================================
# ANALYSE VARIABLES SMD
# =========================================================

# 
Z <- scale(Z)
labels <- as.numeric(labels)

n <- min(nrow(Z), length(labels))
Z <- Z[1:n, ]
labels <- labels[1:n]


effect_size <- sapply(1:ncol(Z), function(j) {
  
  mu1 <- mean(Z[labels == 1, j])
  mu0 <- mean(Z[labels == 0, j])
  sdj <- sd(Z[, j])
  
  abs(mu1 - mu0) / sdj
})


rank1 <- order(effect_size, decreasing = TRUE)


cat("\nTop variables (effect size):\n")
print(rank1[1:10])



###################Extraction de la matrice de covariance sans outliers###########################


#Enlèvement des colonnes qui posent problème d'inversibilité
Z <- Z[, !(colnames(Z) %in% c("V5","V8","V9","V10","V11", "V13","V17","V18","V23","V25","V27", "V29", "V30","V32","V33","V35","V36","V37","V38")), drop = FALSE]
Z_clean = Z[labels == 0,]
Sigma_true = cov(Z_clean)
resSamplecov= SampleCovOnline(Z)
resStrm = onlineRobustVariance(Z,batch = ncol(Z),computeOutliers = TRUE)
resOnline = onlineRobustVariance(Z,batch = 1,computeOutliers = TRUE)
resOffline = offlineRobustVariance(Z,computeOutliers = TRUE)
resmcd = covMcd(Z)

resOfflineClean = offlineRobustVariance(Z_clean,computeOutliers = TRUE)

resogk = covOGK(Z, sigmamu = scaleTau2)
distogk = rep(0,nrow(Z))
distoffl = rep(0,nrow(Z))
distmcd = rep(0,nrow(Z))
invSigmaOGK = solve(resogk$cov)
invSigmaMCD = solve(resmcd$cov)
invSigmaOffl = solve(resOffline$variance)
invSigmaTrue = solve(Sigma_true)
med = t(WeiszfeldMedian(Z))
quantemp = 0
distrue = rep(0,nrow(Z_clean))
meanTrue = colMeans(Z_clean)
#Calcul destinations de Mahalanobis


for(i in 1:nrow(Z))
{
  distmcd[i] = t(Z[i,] - (resmcd$center))%*%invSigmaMCD%*%(Z[i,] - (resmcd$center))
  distogk[i] = t(Z[i,] - resogk$center)%*%invSigmaOGK%*%(Z[i,] - resogk$center)
 distoffl[i] =t(Z[i,] - med)%*%invSigmaOffl%*%(Z[i,] - med)

}
oracle = rep(0,nrow(Z_clean))
quantemp = quantile(distoffl,.95)

table(oracle,labels[labels == 0])

quantemp = median(distoffl)
pred = rep(0,nrow(Z))

disttrue = rep(0,nrow(Z))
for(i in (1:nrow(Z))){
quantemp <- quantemp - (1 + i )^(-0.75) * (as.numeric(distoffl[i] <= quantemp) - 
                                              0.8)
if (distoffl[i] > quantemp){pred[i] = 1}
distrue = t(Z[i,] - meanTrue)%*%invSigmaTrue%*%(Z[i,] - meanTrue)
if(distrue[i] > quantemp) {oracle[i] = 1}
}
#pred = distoffl > quantemp

table(pred,labels)

# =========================================================
# Recherche automatique du meilleur n0
# =========================================================

beta <- 0.75
target_quantile <- 0.9

n0_grid <- seq(1, 5000, by = 10)

results <- data.frame(
  n0 = numeric(),
  TP = numeric(),
  FP = numeric(),
  TN = numeric(),
  FN = numeric(),
  Recall = numeric(),
  Precision = numeric(),
  FPR = numeric(),
  F1 = numeric()
)

for(n0 in n0_grid){
  
  quantemp <- median(distoffl)
  
  pred <- rep(0, length(distoffl))
  
  for(i in 1:length(distoffl)){
    
    gamma <- (i + n0)^(-beta)
    
    quantemp <- quantemp -
      gamma *
      (as.numeric(distoffl[i] <= quantemp) - target_quantile)
    
    if(distoffl[i] > quantemp){
      pred[i] <- 1
    }
  }
  
  # matrice confusion
  TP <- sum(pred == 1 & labels == 1)
  FP <- sum(pred == 1 & labels == 0)
  
  TN <- sum(pred == 0 & labels == 0)
  FN <- sum(pred == 0 & labels == 1)
  
  # métriques
  Recall <- TP / (TP + FN + 1e-12)
  
  Precision <- TP / (TP + FP + 1e-12)
  
  FPR <- FP / (FP + TN + 1e-12)
  
  F1 <- 2 * Precision * Recall /
    (Precision + Recall + 1e-12)
  
  results <- rbind(
    results,
    data.frame(
      n0,
      TP,
      FP,
      TN,
      FN,
      Recall,
      Precision,
      FPR,
      F1
    )
  )
}

# =========================================================
# Meilleur n0 selon F1
# =========================================================

best <- results[which.max(results$F1), ]

print(best)

# =========================================================
# top résultats
# =========================================================

results[order(-results$F1), ][1:10, ]

# 
# library(pROC)
# 
# roc_obj <- roc(res$outlier_labels,distoffl)
# auc_value <- auc(roc_obj)
# 
# cat("AUC :", auc_value, "\n")
# 

# =========================================================
# 14. COURBE ROC
# =========================================================

# plot(roc_obj,
#      xlim = c(0,1),
#      ylim = c(0,1),
#      )

# labels : vecteur 0/1
# score  : score anomalie (ex: Mahalanobis)


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

####Calcul pour les 3 méthodes online######


rates_samplecov = compute_rates(resSamplecov$outliers_labels, labels)
rates_online  <- compute_rates(resOnline$outlier_labels, labels)
rates_strm <- compute_rates(resStrm$outlier_labels, labels)


##################################Trajectoire des seuils################################"

distonl = resOnline$distances

distStrm = resStrm$distances

scal_factor_onl = rep(0,nrow(Z))
scal_factor_str = rep(0,nrow(Z))

for (i in 1:nrow(Z)){
  scal_factor_onl[i] = qchisq(.5,df = ncol(Z))/median(distonl[1:i])
  scal_factor_str[i] = qchisq(.5,df = ncol(Z))/median(distStrm[1:i])
  
}
setwd("~")
pdf("trajectories_SMD.pdf", width = 14, height = 10)

par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

x_vals = 1:nrow(Z)

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

plot(x_vals, scal_factor_onl,
     type = "l", lwd = 3, col = "blue",
     xlab = "", ylab = "",
     yaxt = "n", xaxt = "n",
     main = "Scale factors",
     ylim = range(c(scal_factor_onl, scal_factor_str))
)

lines(x_vals, scal_factor_str,
      col = "red", lwd = 3)

axis(2, las = 1, cex.axis = 1.2)
axis(1, at = seq(1000, max(x_vals), by = 1000),
     las = 1, cex.axis = 1.2)
box()

# =====================================================
# 3. FALSE NEGATIVE RATE
# =====================================================

plot(x_vals, rates_strm$FN_rate,
     type = "l", lwd = 3, col = "red",
     xlab = "", ylab = "",
     yaxt = "n", xaxt = "n",
     main = "False Negative Rate",
     ylim = range(c(rates_strm$FN_rate,
                    rates_samplecov$FN_rate,
                    rates_online$FN_rate))
)

lines(x_vals, rates_samplecov$FN_rate,
      lty = "dotted", col = "darkgreen", lwd = 3)

lines(x_vals, rates_online$FN_rate,
      lty = "dashed", col = "blue", lwd = 3)

axis(2, las = 1, cex.axis = 1.2)
axis(1, at = seq(1000, max(x_vals), by = 1000),
     las = 1, cex.axis = 1.2)
box()

# =====================================================
# 4. FALSE POSITIVE RATE
# =====================================================

plot(x_vals, rates_strm$FP_rate,
     type = "l", lwd = 3, col = "red",
     xlab = "", ylab = "",
     yaxt = "n", xaxt = "n",
     main = "False Positive Rate",
     ylim = range(c(rates_strm$FP_rate,
                    rates_samplecov$FP_rate,
                    rates_online$FP_rate))
)

lines(x_vals, rates_samplecov$FP_rate,
      lty = "dotted", col = "darkgreen", lwd = 3)

lines(x_vals, rates_online$FP_rate,
      lty = "dashed", col = "blue", lwd = 3)

axis(2, las = 1, cex.axis = 1.2)
axis(1, at = seq(1000, max(x_vals), by = 1000),
     las = 1, cex.axis = 1.2)
box()

dev.off()
