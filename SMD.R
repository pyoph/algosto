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



#Enlèvement des colonnes qui posent problème d'inversibilité
Z <- Z[, !(colnames(Z) %in% c("V5","V8","V9","V10","V11", "V13","V17","V18","V23","V25","V27", "V29", "V30","V32","V33","V35","V36","V37", "V38")), drop = FALSE]



resStrm = onlineRobustVariance(Z,batch = ncol(Z),computeOutliers = TRUE)
resOnline = onlineRobustVariance(Z,batch = 1,computeOutliers = TRUE)
resOffline = offlineRobustVariance(Z,computeOutliers = TRUE)
resmcd = covMcd(Z)

resogk = covOGK(Z, sigmamu = scaleTau2)
distogk = rep(0,nrow(Z))
distoffl = rep(0,nrow(Z))
distmcd = rep(0,nrow(Z))
invSigmaOGK = solve(resogk$cov)
invSigmaMCD = solve(resmcd$cov)

#Calcule destinations de Mahalanobis

for(i in 1:nrow(Z))
{
  distmcd[i] = t(Z[i,] - resmcd$center)%*%invSigmaMCD%*%(Z[i,] - resmcd$center)
  distogk[i] = t(Z[i,] - resogk$center)%*%invSigmaOGK%*%(Z[i,] - resogk$center)
 # distoffl[i] = t(Z[i,] - t(resOffline$median))%*%invSigmaOffl%*%(Z[i,] - t(resOffline$median))
  }


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


roc_strm    <- compute_roc(resStrm$distances, labels)
roc_online  <- compute_roc(resOnline$distances, labels)
roc_offline <- compute_roc(resOffline$distances, labels)
roc_mcd     <- compute_roc(distmcd, labels)
roc_ogk     <- compute_roc(distogk, labels)
dev.off()
plot(roc_strm$FPR, roc_strm$TPR,
     type = "l", col = "red", lwd = 2,
     xlim = c(0,1), ylim = c(0,1),
     xlab = "FPR", ylab = "TPR",
     main = "ROC Comparison")

lines(roc_online$FPR, roc_online$TPR, col = "blue", lwd = 2)
lines(roc_offline$FPR, roc_offline$TPR, col = "orange", lwd = 2)
lines(roc_mcd$FPR, roc_mcd$TPR, col = "black", lwd = 2)
lines(roc_ogk$FPR, roc_ogk$TPR, col = "brown", lwd = 2)

abline(0,1,lty=2)


save(roc_strm,roc_online,roc_offline,roc_mcd,roc_ogk,file = "roc.RData")


################Choix du meilleur seuil##############################################
# =========================
# GRID DE SEUILS
# =========================

thresholds <- seq(
  min(scores),
  max(scores),
  length.out = 200
)

# =========================
# FONCTION F1
# =========================

f1_score <- function(pred, true) {
  
  TP <- sum(pred == 1 & true == 1)
  FP <- sum(pred == 1 & true == 0)
  FN <- sum(pred == 0 & true == 1)
  
  precision <- ifelse(TP + FP == 0, 0, TP / (TP + FP))
  recall    <- ifelse(TP + FN == 0, 0, TP / (TP + FN))
  
  if (precision + recall == 0) return(0)
  
  2 * precision * recall / (precision + recall)
}

# =========================
# OPTIMISATION SEUIL
# =========================

f1_values <- numeric(length(thresholds))

for (i in seq_along(thresholds)) {
  
  tau <- thresholds[i]
  
  pred <- as.integer(scores > tau)
  
  f1_values[i] <- f1_score(pred, labels)
}

# =========================
# MEILLEUR SEUIL
# =========================

best_idx <- which.max(f1_values)

best_threshold <- thresholds[best_idx]

best_f1 <- f1_values[best_idx]

best_threshold
best_f1
