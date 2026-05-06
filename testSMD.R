url <- "https://github.com/NetManAIOps/OmniAnomaly/archive/refs/heads/master.zip"
dest <- "smd.zip"

download.file(url, destfile = dest, mode = "wb")

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

Z <- smd_data[[2]]$Z
labels <- smd_data[[2]]$labels
Z <- as.matrix(Z)

#Enlèvement des colonnes qui posent problème
Z <- Z[, !(colnames(Z) %in% c("V5","V8","V9","V10","V11", "V13","V17","V18","V23","V25","V27", "V29", "V30","V32","V33","V35","V36","V37", "V38")), drop = FALSE]

res = onlineRobustVariance(Z,batch = 4,computeOutliers = TRUE)
resmcd = covMcd(Z)
distmcd = rep(0,nrow(Z))

resogk = covOGK(Z, sigmamu = scaleTau2)
distogk = rep(0,nrow(Z))
distoffl = rep(0,nrow(Z))


for(i in 1:nrow(Z))
{
  distmcd[i] = t(Z[i,] - resmcd$center)%*%resmcd$cov%*%(Z[i,] - resmcd$center)
  distogk[i] = t(Z[i,] - resogk$center)%*%resogk$cov%*%(Z[i,] - resogk$center)
  
  distoffl[i] = t(Z[i,] - resogk$center)%*%res$variance%*%(Z[i,] - resogk$center)
  }



library(pROC)

roc_obj <- roc(res$outlier_labels,distoffl)
auc_value <- auc(roc_obj)

cat("AUC :", auc_value, "\n")


# =========================================================
# 14. COURBE ROC
# =========================================================

plot(roc_obj,
     xlim = c(0,1),
     ylim = c(0,1),
     )

