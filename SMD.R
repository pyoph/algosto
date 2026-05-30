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
         nDataInit = 1e3,
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
      # 
       onlineRobustVariance(
         Z[, c(keep, j), drop = FALSE],
         computeOutliers = TRUE,
         cutoff = cut,
         cutinit = 0.6,
         nDataInit = 1e3,
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

setwd("~")

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
temps = array(0,dim=c(6,nbmachines))
keep_columns_ok = list()

for(j in 1:nbmachines){
  
  cat("\n====================\n")
  cat("Machine", j, "\n")
  cat("====================\n")
  
  
  Z <- smd_data[[j]]$Z
  labels <- smd_data[[j]]$labels
  Z <- as.matrix(Z)
  
  #Z = Z[1:2.3e4,]
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
  
  
  Z = Z[, keep_columns_ok[[j]], drop = FALSE]
  
  
  
  #Z = Z[,keep_columns_ok[[j]],drop = FALSE]
  col_before <- colnames(Z)
  
  Z = remove_high_corr(Z,threshold = 0.99)
  col_after <- colnames(Z)
  keep_idx <- match(col_after, col_before)
  
  # mise à jour des colonnes globales gardées
  keep_columns_ok[[j]] <- keep_columns_ok[[j]][keep_idx]
  #Z = Z[,keep_columns_ok[[j]],drop = FALSE]
  Z_clean = Z[labels == 0,]
  
  ##############Calcul taille batch
  
  #table(res0$outlier_labels,labels)
  possible_div <- which(nrow(Z) %% 1:nrow(Z) == 0)
  
  # ne garder que ceux <= ncol(Z)
  possible_div <- possible_div[possible_div <= ncol(Z)]
  
  # prendre le plus grand
  batchStrm <- max(possible_div)

  ###########Extraction moyenne et covariance sans outliers
  meanTrueCov = colMeans(Z_clean)
  
  Sigma_trueCov = cov(Z_clean)
  
  Sigma_trueMCD = covMcd(Z_clean)$cov

  Sigma_trueOGK = covOGK(Z_clean,sigmamu = scaleTau2)$cov
  
  Sigma_trueOffl = offlineRobustVariance(Z_clean)$variance
  
  Sigma_trueOnl = onlineRobustVariance(Z_clean,computeOutliers = TRUE,cutoff=0.95,cutinit=0.6,nDataInit = 2e3,c_m=1,batch = 1)$variance
  
  Sigma_trueStrm = onlineRobustVariance(Z_clean,computeOutliers = TRUE,cutoff=0.95,cutinit=0.6,nDataInit = 2e3,c_m=1,batch = batchStrm)$variance
  
  temps_covonline = system.time(
  {
    resSamplecov= SampleCovOnline(Z,quantcutoff = TRUE,nDataInit = 2e3,cutoffquant =  .95,c_m = .1)
  })
  temps[1,j]= temps_covonline[3]
  
  erreursSigmaSMD[1,j] = norm(resSamplecov$Sigma - Sigma_trueCov,"F")
  
  #table(resSamplecov$outliers_labels,labels)
  
  # resStrm = onlineRobustVariance(Z,batch = ncol(Z),computeOutliers = TRUE)
  # resOnline = onlineRobustVariance(Z,batch = 1,computeOutliers = TRUE)
  # resOffline = offlineRobustVariance(Z,computeOutliers = TRUE)
  distmcd = rep(0,nrow(Z))
  
  temps_mcd = system.time(
  
    {
      
    resmcd = covMcd(Z)
    invSigmaMCD = inv_safe(resmcd$cov)
    
    for(i in 1:nrow(Z))
    {
      distmcd[i] = t(Z[i,] - (resmcd$center))%*%invSigmaMCD%*%(Z[i,] - (resmcd$center))
      #distogk[i] = t(Z[i,] - resogk$center)%*%invSigmaOGK%*%(Z[i,] - resogk$center)
     # distrue[i] = t(Z[i,] - meanTrueCov)%*%invSigmaTrueCov%*%(Z[i,] - meanTrueCov)
      
    }
    }
  )
  temps[2,j]= temps_mcd[3]
  
  
  erreursSigmaSMD[2,j] = norm(resmcd$cov - Sigma_trueMCD,"F")
  
  #resOfflineClean = offlineRobustVariance(Z_clean,computeOutliers = TRUE)
  temps_ogk = system.time(
  {resogk = covOGK(Z, sigmamu = scaleTau2)}
  )
  temps[3,j] = temps_ogk[3] 
  
  erreursSigmaSMD[3,j] = norm(resogk$cov - Sigma_trueOGK,"F")
  
  #distogk = rep(0,nrow(Z))
  #distoffl = rep(0,nrow(Z))
  #invSigmaOGK = inv_safe(resogk$cov)
  
  invSigmaTrueCov = inv_safe(Sigma_trueCov)
  #med = t(WeiszfeldMedian(Z))
  #quantemp = 0
  outliersLabelsSMD[[j]][, 1] <- resSamplecov$outliers_labels
  distrue = rep(0,nrow(Z))
  
  for(i in 1:nrow(Z))
  {
    distmcd[i] = t(Z[i,] - (resmcd$center))%*%invSigmaMCD%*%(Z[i,] - (resmcd$center))
    #distogk[i] = t(Z[i,] - resogk$center)%*%invSigmaOGK%*%(Z[i,] - resogk$center)
    distrue[i] = t(Z[i,] - meanTrueCov)%*%invSigmaTrueCov%*%(Z[i,] - meanTrueCov)
    
  }
  distrueClean = rep(0,nrow(Z_clean))
  for (s in 1:nrow(Z_clean))
  {
    distrueClean[s] = t(Z_clean[s,] - meanTrueCov)%*%invSigmaTrueCov%*%(Z_clean[s,] - meanTrueCov)
  }
  
  ############################Inliers and outliers###################
  setwd("~/figures")
  
  file = paste0("boxplotinliersoutliers-",j,".pdf")
  
  pdf(file, width = 18, height = 6) 
  
  
  dist_outliers <- distrue[labels == 1]
  dist_inliers  <- distrue[labels == 0]
  
  
  boxplot(
    list(
      Inliers = dist_inliers,
      Outliers = dist_outliers
    ),
    col = c("blue", "red"),
    border = "black",
    log = "y",   # échelle log10 sur l'axe Y
    yaxt = "n",   # on supprime l’axe par défaut
    ylab = "",
    main = ""
  )
  
  # ticks en puissances de 10
  y_ticks <- 10^seq(
    floor(log10(min(c(dist_inliers, dist_outliers)))),
    ceiling(log10(max(c(dist_inliers, dist_outliers)))),
    by = 1
  )
  
  axis(2,
       at = y_ticks,
       labels = y_ticks)
  
  dev.off()
  
  
  
  distogk = resogk$distances
  #############################################Quantile oracle###############################
  ###########################################################################################
  quantoracle = quantile(distrueClean,.95)
  
  #oracle = rep(0,nrow(Z))
  
  #quantemp = median(distoffl)
  #pred = rep(0,nrow(Z))
  
  #cutoff_true = calcule_cutoff(distrue,c_m = 2,cutinit=0.6,cutoff = 0.9,n = nrow(Z))
  
  
  cut_off_ogk = qchisq(.95,df = ncol(Z))
  cut_off_mcd = qchisq(.95,df = ncol(Z))
  #cut_off_ogk = calcule_cutoff(distogk,c_m = 2,cutinit=0.6,cutoff = 0.9,n = nrow(Z))
  #cut_off_mcd = calcule_cutoff(distmcd,c_m = 2,cutinit=0.6,cutoff = 0.9,n = nrow(Z))
  
  outliersLabelsSMD[[j]][, 7] = as.integer(distrue > quantoracle)
  outliersLabelsSMD[[j]][, 3] =  as.integer(distogk > cut_off_ogk)
  outliersLabelsSMD[[j]][, 2] =  as.integer(distmcd > cut_off_mcd)
  # for (i in 1:nrow(Z)) {
  #   
  #   if (distrue[i] > quantoracle) {
  #     outliersLabelsSMD[[j]][i, 7] <- 1
  #   }
  #   
  #   if (distogk[i] > cut_off_ogk) {
  #     outliersLabelsSMD[[j]][i, 3] <- 1
  #   }
  #   
  #   if (distmcd[i] > cut_off_mcd) {
  #     outliersLabelsSMD[[j]][i, 2] <- 1
  #   }
  #   
  # }
  #pred = distoffl > quantemp
  
  #table(pred,labels)
  
  #table(oracle,labels)
  cut=0.95
  
  temps_offline = system.time(
    {
    res0 <- offlineRobustVariance(Z,computeOutliers = TRUE,cutoff=cut,cutinit=0.6,c_m=2)
     # res0 <- offlineRobustVariance_old(Z,computeOutliers = TRUE)
    }
  )
  temps[4,j] = temps_offline[3]
  
  erreursSigmaSMD[4,j] = norm(res0$variance - Sigma_trueOffl,"F")
  outliersLabelsSMD[[j]][, 4] = res0$outlier_labels
  
  
  temps_strm = system.time(
  resStrm <- onlineRobustVariance(Z,computeOutliers = TRUE,cutoff=cut,cutinit=0.6,nDataInit = 2e3,c_m= 2,batch = batchStrm)
  #  resStrm <- onlineRobustVariance_old(Z,computeOutliers = TRUE,batch = batchStrm)
    )
  temps[5,j] = temps_strm[3]
  erreursSigmaSMD[5,j] = norm(resStrm$variance - Sigma_trueStrm,"F")
  outliersLabelsSMD[[j]][, 5] = resStrm$outlier_labels
  #table(resStrm$outlier_labels,labels)
  temps_online = system.time(
  resOnline <- onlineRobustVariance(Z,computeOutliers = TRUE,cutoff=cut,cutinit=0.5,nDataInit = 2e3,c_m=1,batch = 1)
    #resOnline <- onlineRobustVariance_old(Z,computeOutliers = TRUE,batch = 1)
    
  )
  temps[6,j] = temps_online[3]
  outliersLabelsSMD[[j]][, 6] = resOnline$outlier_labels
  #table(resOnline$outlier_labels,labels)
  erreursSigmaSMD[6,j] = norm(resOnline$variance - Sigma_trueOnl,"F")
  setwd("~/figures")
  
  file <- paste0("histogram_densityStrm-mach-",j,".pdf")
  pdf(file, width = 18, height = 6)
  
  x <- resStrm$distances
  x <- x[is.finite(x)]
  
  df <- ncol(Z)
  
  # quantile 95%
  q95 <- qchisq(0.95, df = df)
  
  # histogramme
  hist(x,
       probability = TRUE,
       breaks = 80,
       col = "lightblue",
       border = "white",
       main = "",
       xlab = "",
       ylab = "",
       xlim = c(0, 1e3),
       ylim = c(0,0.003))
  
  # densité χ²
  curve(dchisq(x, df = df),
        add = TRUE,
        col = "red",
        lwd = 3)
  
  # seuil théorique 95%
  abline(v = q95,
         col = "black",
         lwd = 3,
         lty = 2)
  
  dev.off()
  
  #table(resOnline$outlier_labels,labels)
  fileres = paste0("res-mach-",j,".RData")
  setwd("~/res")
  save(resSamplecov,resogk,resmcd,res0,resOnline,resStrm,file = fileres)
  
}



################Avec la correction#############

distances = res0$distances

#coef_cor <- qchisq(.5, df = ncol(Z)) / sapply(seq_along(distances), function(i) median(distances[1:i]))

#distances_cor <- distances * coef_cor

x = distances

outliersOffl_old = as.integer(x > qchisq(.95,df = 19))

table(outliersOffl_old,labels)


outliersOnl_old = as.integer(x > qchisq(.95,df = 19))

distances = resOnline$distances

#coef_cor <- qchisq(.5, df = ncol(Z)) / sapply(seq_along(distances), function(i) median(distances[1:i]))

#distances_cor <- distances * coef_cor

x = distances

outliersOnl_old = as.integer(x > qchisq(.95,df = 19))

distances = resStrm$distances

#coef_cor <- qchisq(.5, df = ncol(Z)) / sapply(seq_along(distances), function(i) median(distances[1:i]))

#distances_cor <- distances * coef_cor

x = distances

outliersStrm_old = as.integer(x > qchisq(.95,df = 19))

table(outliersStrm_old,labels)

distances = resSamplecov$distances

#coef_cor <- qchisq(.5, df = ncol(Z)) / sapply(seq_along(distances), function(i) median(distances[1:i]))

#distances_cor <- distances * coef_cor

x = distances

outliersSampleCov_old= as.integer(x > qchisq(.95,df = 19))

table(outliersSampleCov_old,labels)



file <- "histogram_densityStrmCorr.pdf"
pdf(file, width = 18, height = 6)

x <- distances_cor
x <- x[is.finite(x)]

df <- ncol(Z)

# quantile 95%
q95 <- qchisq(0.95, df = df)

# histogramme
hist(x,
     probability = TRUE,
     breaks = 80,
     col = "lightblue",
     border = "white",
     main = "",
     xlab = "",
     ylab = "",
     xlim = c(0, 1e3),
     ylim = c(0,0.003))

# densité χ²
curve(dchisq(x, df = df),
      add = TRUE,
      col = "red",
      lwd = 3)

# seuil théorique 95%
abline(v = q95,
       col = "black",
       lwd = 3,
       lty = 2)

dev.off()




file <- "histogram_densityOffline.pdf"
pdf(file, width = 18, height = 6)

x <- res0$distances
x <- x[is.finite(x)]

df <- ncol(Z)

# quantile 95%
q95 <- qchisq(0.95, df = df)

# histogramme
hist(x,
     probability = TRUE,
     breaks = 80,
     col = "lightblue",
     border = "white",
     main = "",
     xlab = "",
     ylab = "",
     xlim = c(0, 1e3),
     ylim = c(0,0.003))

# densité χ²
curve(dchisq(x, df = df),
      add = TRUE,
      col = "red",
      lwd = 3)

# seuil théorique 95%
abline(v = q95,
       col = "black",
       lwd = 3,
       lty = 2)

dev.off()




file <- "histogram_densityMCD.pdf"
pdf(file, width = 18, height = 6)

x <- distmcd
x <- x[is.finite(x)]

df <- ncol(Z)

# quantile 95%
q95 <- qchisq(0.95, df = df)

# histogramme
hist(x,
     probability = TRUE,
     breaks = 80,
     col = "lightblue",
     border = "white",
     main = "",
     xlab = "",
     ylab = "",
     xlim = c(0, 1e3),
     ylim = c(0,0.003))

# densité χ²
curve(dchisq(x, df = df),
      add = TRUE,
      col = "red",
      lwd = 3)

# seuil théorique 95%
abline(v = q95,
       col = "black",
       lwd = 3,
       lty = 2)

dev.off()


file <- "histogram_densityOGK.pdf"
pdf(file, width = 18, height = 6)

x <- resogk$distances
x <- x[is.finite(x)]

df <- ncol(Z)

# quantile 95%
q95 <- qchisq(0.95, df = df)

# histogramme
hist(x,
     probability = TRUE,
     breaks = 80,
     col = "lightblue",
     border = "white",
     main = "",
     xlab = "",
     ylab = "",
     xlim = c(0, 2*1e3),
     ylim = c(0,0.001))

# densité χ²
curve(dchisq(x, df = df),
      add = TRUE,
      col = "red",
      lwd = 3)

# seuil théorique 95%
abline(v = q95,
       col = "black",
       lwd = 3,
       lty = 2)

dev.off()




file <- paste0("boxplottemps_SMD", ".pdf")

pdf(file, width = 18, height = 6) 

boxplot(
  temps[1,],temps[2,],temps[3,],temps[4,],temps[5,],temps[6,],
  log = "y",                     # échelle logarithmique sur Y
  names = c("sample covariance online","mcd","ogk","offline","streaming","online"),
  main = "",
  ylab = "",
  col = rainbow(ncol(temps))
)
dev.off()


file <- paste0("boxploterreursSigma_SMD", ".pdf")



pdf(file, width = 18, height = 6) 

boxplot(
  erreursSigmaSMD[1,],erreursSigmaSMD[2,],erreursSigmaSMD[3,],erreursSigmaSMD[4,],erreursSigmaSMD[5,],erreursSigmaSMD[6,],
  log = "y",                     # échelle logarithmique sur Y
  names = c("sample covariance online","mcd","ogk","offline","streaming","online"),
  main = "",
  ylab = "",
  col = rainbow(ncol(erreursSigmaSMD))
)
dev.off()

#############################ROC#######################################

# R code: list of `Z` matrices + outlier visualization for each machine


# =====================================================
# 1. LOAD ALL MACHINES INTO A LIST OF Z MATRICES
# =====================================================

load_smd_machine <- function(test_path, label_path) {
  
  Z <- as.matrix(read.table(test_path, sep = ","))
  labels <- scan(label_path)
  
  list(
    Z = Z,
    labels = labels,
    name = basename(test_path)
  )
}

machines <- vector("list", length(test_files))

for (j in seq_along(test_files)) {
  machines[[j]] <- load_smd_machine(test_files[j], label_files[j])
}

Z_list <- lapply(machines, function(m) m$Z)
labels_list = lapply(machines, function(m) m$labels)
# Example:
# Z_list[[1]] = matrix for machine 1
# Z_list[[2]] = matrix for machine 2



robust_outlier_axes <- function(Z, outlier_idx) {
  
  Z <- as.matrix(Z)
  p <- ncol(Z)
  
  contrib <- numeric(p)
  
  for (j in 1:p) {
    
    medj <- median(Z[, j], na.rm = TRUE)
    madj <- mad(Z[, j], constant = 1, na.rm = TRUE)
    
    if (madj < 1e-12) {
      madj <- 1e-12
    }
    
    contrib[j] <- sum(abs(Z[outlier_idx, j] - medj) / madj)
  }
  
  order(contrib, decreasing = TRUE)[1:2]
}




plot_outliers_machine <- function(Z,
                                  labels,
                                  machine_name = "machine",
                                  save_plot = TRUE,
                                  output_dir = "plots_outliers") {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # indices of outliers
  out_idx <- which(labels == 1)
  
  # top 2 robust directions
  top2 <- robust_outlier_axes(Z, out_idx)

  file =paste0("outlier-mach-",j,".pdf")
  
  pdf(file, width = 18, height = 6) 
  
    
  # colors and symbols
  cols <- ifelse(labels == 1, "red", "blue")
  pch_vals <- ifelse(labels == 1, 4, 19)
  
  # positive values for log10
  Z_plot <- Z[, top2]
  Z_plot[Z_plot <= 0] <- 1e-6
  
  # log scale bounds
  xmin <- floor(log10(min(Z_plot[,1])))
  xmax <- ceiling(log10(max(Z_plot[,1])))
  
  ymin <- floor(log10(min(Z_plot[,2])))
  ymax <- ceiling(log10(max(Z_plot[,2])))
  
  # major ticks
  x_ticks <- xmin:xmax
  y_ticks <- ymin:ymax
  
  # minor log grid
  x_minor <- unlist(lapply(xmin:(xmax-1), function(k) log10((2:9) * 10^k)))
  y_minor <- unlist(lapply(ymin:(ymax-1), function(k) log10((2:9) * 10^k)))
  

  plot(log10(Z_plot[,1]),
       log10(Z_plot[,2]),
       col = cols,
       pch = pch_vals,
       cex = 1.1,
       xaxt = "n",
       yaxt = "n",
       #xlab = paste0("Variable ", top2[1]),
       #ylab = paste0("Variable ", top2[2]),
       main = "")
  
  # minor grid
  abline(v = x_minor, col = "grey90", lty = 3)
  abline(h = y_minor, col = "grey90", lty = 3)
  
  # major grid
  abline(v = x_ticks, col = "grey75", lty = 2)
  abline(h = y_ticks, col = "grey75", lty = 2)
  
  # axes
  axis(1,
       at = x_ticks,
       labels = parse(text = paste0("10^", x_ticks)),
       cex.axis = 1.8)
  
  axis(2,
       at = y_ticks,
       labels = parse(text = paste0("10^", y_ticks)),
       las = 1,
       cex.axis = 1.8)
  # 
  # legend("topright",
  #        legend = c("Normal", "Outlier"),
  #        col = c("blue", "red"),
  #        pch = c(19, 4),
  #        bty = "n")
  
}

setwd("~/figures")
# =====================================================
# 5. GENERATE PLOTS FOR ALL MACHINES
# =====================================================
for (j in seq_along(machines)) {
  
  cat("Machine:", j, "\n")
  
  Z <- machines[[j]]$Z
  labels <- machines[[j]]$labels
  name <- machines[[j]]$name
  
  tryCatch({
    
    plot_outliers_machine(Z,
                          labels,
                          machine_name = name,
                          save_plot = FALSE)
    
  }, error = function(e) {
    cat("Error for machine", j, ":", e$message, "\n")
  })
  dev.off()
}

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
  #rates_samplecov = compute_rates(outliersSampleCov_old, labels)
  rates_strm  <- compute_rates(outlmach[,5], labels)
  #rates_strm = compute_rates(outliersStrm_old,labels)
  rates_online <- compute_rates(outlmach[,6], labels)
  #rates_online <- compute_rates(outliersOnl_old, labels)
  rates_oracle <- compute_rates(outlmach[,7], labels)
  
  ##################################Trajectoires################################"
  setwd("~/figures")
  nom_fichier = paste0("trajectories_SMD_mach-",j,"-Ninit2e3cm1",".pdf")
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
       ylim = range(c(0,40))
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
  
  setwd("~/res")
  fileres= paste0("res-mach-",j,".RData")
  load(fileres)
  cutoff = qchisq(.95,df =ncol(resStrm$variance))
  outliersSampleCov_old = as.integer(resSamplecov$distances > cutoff )
  outliersStrm_old = as.integer(resStrm$distances > cutoff )
  outliersOnl_old = as.integer(resOnline$distances > cutoff )
  
  #rates_samplecov = compute_rates(outlmach[,1], labels)
  rates_samplecov = compute_rates(outliersSampleCov_old, labels)
  #rates_strm  <- compute_rates(outlmach[,5], labels)
  rates_strm = compute_rates(outliersStrm_old,labels)
  #rates_online <- compute_rates(outlmach[,6], labels)
  rates_online <- compute_rates(outliersOnl_old, labels)
  rates_oracle <- compute_rates(outlmach[,7], labels)
  
  ##################################Trajectoires################################"
  setwd("~/figures")
  nom_fichier = paste0("trajectories_SMD_mach-",j,"-old",".pdf")
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
       ylim = range(c(0,40))
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
