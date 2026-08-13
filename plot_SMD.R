########################Boxplot Sigma errors##############################


methodes = c("SampleNaiveQuantonlinecorr","OnlineUsQuantonlinecorr","StreamingUsonlineQuantcorr","Offline-withonlinequantile","OGK","MCD")
methodes_add  = c("SampleNaivewithoutonlinequantilecorr","OnlineUswithoutQuantonlinecorr","StreamingUswithoutQuantonlinecorr","OfflineUswithoutQuantcorr","OracleRD","OracleQC","SampleRaw","OnlRaw","StrmRaw","OfflRaw","OGKRD","OGKQC","MCDRD","MCDQC")
methode_oracle = c("Oracle")

erreursSigmaSMD <- array(0, dim = c(length(methodes), nbmachines))

setwd(crit_SMD)

for (j in seq_len(nbmachines)) {
  
  for (s in seq_along(methodes)) {
    
    critFile <- paste0("Crit-", methodes[s], "-machine-", j, ".RData")
    
    load(critFile)
    
    erreursSigmaSMD[s, j] <- crit$erreurFrob
  }
}

setwd(fig_SMD)

pdf("boxploterreursSigma_SMD.pdf", width = 18, height = 6)

# Transformation logarithmique
erreurs_log <- log1p(erreursSigmaSMD)

# Données à tracer
data_plot <- list(
  Streaming   = erreurs_log[3, ],
  MCD         = erreurs_log[6, ],
  OGK         = erreurs_log[5, ],
  Offline     = erreurs_log[4, ],
  Sample_cov  = erreurs_log[1, ],
  Online      = erreurs_log[2, ]
)

# Couleurs
cols <- c(
  Streaming  = "red",
  MCD        = "black",
  OGK        = "brown",
  Offline    = "orange",
  Sample_cov = "darkgreen",
  Online     = "blue"
)

# Graduation de l'axe Y
ticks_raw <- c(0, 1e-2, 1e-1)
ticks <- log1p(ticks_raw)

boxplot(
  data_plot,
  col = cols[names(data_plot)],
  names = c(
    "Streaming",
    "MCD",
    "OGK",
    "Offline",
    "Sample covariance",
    "Online"
  ),
  yaxt = "n",
  ylab = "",
  main = "",
  ylim = c(0, max(c(unlist(erreurs_log), ticks)))
)

axis(
  side = 2,
  at = ticks,
  labels = c("0", expression(10^{-2}), expression(10^{-1})),
  las = 1,
  cex.axis = 0.9
)
dev.off()

####################AUC ARI#######################

all_methodes = c(methodes,methodes_add,methode_oracle)

############################# STOCKAGE #################################
ariPlot <- matrix(0, nrow = nbmachines, ncol = length(all_methodes))
aucPlot <- matrix(0, nrow = nbmachines, ncol = length(methodes))
accPlot = matrix(0, nrow = nbmachines, ncol = length(all_methodes))
temps   <- matrix(0, nrow = nbmachines, ncol = length(methodes))

for (j in seq_len(nbmachines)) {
  
  setwd(crit_SMD)
  
  ## Critères des 6 méthodes principales
  for (s in seq_along(methodes)) {
    
    load(paste0("Crit-", methodes[s], "-machine-", j, ".RData"))
    
    aucPlot[j, s] <- crit$AUC
  }
  
  ## Temps
  setwd(res_SMD)
  
  for (s in seq_along(methodes)) {
    
    load(paste0("Fit-", methodes[s], "-machine-", j, ".RData"))
    
    temps[j, s] <- resultats$temps[3]
  }
  
  ## ARI et Accuracy pour toutes les méthodes
  setwd(crit_SMD)
  
  for (s in seq_along(all_methodes)) {
    
    load(paste0("Crit-", all_methodes[s], "-machine-", j, ".RData"))
    
    ariPlot[j, s] <- crit$ARI
    accPlot[j, s] <- 1 - crit$prop_hors_diag
  }
}
############################# LABELS #################################
setwd(fig_SMD)
names_plot_all <- c(
  "Sple QC",
  "Onl QC",
  "Strm QC",
  "Offl QC",
  "OGK QC",
  "MCD QC",
  "Oracle QC",
  
  "Sple RD",
  "Onl RD",
  "Strm RD",
  "Offl RD",
  "OGK RD",
  "MCD RD",
  "Oracle RD",
  
  "Sple Raw",
  "Onl Raw",
  "Strm Raw",
  "Offl Raw",
  "OGK Raw",
  "MCD Raw",
  "Oracle Raw"
)

############################# COULEURS #################################

cols_all <- c(
  # QC
  "darkgreen",
  "blue",
  "red",
  "orange",
  "brown",
  "black",
  "purple4",
  
  # RD
  "darkgreen",
  "blue",
  "red",
  "orange",
  "brown",
  "black",
  "purple4",
  
  # Raw
  "darkgreen",
  "blue",
  "red",
  "orange",
  "brown",
  "black",
  "purple4"
)

############################# ARI #################################

pdf("boxplotARI_SMD.pdf", width = 20, height = 7)

boxplot(
  ariPlot,
  names = names_plot_all,
  las = 2,
  col = cols_all,
  main = "ARI",
  ylab = "ARI",
  cex.axis = 1.2
)

dev.off()

############################# ACCURACY #################################

pdf("boxplotACC_SMD.pdf", width = 20, height = 7)

boxplot(
  accPlot,
  names = names_plot_all,
  las = 2,
  col = cols_all,
  main = "Accuracy",
  ylab = "Accuracy",
  cex.axis = 1.2
)

dev.off()

######################AUC##################################


cols_auc <- c(
  "darkgreen",
  "blue",
  "red",
  "orange",
  "brown",
  "black"
)

############################# AUC #################################
# Méthodes Frobenius / AUC / Temps
names_plot_auc <- c(
  "Sple",
  "Onl",
  "Strm",
  "Offl",
  "OGK",
  "MCD"
)

cols_auc <- c(
  "darkgreen",
  "blue",
  "red",
  "orange",
  "brown",
  "black"
)
pdf("boxplotAUC_SMD.pdf", width = 12, height = 6)

boxplot(
  aucPlot,
  names = names_plot_auc,
  las = 2,
  col = cols_auc,
  main = "AUC",
  ylab = "AUC",
  cex.axis = 1.2
)

dev.off()

######################Temps##################


temps_log <- log1p(temps)

tick_vals <- c(0, 1e-3, 1e-2, 1e-1, 1, 10)
tick_pos  <- log1p(tick_vals)
tick_lab  <- c("0","1e-3","1e-2","1e-1","1","10")

pdf("boxplotTemps_SMD.pdf", width = 12, height = 6)

boxplot(
  temps_log,
  names = names_plot_auc,
  las = 2,
  col = cols_auc,
  main = "Computation times",
  ylab = "",
  yaxt = "n",
  cex.axis = 1.2
)

axis(
  2,
  at = tick_pos,
  labels = tick_lab,
  las = 1,
  cex.axis = 1.2
)

dev.off()
############################Trajectoires######################
#===============================
# Méthodes online
#===============================

methodes_online_quantile <- c(
  "SampleNaiveQuantonlinecorr",
  "OnlineUsQuantonlinecorr",
  "StreamingUsonlineQuantcorr",
  "OracleQC"
)

methodes_online_rescale <- c(
  "SampleNaivewithoutonlinequantilecorr",
  "OnlineUswithoutQuantonlinecorr",
  "StreamingUswithoutQuantonlinecorr",
  "OracleRD"
)

methodes_online_raw <- c(
  "SampleRaw",
  "OnlRaw",
  "StrmRaw",
  "Oracle"
)

methodes_online <- c(
  methodes_online_quantile,
  methodes_online_rescale,
  methodes_online_raw
)


#===============================
# Symboles associés aux familles
#===============================

pch_methodes <- rep(NA, length(methodes_online))
names(pch_methodes) <- methodes_online

# Correction quantile : étoile
pch_methodes[methodes_online_quantile] <- 8

# Rescale distance : carré
pch_methodes[methodes_online_rescale] <- 15

# Raw : triangle
pch_methodes[methodes_online_raw] <- 17



#### Calcul des trajectoires ####

for(j in 1:nbmachines){
  
  setwd(smd_data_dir)
  
  data_smd_mach <- paste0(
    "data_machine-", j, ".RData"
  )
  
  load(data_smd_mach)
  
  outlmach <- matrix(
    0,
    nrow = nrow(Z),
    ncol = length(methodes_online)
  )
  
  colnames(outlmach) <- methodes_online
  
  setwd(res_SMD)
  
  distoracle <- rep(0, nrow(Z))
  
  
  #=========================================================
  # Chargement des résultats
  #=========================================================
  
  for(s in seq_along(methodes_online)){
    
    methode <- methodes_online[s]
    
    fitFile <- paste0(
      "Fit-",
      methode,
      "-machine-",
      j,
      ".RData"
    )
    
    load(fitFile)
    
    outlmach[, s] <- resultats$outliers_labels
    
    # Récupération des distances de l'Oracle
    if(methode == "Oracle"){
      distoracle <- resultats$distances
    }
  }
  
  
  #=========================================================
  # Calcul des taux
  #=========================================================
  
  
  #-------------------------
  # QC
  #-------------------------
  
  rates_samplecov_quantcorr <-
    compute_rates(
      outlmach[, "SampleNaiveQuantonlinecorr"],
      labels
    )
  
  rates_online_with_quantcorr <-
    compute_rates(
      outlmach[, "OnlineUsQuantonlinecorr"],
      labels
    )
  
  rates_Strm_with_quantcorr <-
    compute_rates(
      outlmach[, "StreamingUsonlineQuantcorr"],
      labels
    )
  
  rates_oracle_qc <-
    compute_rates(
      outlmach[, "OracleQC"],
      labels
    )
  
  
  #-------------------------
  # RD
  #-------------------------
  
  rates_samplecov_without_quantcorr <-
    compute_rates(
      outlmach[, "SampleNaivewithoutonlinequantilecorr"],
      labels
    )
  
  rates_online_without_quantcorr <-
    compute_rates(
      outlmach[, "OnlineUswithoutQuantonlinecorr"],
      labels
    )
  
  rates_Strm_without_quantcorr <-
    compute_rates(
      outlmach[, "StreamingUswithoutQuantonlinecorr"],
      labels
    )
  
  rates_oracle_rd <-
    compute_rates(
      outlmach[, "OracleRD"],
      labels
    )
  
  
  #-------------------------
  # Raw
  #-------------------------
  
  rates_samplecov_raw <-
    compute_rates(
      outlmach[, "SampleRaw"],
      labels
    )
  
  rates_online_raw <-
    compute_rates(
      outlmach[, "OnlRaw"],
      labels
    )
  
  rates_Strm_raw <-
    compute_rates(
      outlmach[, "StrmRaw"],
      labels
    )
  
  rates_oracle_raw <-
    compute_rates(
      outlmach[, "Oracle"],
      labels
    )
  
  
  
  ##################################
  # Trajectoires
  ##################################
  
  setwd(fig_SMD)
  
  nom_fichier <- paste0(
    "trajectories_SMD_mach-",
    j,
    ".pdf"
  )
  
  pdf(
    nom_fichier,
    width = 14,
    height = 10
  )
  
  
  par(
    mfrow = c(2,2),
    mar = c(4,4,2,1)
  )
  
  
  x_vals <- 1:length(
    rates_Strm_with_quantcorr$FN_rate
  )
  
  # Position des symboles
  idx_symbols <- seq(
    1,
    length(x_vals),
    by = 500
  )
  
  
  
  ##################################
  # 1. BOXPLOT DISTANCES
  ##################################
  
  Z_clean <- Z[labels == 0, , drop = FALSE]
  Z_outliers <- Z[labels == 1, , drop = FALSE]
  
  
  distinliers <- rep(
    0,
    nrow(Z_clean)
  )
  
  distoutliers <- rep(
    0,
    nrow(Z_outliers)
  )
  
  
  invSigmaTrueCov <- solve(
    cov(Z_clean)
  )
  
  mu <- colMeans(Z_clean)
  
  
  for(m in 1:nrow(Z_clean)){
    
    diff <- Z_clean[m,] - mu
    
    distinliers[m] <-
      t(diff) %*%
      invSigmaTrueCov %*%
      diff
  }
  
  
  for(m in 1:nrow(Z_outliers)){
    
    diff <- Z_outliers[m,] - mu
    
    distoutliers[m] <-
      t(diff) %*%
      invSigmaTrueCov %*%
      diff
  }
  
  
  boxplot(
    distinliers,
    distoutliers,
    col = c("lightblue","red"),
    names = c("0","1"),
    ylim = c(0,400)
  )
  
  
  
  ##################################
  # 2. FALSE NEGATIVE RATE
  ##################################
  
  plot(
    x_vals,
    rates_Strm_with_quantcorr$FN_rate * 100,
    type = "l",
    lwd = 3,
    col = "red",
    ylim = c(0,100),
    xlab = "",
    ylab = "",
    main = "False Negative Rate",
    xaxt = "n",
    yaxt = "n"
  )
  
  
  #=========================================================
  # QC
  #=========================================================
  
  
  # Sample QC
  lines(
    x_vals,
    rates_samplecov_quantcorr$FN_rate * 100,
    lty = "dotted",
    col = "darkgreen",
    lwd = 3
  )
  
  points(
    x_vals[idx_symbols],
    rates_samplecov_quantcorr$FN_rate[idx_symbols] * 100,
    pch = pch_methodes[
      "SampleNaiveQuantonlinecorr"
    ],
    col = "darkgreen",
    cex = 1.4
  )
  
  
  # Online QC
  lines(
    x_vals,
    rates_online_with_quantcorr$FN_rate * 100,
    lty = "dashed",
    col = "blue",
    lwd = 3
  )
  
  points(
    x_vals[idx_symbols],
    rates_online_with_quantcorr$FN_rate[idx_symbols] * 100,
    pch = pch_methodes[
      "OnlineUsQuantonlinecorr"
    ],
    col = "blue",
    cex = 1.4
  )
  
  
  # Streaming QC
  lines(
    x_vals,
    rates_Strm_with_quantcorr$FN_rate * 100,
    lty = "solid",
    col = "red",
    lwd = 3
  )
  
  points(
    x_vals[idx_symbols],
    rates_Strm_with_quantcorr$FN_rate[idx_symbols] * 100,
    pch = pch_methodes[
      "StreamingUsonlineQuantcorr"
    ],
    col = "red",
    cex = 1.4
  )
  
  
  # Oracle QC
  lines(
    x_vals,
    rates_oracle_qc$FN_rate * 100,
    lty = "dashed",
    col = "purple",
    lwd = 3
  )
  
  points(
    x_vals[idx_symbols],
    rates_oracle_qc$FN_rate[idx_symbols] * 100,
    pch = pch_methodes["OracleQC"],
    col = "purple",
    cex = 1.4
  )
  
  
  
  #=========================================================
  # RD
  #=========================================================
  
  
  # Sample RD
  lines(
    x_vals,
    rates_samplecov_without_quantcorr$FN_rate * 100,
    lty = "dotdash",
    col = "darkgreen",
    lwd = 3
  )
  
  points(
    x_vals[idx_symbols],
    rates_samplecov_without_quantcorr$FN_rate[idx_symbols] * 100,
    pch = pch_methodes[
      "SampleNaivewithoutonlinequantilecorr"
    ],
    col = "darkgreen",
    cex = 1.4
  )
  
  
  # Online RD
  lines(
    x_vals,
    rates_online_without_quantcorr$FN_rate * 100,
    lty = "twodash",
    col = "blue",
    lwd = 3
  )
  
  points(
    x_vals[idx_symbols],
    rates_online_without_quantcorr$FN_rate[idx_symbols] * 100,
    pch = pch_methodes[
      "OnlineUswithoutQuantonlinecorr"
    ],
    col = "blue",
    cex = 1.4
  )
  
  
  # Streaming RD
  lines(
    x_vals,
    rates_Strm_without_quantcorr$FN_rate * 100,
    lty = "longdash",
    col = "red",
    lwd = 3
  )
  
  points(
    x_vals[idx_symbols],
    rates_Strm_without_quantcorr$FN_rate[idx_symbols] * 100,
    pch = pch_methodes[
      "StreamingUswithoutQuantonlinecorr"
    ],
    col = "red",
    cex = 1.4
  )
  
  
  # Oracle RD
  lines(
    x_vals,
    rates_oracle_rd$FN_rate * 100,
    lty = "longdash",
    col = "purple",
    lwd = 3
  )
  
  points(
    x_vals[idx_symbols],
    rates_oracle_rd$FN_rate[idx_symbols] * 100,
    pch = pch_methodes["OracleRD"],
    col = "purple",
    cex = 1.4
  )
  
  
  
  #=========================================================
  # Raw
  #=========================================================
  
  
  # Sample Raw
  lines(
    x_vals,
    rates_samplecov_raw$FN_rate * 100,
    lty = "solid",
    col = "darkgreen",
    lwd = 3
  )
  
  points(
    x_vals[idx_symbols],
    rates_samplecov_raw$FN_rate[idx_symbols] * 100,
    pch = pch_methodes["SampleRaw"],
    col = "darkgreen",
    cex = 1.4
  )
  
  
  # Online Raw
  lines(
    x_vals,
    rates_online_raw$FN_rate * 100,
    lty = "dashed",
    col = "blue",
    lwd = 3
  )
  
  points(
    x_vals[idx_symbols],
    rates_online_raw$FN_rate[idx_symbols] * 100,
    pch = pch_methodes["OnlRaw"],
    col = "blue",
    cex = 1.4
  )
  
  
  # Streaming Raw
  lines(
    x_vals,
    rates_Strm_raw$FN_rate * 100,
    lty = "solid",
    col = "red",
    lwd = 3
  )
  
  points(
    x_vals[idx_symbols],
    rates_Strm_raw$FN_rate[idx_symbols] * 100,
    pch = pch_methodes["StrmRaw"],
    col = "red",
    cex = 1.4
  )
  
  
  # Oracle Raw
  lines(
    x_vals,
    rates_oracle_raw$FN_rate * 100,
    lty = "solid",
    col = "purple",
    lwd = 3
  )
  
  points(
    x_vals[idx_symbols],
    rates_oracle_raw$FN_rate[idx_symbols] * 100,
    pch = pch_methodes["Oracle"],
    col = "purple",
    cex = 1.4
  )
  
  
  
  axis(
    2,
    las = 1,
    cex.axis = 1.8
  )
  
  axis(
    1,
    at = seq(
      1000,
      max(x_vals),
      by = 1000
    ),
    las = 1,
    cex.axis = 1.8
  )
  
  box()
  
  
  
  ##################################
  # 3. FALSE POSITIVE RATE
  ##################################
  
  plot(
    x_vals,
    rates_Strm_with_quantcorr$FP_rate * 100,
    type = "l",
    lwd = 3,
    col = "red",
    ylim = c(0,20),
    xlab = "",
    ylab = "",
    main = "False Positive Rate",
    xaxt = "n",
    yaxt = "n"
  )
  
  
  #=========================================================
  # QC
  #=========================================================
  
  
  # Sample QC
  lines(
    x_vals,
    rates_samplecov_quantcorr$FP_rate * 100,
    lty = "dotted",
    col = "darkgreen",
    lwd = 3
  )
  
  points(
    x_vals[idx_symbols],
    rates_samplecov_quantcorr$FP_rate[idx_symbols] * 100,
    pch = pch_methodes[
      "SampleNaiveQuantonlinecorr"
    ],
    col = "darkgreen",
    cex = 1.4
  )
  
  
  # Online QC
  lines(
    x_vals,
    rates_online_with_quantcorr$FP_rate * 100,
    lty = "dashed",
    col = "blue",
    lwd = 3
  )
  
  points(
    x_vals[idx_symbols],
    rates_online_with_quantcorr$FP_rate[idx_symbols] * 100,
    pch = pch_methodes[
      "OnlineUsQuantonlinecorr"
    ],
    col = "blue",
    cex = 1.4
  )
  
  
  # Streaming QC
  lines(
    x_vals,
    rates_Strm_with_quantcorr$FP_rate * 100,
    lty = "solid",
    col = "red",
    lwd = 3
  )
  
  points(
    x_vals[idx_symbols],
    rates_Strm_with_quantcorr$FP_rate[idx_symbols] * 100,
    pch = pch_methodes[
      "StreamingUsonlineQuantcorr"
    ],
    col = "red",
    cex = 1.4
  )
  
  
  # Oracle QC
  lines(
    x_vals,
    rates_oracle_qc$FP_rate * 100,
    lty = "dashed",
    col = "purple",
    lwd = 3
  )
  
  points(
    x_vals[idx_symbols],
    rates_oracle_qc$FP_rate[idx_symbols] * 100,
    pch = pch_methodes["OracleQC"],
    col = "purple",
    cex = 1.4
  )
  
  
  
  #=========================================================
  # RD
  #=========================================================
  
  
  # Sample RD
  lines(
    x_vals,
    rates_samplecov_without_quantcorr$FP_rate * 100,
    lty = "dotdash",
    col = "darkgreen",
    lwd = 3
  )
  
  points(
    x_vals[idx_symbols],
    rates_samplecov_without_quantcorr$FP_rate[idx_symbols] * 100,
    pch = pch_methodes[
      "SampleNaivewithoutonlinequantilecorr"
    ],
    col = "darkgreen",
    cex = 1.4
  )
  
  
  # Online RD
  lines(
    x_vals,
    rates_online_without_quantcorr$FP_rate * 100,
    lty = "twodash",
    col = "blue",
    lwd = 3
  )
  
  points(
    x_vals[idx_symbols],
    rates_online_without_quantcorr$FP_rate[idx_symbols] * 100,
    pch = pch_methodes[
      "OnlineUswithoutQuantonlinecorr"
    ],
    col = "blue",
    cex = 1.4
  )
  
  
  # Streaming RD
  lines(
    x_vals,
    rates_Strm_without_quantcorr$FP_rate * 100,
    lty = "longdash",
    col = "red",
    lwd = 3
  )
  
  points(
    x_vals[idx_symbols],
    rates_Strm_without_quantcorr$FP_rate[idx_symbols] * 100,
    pch = pch_methodes[
      "StreamingUswithoutQuantonlinecorr"
    ],
    col = "red",
    cex = 1.4
  )
  
  
  # Oracle RD
  lines(
    x_vals,
    rates_oracle_rd$FP_rate * 100,
    lty = "longdash",
    col = "purple",
    lwd = 3
  )
  
  points(
    x_vals[idx_symbols],
    rates_oracle_rd$FP_rate[idx_symbols] * 100,
    pch = pch_methodes["OracleRD"],
    col = "purple",
    cex = 1.4
  )
  
  
  
  #=========================================================
  # Raw
  #=========================================================
  
  
  # Sample Raw
  lines(
    x_vals,
    rates_samplecov_raw$FP_rate * 100,
    lty = "solid",
    col = "darkgreen",
    lwd = 3
  )
  
  points(
    x_vals[idx_symbols],
    rates_samplecov_raw$FP_rate[idx_symbols] * 100,
    pch = pch_methodes["SampleRaw"],
    col = "darkgreen",
    cex = 1.4
  )
  
  
  # Online Raw
  lines(
    x_vals,
    rates_online_raw$FP_rate * 100,
    lty = "dashed",
    col = "blue",
    lwd = 3
  )
  
  points(
    x_vals[idx_symbols],
    rates_online_raw$FP_rate[idx_symbols] * 100,
    pch = pch_methodes["OnlRaw"],
    col = "blue",
    cex = 1.4
  )
  
  
  # Streaming Raw
  lines(
    x_vals,
    rates_Strm_raw$FP_rate * 100,
    lty = "solid",
    col = "red",
    lwd = 3
  )
  
  points(
    x_vals[idx_symbols],
    rates_Strm_raw$FP_rate[idx_symbols] * 100,
    pch = pch_methodes["StrmRaw"],
    col = "red",
    cex = 1.4
  )
  
  
  # Oracle Raw
  lines(
    x_vals,
    rates_oracle_raw$FP_rate * 100,
    lty = "solid",
    col = "purple",
    lwd = 3
  )
  
  points(
    x_vals[idx_symbols],
    rates_oracle_raw$FP_rate[idx_symbols] * 100,
    pch = pch_methodes["Oracle"],
    col = "purple",
    cex = 1.4
  )
  
  
  
  axis(
    2,
    las = 1,
    cex.axis = 1.8
  )
  
  axis(
    1,
    at = seq(
      1000,
      max(x_vals),
      by = 1000
    ),
    las = 1,
    cex.axis = 1.8
  )
  
  box()
  
  
  dev.off()
}
