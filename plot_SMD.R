########################Boxplot Sigma errors##############################


methodes = c("SampleNaiveQuantonlinecorr","OnlineUsQuantonlinecorr","StreamingUsonlineQuantcorr","Offline-withonlinequantile","OGK","MCD")
methodes_add  = c("SampleNaivewithoutonlinequantilecorr","OnlineUswithoutQuantonlinecorr","StreamingUswithoutQuantonlinecorr","OfflineUswithoutQuantcorr","OracleRD","OracleQC","SampleRaw","OnlRaw","StrmRaw","OfflRaw","OGKRD","OGKQC","MCDRD","MCDQC")
methode_oracle = c("Oracle")

erreursSigmaSMD <- array(0, dim = c(length(methodes), nbmachines))

setwd(crit_SMD)

for (j in machines_ok_final_boucle) {
  
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
nbmachines = 56
############################# STOCKAGE #################################
ariPlot <- matrix(0, nrow = nbmachines, ncol = length(all_methodes))
aucPlot <- matrix(0, nrow = nbmachines, ncol = length(methodes))
accPlot = matrix(0, nrow = nbmachines, ncol = length(all_methodes))
temps   <- matrix(0, nrow = nbmachines, ncol = length(methodes))

for (j in machines_ok_final_boucle) {
  
  
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
}############################# LABELS #################################
 ariPlot <- ariPlot[, -c(7, 14, 21), drop = FALSE]
 accPlot <- accPlot[, -c(7, 14, 21), drop = FALSE]
 aucPlot <- aucPlot[, -7, drop = FALSE]
setwd(fig_SMD)

names_plot_all <- c(
  # QC
  "Sple QC",
  "Onl QC",
  "Strm QC",
  "Offl QC",
  "OGK QC",
  "MCD QC",
  
  # RD
  "Sple RD",
  "Onl RD",
  "Strm RD",
  "Offl RD",
  "OGK RD",
  "MCD RD",
  
  # Raw
  "Sple Raw",
  "Onl Raw",
  "Strm Raw",
  "Offl Raw",
  "OGK Raw",
  "MCD Raw"
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
  
  # RD
  "darkgreen",
  "blue",
  "red",
  "orange",
  "brown",
  "black",
  
  # Raw
  "darkgreen",
  "blue",
  "red",
  "orange",
  "brown",
  "black"
)


############################# ARI #################################

pdf(
  "boxplotARI_SMD.pdf",
  width = 20,
  height = 7
)

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

pdf(
  "boxplotACC_SMD.pdf",
  width = 20,
  height = 7
)

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

############################# AUC #################################

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

pdf(
  "boxplotAUC_SMD.pdf",
  width = 12,
  height = 6
)

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

setwd(fig_SMD)

#=========================================================
# 2. Variables dédiées : méthodes, conditions, ordre, couleurs
#=========================================================

methodes_base <- c("Sple", "Onl", "Strm", "Offl", "OGK", "MCD")
conditions    <- c("QC", "RD", "Raw")

# Indices des colonnes pour chaque condition (après suppression 7,14,21)
idx_qc  <- 1:6      # Sple QC, Onl QC, Strm QC, Offl QC, OGK QC, MCD QC
idx_rd  <- 7:12     # Sple RD, Onl RD, Strm RD, Offl RD, OGK RD, MCD RD
idx_raw <- 13:18    # Sple Raw, Onl Raw, Strm Raw, Offl Raw, OGK Raw, MCD Raw

# Ordre voulu : Sple QC, Sple RD, Sple Raw, Onl QC, Onl RD, Onl Raw, ...
ordre_methodes <- as.vector(rbind(idx_qc, idx_rd, idx_raw))

# Noms correspondants
names_plot_all <- as.vector(t(outer(methodes_base, conditions, paste)))
# -> "Sple QC" "Sple RD" "Sple Raw" "Onl QC" "Onl RD" "Onl Raw" ...

# Couleurs d'origine par MÉTHODE, répétées pour QC / RD / Raw
cols_methodes <- c("darkgreen", "blue", "red", "orange", "brown", "black")
cols_all      <- rep(cols_methodes, each = length(conditions))

# Positions avec un petit gap entre chaque trio de méthode
at_positions <- as.vector(
  sapply(seq_along(methodes_base),
         function(k) (k - 1) * 3.5 + 1:3)
)


#=========================================================
# 3. Réordonnancement des matrices
#=========================================================

ariPlot_ord <- ariPlot[, ordre_methodes, drop = FALSE]
accPlot_ord <- accPlot[, ordre_methodes, drop = FALSE]

#=========================================================
# 4. Boxplot ARI
#=========================================================

pdf("boxplotARI_SMD.pdf", width = 20, height = 7)

boxplot(
  ariPlot_ord,
  at       = at_positions,
  names    = names_plot_all,
  las      = 2,
  col      = cols_all,
  main     = "ARI",
  ylab     = "ARI",
  cex.axis = 1.2,
  xlim     = c(0.5, max(at_positions) + 0.5)
)

legend(
  "bottomright",
  legend = methodes_base,
  fill   = cols_methodes,
  bty    = "n",
  cex    = 1.1
)

dev.off()


#=========================================================
# 5. Boxplot Accuracy
#=========================================================

pdf("boxplotACC_SMD.pdf", width = 20, height = 7)

boxplot(
  accPlot_ord,
  at       = at_positions,
  names    = names_plot_all,
  las      = 2,
  col      = cols_all,
  main     = "Accuracy",
  ylab     = "Accuracy",
  cex.axis = 1.2,
  xlim     = c(0.5, max(at_positions) + 0.5)
)

# legend(
#   "bottomright",
#   legend = "",
#   fill   = cols_methodes,
#   bty    = "n",
#   cex    = 1.1
# )

dev.off()

#### Calcul des trajectoires ####
machines_ok_final_boucle <- machines_ok_final[!machines_ok_final %in% c(6, 7,15,27,34,35,43,44,45,47,55)]
methodes_online_quantile = c(
  "SampleNaiveQuantonlinecorr",
  "OnlineUsQuantonlinecorr",
  "StreamingUsonlineQuantcorr",
  "OracleQC","MCDonlineQC","OGKonlineQC"
)

methodes_online_rescale = c(
  "SampleNaivewithoutonlinequantilecorr",
  "OnlineUswithoutQuantonlinecorr",
  "StreamingUswithoutQuantonlinecorr",
  "OracleRD"
)

methodes_online_raw = c(
  "SampleRaw",
  "OnlRaw",
  "StrmRaw",
  "Oracle"
)

methodes_online = c(
  methodes_online_quantile,
  methodes_online_rescale,
  methodes_online_raw
)
# Symboles associés aux familles
pch_methodes = rep(NA,length(methodes_online))

# Correction quantile : étoile
pch_methodes[methodes_online %in% methodes_online_quantile] = 8

# Rescale distance : carré
pch_methodes[methodes_online %in% methodes_online_rescale] = 15

# Raw : triangle
pch_methodes[methodes_online %in% methodes_online_raw] = 17





for(j in machines_ok_final_boucle){
  
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
    
    
    if(methode == "StreamingUsonlineQuantcorr") {
         
        cum_times_strm = resultats$cum_times 
        cum_times_strm = cum_times_strm[1:(length(cum_times_strm) - 1)]
    }
    
    if(methode == "MCDonlineQC") {
      cum_times_mcd = resultats$cum_times
    }
    
    if(methode == "OGKonlineQC") {
      cum_times_ogk = resultats$cum_times
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
  
  
  rates_mcd_online_qc <-
    compute_rates(
      outlmach[, "MCDonlineQC"],
      labels
    )
  
  rates_ogk_online_qc <-
    compute_rates(
      outlmach[, "OGKonlineQC"],
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
  
  
  nom_fichier_box <- paste0("boxplot_SMD_mach-", j, ".pdf")
  
  pdf(nom_fichier_box, width = 6, height = 5)
  par(mar = c(4, 4, 2, 1))

  
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
    col = c("lightblue", "red"),
    names = c("inliers", "outliers"),
    log = "y",
    yaxt = "n",
    ylab = ""
  )
  
  # Détermination automatique des puissances nécessaires
  ymin <- min(c(distinliers, distoutliers), na.rm = TRUE)
  ymax <- max(c(distinliers, distoutliers), na.rm = TRUE)
  
  pmin <- floor(log10(ymin))
  pmax <- ceiling(log10(ymax))
  
  axis(
    2,
    at = 10^(pmin:pmax),
    labels = parse(
      text = paste0("10^", pmin:pmax)
    ),
    las = 1
  )
  
  dev.off()
  
  
  ##################################
  # 2. FALSE NEGATIVE RATE
  ##################################
  nom_fichier_fn <- paste0("FN_SMD_mach-", j, ".pdf")
  
  pdf(nom_fichier_fn, width = 8, height = 6)
  par(mar = c(4, 4, 2, 1))
  
  plot(
    x_vals,
    rates_Strm_with_quantcorr$FN_rate * 100,
    type = "l",
    lwd = 3,
    col = "red",
    ylim = c(0,100),
    xlab = "",
    ylab = "",
    main = "",
    xaxt = "n",
    yaxt = "n"
  )
  
  
  #=========================================================
  # QC
  #=========================================================
  
  # 
  # # Sample QC
  # lines(
  #   x_vals,
  #   rates_samplecov_quantcorr$FN_rate * 100,
  #   lty = "dotted",
  #   col = "darkgreen",
  #   lwd = 3
  # )
  # 
  # points(
  #   x_vals[idx_symbols],
  #   rates_samplecov_quantcorr$FN_rate[idx_symbols] * 100,
  #   pch = pch_methodes[
  #     "SampleNaiveQuantonlinecorr"
  #   ],
  #   col = "darkgreen",
  #   cex = 1.4
  # )
  # 
  
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
    pch = 8,
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
    pch = 8,
    col = "red",
    cex = 1.4
  )
  # 
  # # # -------------------------
  # # # MCD Online raw
  # # # -------------------------
  # # 
   lines(
      x_vals,
     rates_mcd_online_qc$FN_rate * 100,
      lty = "dotdash",
      col = "black",
      lwd = 3
    )
   # 
    points(
      x_vals[idx_symbols],
      rates_mcd_online_qc$FN_rate[idx_symbols] * 100,
      pch = 17,
      col = "black",
      cex = 1.8
    )
  # # 
  # # 
  # # # -------------------------
  # # # OGK Online Raw
  # # # -------------------------
  # # 
    lines(
      x_vals,
      rates_ogk_online_qc$FN_rate * 100,
      lty = "twodash",
      col = "brown",
      lwd = 3
    )
  # # 
    points(
      x_vals[idx_symbols],
      rates_ogk_online_qc$FN_rate[idx_symbols] * 100,
      pch = 17,
      col = "brown",
      cex = 1.4
    )
  # # 
  
  # 
  # 
  # # Oracle QC
  # lines(
  #   x_vals,
  #   rates_oracle_qc$FN_rate * 100,
  #   lty = "dashed",
  #   col = "purple",
  #   lwd = 3
  # )
  # 
  # points(
  #   x_vals[idx_symbols],
  #   rates_oracle_qc$FN_rate[idx_symbols] * 100,
  #   pch = pch_methodes["OracleQC"],
  #   col = "purple",
  #   cex = 1.4
  # )
  # 
  # 
  # 
  #=========================================================
  # RD
  #=========================================================
  # 
  # 
  # # Sample RD
  # lines(
  #   x_vals,
  #   rates_samplecov_without_quantcorr$FN_rate * 100,
  #   lty = "dotdash",
  #   col = "darkgreen",
  #   lwd = 3
  # )
  # 
  # points(
  #   x_vals[idx_symbols],
  #   rates_samplecov_without_quantcorr$FN_rate[idx_symbols] * 100,
  #   pch = pch_methodes[
  #     "SampleNaivewithoutonlinequantilecorr"
  #   ],
  #   col = "darkgreen",
  #   cex = 1.4
  # )
  # 
  # 
  # # Online RD
   lines(
     x_vals,
    rates_online_without_quantcorr$FN_rate * 100,
    lty = "twodash",
     col = "blue",
     lwd = 3
   )
  # 
   points(
     x_vals[idx_symbols],
     rates_online_without_quantcorr$FN_rate[idx_symbols] * 100,
     pch = 15,
     col = "blue",
     cex = 1.8
   )
  # 
  # 
  # # Streaming RD
   lines(
     x_vals,
     rates_Strm_without_quantcorr$FN_rate * 100,
     lty = "longdash",
     col = "red",
     lwd = 3
   )
  # 
   points(
     x_vals[idx_symbols],
     rates_Strm_without_quantcorr$FN_rate[idx_symbols] * 100,
     pch = 15,
     col = "red",
   cex = 1.8
   )
  # # 
  # # 
  # # # Oracle RD
  # # lines(
  # #   x_vals,
  # #   rates_oracle_rd$FN_rate * 100,
  # #   lty = "longdash",
  # #   col = "purple",
  # #   lwd = 3
  # # )
  # # 
  # # points(
  # #   x_vals[idx_symbols],
  # #   rates_oracle_rd$FN_rate[idx_symbols] * 100,
  # #   pch = pch_methodes["OracleRD"],
  # #   col = "purple",
  # #   cex = 1.4
  # # )
  # # 
  # # 
  # 
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
    pch = 17,
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
    pch = 17,
    col = "blue",
    cex = 1.8
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
    pch = 17,
    col = "red",
    cex = 1.8
  )
  
  # 
  # # Oracle Raw
  # lines(
  #   x_vals,
  #   rates_oracle_raw$FN_rate * 100,
  #   lty = "solid",
  #   col = "purple",
  #   lwd = 3
  # )
  # 
  # points(
  #   x_vals[idx_symbols],
  #   rates_oracle_raw$FN_rate[idx_symbols] * 100,
  #   pch = pch_methodes["Oracle"],
  #   col = "purple",
  #   cex = 1.4
  # )
  # 
  # 
  
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
  
  
  ##################################
  # 3. FALSE POSITIVE RATE
  ##################################
  
  nom_fichier_fp <- paste0("FP_SMD_mach-", j, ".pdf")
  
  pdf(nom_fichier_fp, width = 8, height = 6)
  
  par(mar = c(4, 4, 2, 1))
  
  plot(
    x_vals,
    rates_Strm_with_quantcorr$FP_rate * 100,
    type = "l",
    lwd = 3,
    col = "red",
    ylim = c(0,100),
    xlab = "",
    ylab = "",
    main = "",
    xaxt = "n",
    yaxt = "n"
  )
  
  
  #=========================================================
  # QC
  #=========================================================
  # 
  # 
  # # Sample QC
  # lines(
  #   x_vals,
  #   rates_samplecov_quantcorr$FP_rate * 100,
  #   lty = "dotted",
  #   col = "darkgreen",
  #   lwd = 3
  # )
  # 
  # points(
  #   x_vals[idx_symbols],
  #   rates_samplecov_quantcorr$FP_rate[idx_symbols] * 100,
  #   pch = pch_methodes[
  #     "SampleNaiveQuantonlinecorr"
  #   ],
  #   col = "darkgreen",
  #   cex = 1.4
  # )
  # 
  
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
    pch = 8,
    col = "blue",
    cex = 1.8
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
    pch = 8,
    col = "red",
    cex = 1.8
  )
  
  # -------------------------
  # MCD Online Raw
  # -------------------------
  # 
  lines(
     x_vals,
     rates_mcd_online_qc$FP_rate * 100,
     lty = "dotdash",
     col = "black",
     lwd = 3
   )
  # 
   points(
     x_vals[idx_symbols],
     rates_mcd_online_qc$FP_rate[idx_symbols] * 100,
     pch = 17,
     col = "black",
     cex = 1.8
   )
  # 
  # 
  # # -------------------------
  # # OGK Online QC
  # # -------------------------
  # 
  lines(
     x_vals,
     rates_ogk_online_qc$FP_rate * 100,
     lty = "twodash",
     col = "brown",
     lwd = 3
   )
  # 
   points(
     x_vals[idx_symbols],
     rates_ogk_online_qc$FP_rate[idx_symbols] * 100,
     pch = 17,
     col = "brown",
     cex = 1.8
   )
  # 
  # 
  
  # 
  # # Oracle QC
  # lines(
  #   x_vals,
  #   rates_oracle_qc$FP_rate * 100,
  #   lty = "dashed",
  #   col = "purple",
  #   lwd = 3
  # )
  # 
  # points(
  #   x_vals[idx_symbols],
  #   rates_oracle_qc$FP_rate[idx_symbols] * 100,
  #   pch = pch_methodes["OracleQC"],
  #   col = "purple",
  #   cex = 1.4
  # )
  # 
  # 
  
  #=========================================================
  # RD
  #=========================================================
  # 
  # 
  # # Sample RD
  # lines(
  #   x_vals,
  #   rates_samplecov_without_quantcorr$FP_rate * 100,
  #   lty = "dotdash",
  #   col = "darkgreen",
  #   lwd = 3
  # )
  # 
  # points(
  #   x_vals[idx_symbols],
  #   rates_samplecov_without_quantcorr$FP_rate[idx_symbols] * 100,
  #   pch = pch_methodes[
  #     "SampleNaivewithoutonlinequantilecorr"
  #   ],
  #   col = "darkgreen",
  #   cex = 1.4
  # )
  # 
  # 
  # # Online RD
  lines(
     x_vals,
     rates_online_without_quantcorr$FP_rate * 100,
     lty = "twodash",
     col = "blue",
     lwd = 3
   )
  # 
   points(
     x_vals[idx_symbols],
     rates_online_without_quantcorr$FP_rate[idx_symbols] * 100,
     pch = 15,
    col = "blue",
     cex = 1.8 )
  # 
  # 
  # # Streaming RD
  lines(
     x_vals,
     rates_Strm_without_quantcorr$FP_rate * 100,
     lty = "longdash",
     col = "red",
     lwd = 3
   )
  # 
  points(
     x_vals[idx_symbols],
     rates_Strm_without_quantcorr$FP_rate[idx_symbols] * 100,
     pch = 15,
  col = "red",
     cex = 1.8
 )
  # 
  # 
  # # Oracle RD
  # lines(
  #   x_vals,
  #   rates_oracle_rd$FP_rate * 100,
  #   lty = "longdash",
  #   col = "purple",
  #   lwd = 3
  # )
  # 
  # points(
  #   x_vals[idx_symbols],
  #   rates_oracle_rd$FP_rate[idx_symbols] * 100,
  #   pch = pch_methodes["OracleRD"],
  #   col = "purple",
  #   cex = 1.4
  # )
  # 
  # 
  
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
    pch = 17,
    col = "darkgreen",
    cex = 1.8
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
    pch = 17,
    col = "blue",
    cex = 1.8
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
    pch = 17,
    col = "red",
    cex = 1.8
  )
  # 
  # 
  # # Oracle Raw
  # lines(
  #   x_vals,
  #   rates_oracle_raw$FP_rate * 100,
  #   lty = "solid",
  #   col = "purple",
  #   lwd = 3
  # )
  # 
  # points(
  #   x_vals[idx_symbols],
  #   rates_oracle_raw$FP_rate[idx_symbols] * 100,
  #   pch = pch_methodes["Oracle"],
  #   col = "purple",
  #   cex = 1.4
  # )
  # 
  
  
  axis(
    2,
    las = 1,
    at = seq(0, 100, by = 5),
    cex.axis = 1.5
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
  
  ##################################
  # 4. TEMPS CUMULES
  ##################################
  # 
  # #####Ajout du temps d'attente entre chaque itération = 1 minute
  # 
 #  cum_times_strm = cum_times_strm + x_vals
  # 
  # cum_times_mcd = cum_times_mcd + x_vals * 60
  # 
  # cum_times_ogk = cum_times_ogk + x_vals * 60
  #n_init = 1e3
  #batch = 10
  # Nombre d'observations traitées par MCD à chaque itération streaming
  # n_obs_mcd <- rep(0, length(x_vals))
  # mask <- x_vals >= n_init
  # n_obs_mcd[mask] <- n_init + floor((x_vals[mask] - n_init) / batch) * batch
  # 
  # # Itérations où un nouveau batch MCD est déclenché
  # new_batch <- c(FALSE, diff(n_obs_mcd) > 0)
  # 
  # # Attente à chaque itération, en minutes
  # attente_s <- numeric(length(x_vals))
  # attente_s[new_batch] <- n_obs_mcd[new_batch] 
  # 
  # # Cumul de l'attente, en SECONDES
  # attente_cumulee_s <- cumsum(attente_s)
  # 
  # Ajout au temps de calcul MCD / OGK
  #cum_times_mcd <- cum_times_mcd/60 + attente_cumulee_s
  #cum_times_ogk <- cum_times_ogk/60 + attente_cumulee_s
  nom_fichier_cumtimes<- paste0("cum_times_SMD_mach-", j, ".pdf")


  pdf(  nom_fichier_cumtimes, width = 8, height = 6)
  par(mar = c(4, 4, 2, 1))
  
  # Minimum et maximum globaux des ordonnées
  y_min <- min(
    cum_times_strm,
    cum_times_mcd,
    cum_times_ogk,
    na.rm = TRUE
  )
  
  y_max <- max(
    cum_times_strm,
    cum_times_mcd,
    cum_times_ogk,
    na.rm = TRUE
  )
  
  # # Puissances de 10 nécessaires
  # pow_min <- floor(log10(y_min))
  # pow_max <- ceiling(log10(y_max))
  y_max <- max(cum_times_strm, cum_times_mcd, cum_times_ogk, na.rm = TRUE)
  
  #y_ticks <- 10^(pow_min:pow_max)
  
  plot(
    x_vals,
    cum_times_strm,
    type = "l",
    lwd = 3,
    col = "red",
    xlab = "",
    ylab = "",
    main = "",
    xaxt = "n",
    yaxt = "n",
    #log = "y",
    ylim  = c(0, y_max * 1.05)
  )
  
  lines(
    x_vals,
    cum_times_mcd,
    lwd = 3,
    col = "black",
    lty = "dotdash"
  )
  
  lines(
    x_vals,
    cum_times_ogk,
    lwd = 3,
    col = "brown",
    lty = "twodash"
  )
  
  points(
    x_vals[idx_symbols],
    cum_times_strm[idx_symbols],
    pch = 17,
    col = "red",
    cex = 1.8
  )
  
  points(
    x_vals[idx_symbols],
    cum_times_mcd[idx_symbols],
    pch = 17,
    col = "black",
    cex = 1.8
  )
  
  points(
    x_vals[idx_symbols],
    cum_times_ogk[idx_symbols],
    pch = 17,
    col = "brown",
    cex = 1.8
  )
  
  axis(
    1,
    at = seq(1000, max(x_vals), by = 1000),
    las = 1,
    cex.axis = 1.8
  )
  
  # Axe Y avec 10^-2, 10^-1, 10^0, 10^1, ...
  # axis(
  #   2,
  #   at = y_ticks,
  #   labels = parse(text = paste0("10^", pow_min:pow_max)),
  #   las = 1,
  #   cex.axis = 1.8
  # )
  # 
  
  axis(2, las = 1, cex.axis = 1.5)
  
  dev.off()
}
