#### Calcul des trajectoires ####
machines_ok_final = machines_ok_final_boucle
machines_ok_final_boucle <- machines_ok_final[
  !machines_ok_final %in% c(6, 7, 15, 27, 34, 35, 43, 44, 45, 47, 55)
]

#############################
# Config visuelle unique
#############################
setwd("~/work")
load("methods_colors.RData")   # restaure methods_df

# Petit pont SMD -> methods_df (les noms diffèrent pour MCD/OGK en ligne)
normalize_method <- function(m) {
  mapping <- c(
    "MCDonlineQC" = "MCDQC",
    "OGKonlineQC" = "OGKQC",
    "MCDonlineRD" = "MCDRD",
    "OGKonlineRD" = "OGKRD"
  )
  if (m %in% names(mapping)) return(mapping[[m]])
  m
}

# Helpers d'accès : couleur / lty / pch par nom de méthode
get_color <- function(m) {
  m2 <- normalize_method(m)
  if (m2 %in% methods_df$method)
    return(methods_df$color_alpha[match(m2, methods_df$method)])
  "gray50"
}
get_lty <- function(m) {
  m2 <- normalize_method(m)
  if (m2 %in% methods_df$method)
    return(methods_df$lty[match(m2, methods_df$method)])
  1
}
get_pch <- function(m) {
  m2 <- normalize_method(m)
  if (m2 %in% methods_df$method)
    return(methods_df$pch[match(m2, methods_df$method)])
  NA
}

#############################
# Liste des méthodes
#############################
methodes_online_quantile = c(
  "SampleNaiveQuantonlinecorr",
  "OnlineUsQuantonlinecorr",
  "StreamingUsonlineQuantcorr",
  "OracleQC", "MCDonlineQC", "OGKonlineQC"
)

methodes_online_rescale = c(
  "SampleNaivewithoutonlinequantilecorr",
  "OnlineUswithoutQuantonlinecorr",
  "StreamingUswithoutQuantonlinecorr",
  "OracleRD"
)

methodes_online_raw = c(
  "SampleRaw", "OnlRaw", "StrmRaw", "Oracle"
)

methodes_online = c(
  methodes_online_quantile,
  methodes_online_rescale,
  methodes_online_raw
)

#############################
# Boucle machines
#############################
for (j in machines_ok_final_boucle) {
  
  setwd(smd_data_dir)
  load(paste0("data_machine-", j, ".RData"))
  
  outlmach <- matrix(0, nrow = nrow(Z), ncol = length(methodes_online))
  colnames(outlmach) <- methodes_online
  
  setwd(res_SMD)
  distoracle <- rep(0, nrow(Z))
  
  #=========================================================
  # Chargement des résultats
  #=========================================================
  for (s in seq_along(methodes_online)) {
    methode <- methodes_online[s]
    load(paste0("Fit-", methode, "-machine-", j, ".RData"))
    outlmach[, s] <- resultats$outliers_labels
    
    if (methode == "Oracle") distoracle <- resultats$distances
    if (methode == "StreamingUsonlineQuantcorr") {
      cum_times_strm <- resultats$cum_times
      cum_times_strm <- cum_times_strm[1:(length(cum_times_strm) - 1)]
    }
    if (methode == "MCDonlineQC") cum_times_mcd <- resultats$cum_times
    if (methode == "OGKonlineQC") cum_times_ogk <- resultats$cum_times
  }
  
  #=========================================================
  # Calcul des taux
  #=========================================================
  # QC
  rates_samplecov_quantcorr  <- compute_rates(outlmach[, "SampleNaiveQuantonlinecorr"], labels)
  rates_online_with_quantcorr<- compute_rates(outlmach[, "OnlineUsQuantonlinecorr"],    labels)
  rates_Strm_with_quantcorr  <- compute_rates(outlmach[, "StreamingUsonlineQuantcorr"], labels)
  rates_mcd_online_qc        <- compute_rates(outlmach[, "MCDonlineQC"],                labels)
  rates_ogk_online_qc        <- compute_rates(outlmach[, "OGKonlineQC"],                labels)
  rates_oracle_qc            <- compute_rates(outlmach[, "OracleQC"],                   labels)
  
  # RD
  rates_samplecov_without_quantcorr <- compute_rates(outlmach[, "SampleNaivewithoutonlinequantilecorr"], labels)
  rates_online_without_quantcorr    <- compute_rates(outlmach[, "OnlineUswithoutQuantonlinecorr"],       labels)
  rates_Strm_without_quantcorr      <- compute_rates(outlmach[, "StreamingUswithoutQuantonlinecorr"],    labels)
  rates_oracle_rd                   <- compute_rates(outlmach[, "OracleRD"],                             labels)
  
  # Raw
  rates_samplecov_raw <- compute_rates(outlmach[, "SampleRaw"], labels)
  rates_online_raw    <- compute_rates(outlmach[, "OnlRaw"],    labels)
  rates_Strm_raw      <- compute_rates(outlmach[, "StrmRaw"],   labels)
  rates_oracle_raw    <- compute_rates(outlmach[, "Oracle"],    labels)
  
  #=========================================================
  # Préparation du tracé
  #=========================================================
  setwd(fig_SMD)
  x_vals <- 1:length(rates_Strm_with_quantcorr$FN_rate)
  idx_outliers <- which(labels == 1)   # positions des vraies anomalies
  
  ##################################
  # 1. BOXPLOT DISTANCES
  ##################################
  pdf(paste0("boxplot_SMD_mach-", j, ".pdf"), width = 6, height = 5)
  par(mar = c(4, 4, 2, 1))
  
  Z_clean    <- Z[labels == 0, , drop = FALSE]
  Z_outliers <- Z[labels == 1, , drop = FALSE]
  
  distinliers  <- rep(0, nrow(Z_clean))
  distoutliers <- rep(0, nrow(Z_outliers))
  
  invSigmaTrueCov <- solve(cov(Z_clean))
  mu <- colMeans(Z_clean)
  
  for (m in 1:nrow(Z_clean)) {
    diff <- Z_clean[m, ] - mu
    distinliers[m] <- t(diff) %*% invSigmaTrueCov %*% diff
  }
  for (m in 1:nrow(Z_outliers)) {
    diff <- Z_outliers[m, ] - mu
    distoutliers[m] <- t(diff) %*% invSigmaTrueCov %*% diff
  }
  
  boxplot(
    distinliers, distoutliers,
    col   = c("lightblue", "red"),
    names = c("inliers", "outliers"),
    log   = "y", yaxt = "n", ylab = ""
  )
  
  ymin <- min(c(distinliers, distoutliers), na.rm = TRUE)
  ymax <- max(c(distinliers, distoutliers), na.rm = TRUE)
  pmin <- floor(log10(ymin))
  pmax <- ceiling(log10(ymax))
  
  axis(2,
       at = 10^(pmin:pmax),
       labels = parse(text = paste0("10^", pmin:pmax)),
       las = 1)
  dev.off()
  
  ##################################
  # 2. FALSE NEGATIVE RATE
  ##################################
  pdf(paste0("FN_SMD_mach-", j, ".pdf"), width = 8, height = 6)
  par(mar = c(4, 4, 2, 1))
  
  plot(
    x_vals,
    rates_Strm_without_quantcorr$FN_rate * 100,
    type = "l", lwd = 6,
    col  = get_color("StreamingUsonlineQuantcorr"),
    lty  = get_lty("StreamingUsonlineQuantcorr"),
    ylim = c(0, 100),
    xlab = "", ylab = "", main = "",
    xaxt = "n", yaxt = "n"
  )
  
  # Traits verticaux aux positions des vraies anomalies
  abline(v = idx_outliers,
         col = adjustcolor("grey50", alpha.f = 0.25),
         lwd = 0.5)
  
  # QC (optionnel, commenté)
  # lines(x_vals, rates_samplecov_quantcorr$FN_rate * 100,
  #       lty = get_lty("SampleNaiveQuantonlinecorr"),
  #       col = get_color("SampleNaiveQuantonlinecorr"), lwd = 3)
  # lines(x_vals, rates_online_with_quantcorr$FN_rate * 100,
  #       lty = get_lty("OnlineUsQuantonlinecorr"),
  #       col = get_color("OnlineUsQuantonlinecorr"), lwd = 3)
  # lines(x_vals, rates_Strm_with_quantcorr$FN_rate * 100,
  #       lty = get_lty("StreamingUsonlineQuantcorr"),
  #       col = get_color("StreamingUsonlineQuantcorr"), lwd = 3)
  
  lines(x_vals, rates_mcd_online_qc$FN_rate * 100,
        lty = get_lty("MCDonlineQC"),
        col = get_color("MCDonlineQC"), lwd = 6)
  
  lines(x_vals, rates_ogk_online_qc$FN_rate * 100,
        lty = get_lty("OGKonlineQC"),
        col = get_color("OGKonlineQC"), lwd = 6)
  
  # lines(x_vals, rates_oracle_qc$FN_rate * 100,
  #       lty = get_lty("OracleQC"),
  #       col = get_color("OracleQC"), lwd = 3)
  
  # RD
  lines(x_vals, rates_samplecov_without_quantcorr$FN_rate * 100,
        lty = get_lty("SampleNaivewithoutonlinequantilecorr"),
        col = get_color("SampleNaivewithoutonlinequantilecorr"),
        lwd = 6)
  
  lines(x_vals, rates_online_without_quantcorr$FN_rate * 100,
        lty = get_lty("OnlineUswithoutQuantonlinecorr"),
        col = get_color("OnlineUswithoutQuantonlinecorr"),
        lwd = 6)
  
  # lines(x_vals, rates_Strm_without_quantcorr$FN_rate * 100,
  #       lty = get_lty("StreamingUswithoutQuantonlinecorr"),
  #       col = get_color("StreamingUswithoutQuantonlinecorr"), lwd = 6)
  # lines(x_vals, rates_oracle_rd$FN_rate * 100,
  #       lty = get_lty("OracleRD"),
  #       col = get_color("OracleRD"), lwd = 3)
  
  # Raw (commenté)
  # lines(x_vals, rates_samplecov_raw$FN_rate * 100,
  #       lty = get_lty("SampleRaw"),
  #       col = get_color("SampleRaw"), lwd = 3)
  # lines(x_vals, rates_online_raw$FN_rate * 100,
  #       lty = get_lty("OnlRaw"),
  #       col = get_color("OnlRaw"), lwd = 3)
  # lines(x_vals, rates_Strm_raw$FN_rate * 100,
  #       lty = get_lty("StrmRaw"),
  #       col = get_color("StrmRaw"), lwd = 3)
  # lines(x_vals, rates_oracle_raw$FN_rate * 100,
  #       lty = get_lty("Oracle"),
  #       col = get_color("Oracle"), lwd = 3)
  
  axis(2, las = 1, cex.axis = 1.8)
  axis(1, at = seq(1000, max(x_vals), by = 1000), las = 1, cex.axis = 1.8)
  box()
  dev.off()
  
  ##################################
  # 3. FALSE POSITIVE RATE
  ##################################
  pdf(paste0("FP_SMD_mach-", j, ".pdf"), width = 8, height = 6)
  par(mar = c(4, 4, 2, 1))
  
  plot(
    x_vals,
    rates_Strm_without_quantcorr$FP_rate * 100,
    type = "l", lwd = 6,
    col  = get_color("StreamingUswithoutQuantonlinecorr"),
    lty  = get_lty("StreamingUswithoutQuantonlinecorr"),
    ylim = c(0, 100),
    xlab = "", ylab = "", main = "",
    xaxt = "n", yaxt = "n"
  )
  
  # Traits verticaux aux positions des vraies anomalies
  abline(v = idx_outliers,
         col = adjustcolor("grey50", alpha.f = 0.25),
         lwd = 0.5)
  
  # QC (commenté)
  # lines(x_vals, rates_samplecov_quantcorr$FP_rate * 100,
  #       lty = get_lty("SampleNaiveQuantonlinecorr"),
  #       col = get_color("SampleNaiveQuantonlinecorr"), lwd = 3)
  # lines(x_vals, rates_online_with_quantcorr$FP_rate * 100,
  #       lty = get_lty("OnlineUsQuantonlinecorr"),
  #       col = get_color("OnlineUsQuantonlinecorr"), lwd = 3)
  # lines(x_vals, rates_Strm_with_quantcorr$FP_rate * 100,
  #       lty = get_lty("StreamingUsonlineQuantcorr"),
  #       col = get_color("StreamingUsonlineQuantcorr"), lwd = 3)
  # lines(x_vals, rates_oracle_qc$FP_rate * 100,
  #       lty = get_lty("OracleQC"),
  #       col = get_color("OracleQC"), lwd = 3)
  
  # RD
  lines(x_vals, rates_samplecov_without_quantcorr$FP_rate * 100,
        lty = get_lty("SampleNaivewithoutonlinequantilecorr"),
        col = get_color("SampleNaivewithoutonlinequantilecorr"), lwd = 6)
  
  lines(x_vals, rates_online_without_quantcorr$FP_rate * 100,
        lty = get_lty("OnlineUswithoutQuantonlinecorr"),
        col = get_color("OnlineUswithoutQuantonlinecorr"), lwd = 6)
  
  lines(x_vals, rates_mcd_online_qc$FP_rate * 100,
        lty = get_lty("MCDonlineQC"),
        col = get_color("MCDonlineQC"), lwd = 6)
  
  # lines(x_vals, rates_oracle_rd$FP_rate * 100,
  #       lty = get_lty("OracleRD"),
  #       col = get_color("OracleRD"), lwd = 3)
  
  # Raw (commenté)
  # lines(x_vals, rates_samplecov_raw$FP_rate * 100,
  #       lty = get_lty("SampleRaw"),
  #       col = get_color("SampleRaw"), lwd = 3)
  # lines(x_vals, rates_online_raw$FP_rate * 100,
  #       lty = get_lty("OnlRaw"),
  #       col = get_color("OnlRaw"), lwd = 3)
  # lines(x_vals, rates_Strm_raw$FP_rate * 100,
  #       lty = get_lty("StrmRaw"),
  #       col = get_color("StrmRaw"), lwd = 3)
  
  lines(x_vals, rates_ogk_online_qc$FP_rate * 100,
        lty = get_lty("OGKonlineQC"),
        col = get_color("OGKonlineQC"), lwd = 6)
  
  # lines(x_vals, rates_oracle_raw$FP_rate * 100,
  #       lty = get_lty("Oracle"),
  #       col = get_color("Oracle"), lwd = 6)
  
  axis(2, las = 1, at = seq(0, 100, by = 5), cex.axis = 1.5)
  axis(1, at = seq(1000, max(x_vals), by = 1000), las = 1, cex.axis = 1.8)
  box()
  dev.off()
  
  ##################################
  # 4. TEMPS CUMULES
  ##################################
  pdf(paste0("cum_times_SMD_mach-", j, ".pdf"), width = 8, height = 6)
  par(mar = c(4, 4, 2, 1))
  
  y_min <- min(cum_times_strm, cum_times_mcd, cum_times_ogk, na.rm = TRUE)
  y_max <- max(cum_times_strm, cum_times_mcd, cum_times_ogk, na.rm = TRUE)
  
  pow_min <- floor(log10(y_min))
  pow_max <- ceiling(log10(y_max))
  y_ticks <- 10^(pow_min:pow_max)
  
  plot(
    x_vals,
    cum_times_strm,
    type = "l", lwd = 6,
    col  = get_color("StreamingUsonlineQuantcorr"),
    lty  = get_lty("StreamingUsonlineQuantcorr"),
    xlab = "", ylab = "", main = "",
    xaxt = "n", yaxt = "n",
    log  = "y",
    ylim = c(10^pow_min, 10^pow_max)
  )
  
  # Traits verticaux aux positions des vraies anomalies
  # abline(v = idx_outliers,
  #        col = adjustcolor("grey50", alpha.f = 0.25),
  #        lwd = 0.5)
  # 
  lines(x_vals, cum_times_mcd,
        lwd = 6,
        lty = get_lty("MCDonlineQC"),
        col = get_color("MCDonlineQC"))
  
  lines(x_vals, cum_times_ogk,
        lwd = 6,
        lty = get_lty("OGKonlineQC"),
        col = get_color("OGKonlineQC"))
  
  axis(1, at = seq(1000, max(x_vals), by = 1000), las = 1, cex.axis = 1.8)
  
  axis(2,
       at = y_ticks,
       labels = parse(text = paste0("10^", pow_min:pow_max)),
       las = 1,
       cex.axis = 1.8)
  
  box()
  dev.off()
}
