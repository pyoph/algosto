smd_data_dir = "~/smd_data_dir"
res_SMD = "~/res_SMD"
crit_SMD = "~/criteres_SMD"
fig_SMD = "~/fig_SMD/"

########################Boxplot Sigma errors##############################

methodes = c("Oracle","SampleNaiveQuantonlinecorr","SampleNaivewithoutQuantonlinecorr","MCD","Offline-without_onlinequantile","Offline-without_onlinequantile","OnlineUsQuantonlinecorr","OnlineUsWithoutQuantonlinecorr","StreamingUsonlineQuantcorr","StreamingUs_without_onlineQuantcorr","OGK")


erreursSigmaSMD = array(0, dim = c(length(methodes),28))  


for(j in 1:28){

setwd(crit_SMD)
  
  
  
  for(s in seq_along(methodes)){
    
    methode = methodes[s]
    
    critFile <- paste0('Crit-',methode,"-machine-",j,".RData")
    
    load(critFile)
    
    erreursSigmaSMD[s,j] = crit$erreurFrob
    
    
  }
  
}

setwd(fig_SMD)

file <- paste0("boxploterreursSigma_SMD", ".pdf")

pdf(file, width = 18, height = 6)

# -----------------------------
# transformation
# -----------------------------
erreurs_log <- log1p(erreursSigmaSMD)

# -----------------------------
# données
# -----------------------------
data_plot <- list(
  Streaming   = erreurs_log[8,],
  MCD         = erreurs_log[3,],
  OGK         = erreurs_log[length(methodes),],
  Offline     = erreurs_log[4,],
  Sample_cov  = erreurs_log[1,],
  Online      = erreurs_log[6,]
)

# -----------------------------
# couleurs
# -----------------------------
cols <- c(
  Streaming  = "red",
  Online     = "blue",
  Sample_cov = "darkgreen",
  Offline    = "orange",
  OGK        = "brown",
  MCD        = "black"
)

# -----------------------------
# ticks (valeurs originales -> transformées)
# -----------------------------
ticks_raw <- c(0, 1e-2, 1e-1)
ticks <- log1p(ticks_raw)

# -----------------------------
# boxplot
# -----------------------------
boxplot(
  data_plot,
  col = cols[names(data_plot)],
  main = "",
  ylab = "",
  names = c(
    "Streaming",
    "MCD",
    "OGK",
    "Offline",
    "Sample covariance + correction",
    "Online"
  ),
  ylim = c(0, max(unlist(erreurs_log), ticks)),
  yaxt = "n"
)

# -----------------------------
# axe Y propre (lecture en scale originale)
# -----------------------------
axis(
  side = 2,
  at = ticks,
  labels = c("0", "1e-2", "1e-1"),
  las = 1,
  cex.axis = 0.85
)

dev.off()


####################AUC ARI#######################
############################# AUC + ARI #################################
############################# PARAM #################################

methodes <- c(
  "Oracle",
  "SampleNaiveQuantonlinecorr",
  "SampleNaivewithoutQuantonlinecorr",
  "MCD",
  "Offline-without_onlinequantile",
  "Offline-without_onlinequantile",
  "OnlineUsQuantonlinecorr",
  "OnlineUsWithoutQuantonlinecorr",
  "StreamingUsonlineQuantcorr",
  "StreamingUs_without_onlineQuantcorr",
  "OGK"
)

n <- 28

############################# STOCKAGE #################################

ariPlot <- matrix(0, nrow = n, ncol = length(methodes))
aucPlot <- matrix(0, nrow = n, ncol = length(methodes))
temps   <- matrix(0, nrow = n, ncol = length(methodes))

for (j in 1:n) {
  
  setwd(crit_SMD)
  
  for (s in seq_along(methodes)) {
    
    methode <- methodes[s]
    
    critFile <- paste0("Crit-", methode, "-machine-", j, ".RData")
    
    load(critFile)
    
    ariPlot[j, s] <- crit$ARI
    aucPlot[j, s] <- crit$AUC
  }
  
  setwd(res_SMD)
  
  for (s in seq_along(methodes)) {
    
    methode <- methodes[s]
    
    fitFile <- paste0("Fit-", methode, "-machine-", j, ".RData")
    
    load(fitFile)
    
    temps[j, s] <- resultats$temps[3]
  }
}

############################# LABELS #################################

names_plot <- c(
  "Oracle",
  "Sple QC",
  "Sple no QC",
  "MCD",
  "Offl no QC",
  "Offl QC",
  "Onl QC",
  "Onl no QC",
  "Strm QC",
  "Strm no QC",
  "OGK"
)

############################# COULEURS #################################

cols <- c(
  "purple",
  "darkgreen",
  "darkgreen",
  "black",
  "orange",
  "orange",
  "blue",
  "blue",
  "red",
  "red",
  "brown"
)

############################# TICKS LOG (10^x style) ##################

tick_vals <- c(0, 1e-3, 1e-2, 1e-1, 1)

tick_pos <- log1p(tick_vals)

tick_lab <- c("0", "1e-3", "1e-2", "1e-1", "1")

############################# AUC #################################

setwd(fig_SMD)

pdf("boxplotAUC_SMD.pdf", width = 18, height = 6)

boxplot(
  as.data.frame(aucPlot),
  names = names_plot,
  las = 2,
  col = cols,
  main = "AUC",
  ylab = "AUC"
)

dev.off()

############################# ARI #################################

pdf("boxplotARI_SMD.pdf", width = 18, height = 6)

boxplot(
  as.data.frame(ariPlot),
  names = names_plot,
  las = 2,
  col = cols,
  main = "ARI",
  ylab = "ARI"
)

dev.off()
######################### TEMPS ##################################

setwd(fig_SMD)

pdf("boxplotTemps_SMD.pdf", width = 18, height = 6)

# enlever Oracle
temps_no_oracle <- temps[, -1]

# log transformation
temps_log <- log1p(temps_no_oracle)

# ticks en échelle originale -> log1p
tick_vals <- c(0, 1e-3, 1e-2, 1e-1, 1, 10)
tick_pos <- log1p(tick_vals)
tick_lab <- c("0", "1e-3", "1e-2", "1e-1", "1", "10")

boxplot(
  as.data.frame(temps_log),
  names = names_plot[-1],
  las = 2,
  col = cols[-1],
  main = "Time",
  ylab = "",
  yaxt = "n"
)

axis(
  side = 2,
  at = tick_pos,
  labels = tick_lab
)

dev.off()

###############Trajectoires##################

methodes_online = c("Oracle","SampleNaiveQuantonlinecorr","SampleNaivewithoutQuantonlinecorr","OnlineUsQuantonlinecorr","OnlineUsWithoutQuantonlinecorr","StreamingUsonlineQuantcorr","StreamingUs_without_onlineQuantcorr")



####Calcul des trajectoires pour les 3 méthodes online et l'oracle######
for(j in 1:28){
  
  
  setwd(smd_data_dir)
  
  
  data_smd_mach = paste0("data_machine-",j,".RData")
  
  
  load(data_smd_mach)
  outlmach = matrix(0,nrow = nrow(Z),ncol = length(methodes_online))
  
  setwd(res_SMD)
  
    for(s in seq_along(methodes_online)){
  
  methode = methodes_online[s]
  
  fitFile <- paste0('Fit-',methode,"-machine-",j,".RData")
  
  load(fitFile)
  
  outlmach[,s] = resultats$outliers_labels   
    }
  
  rates_samplecov_quantcorr = compute_rates(outlmach[,1], labels)
  
  rates_samplecov_without_quantcorr = compute_rates(outlmach[,2], labels)
  
  rates_online_with_quantcorr = compute_rates(outlmach[,3], labels)
  
  rates_online_without_quantcorr = compute_rates(outlmach[,4], labels)
  
  rates_Strm_with_quantcorr = compute_rates(outlmach[,5], labels)
  
  rates_Strm_without_quantcorr = compute_rates(outlmach[,6], labels)
  
  rates_oracle <- compute_rates(outlmach[,7], labels)
  
  ##################################Trajectoires################################"
  setwd(fig_SMD)
  nom_fichier = paste0("trajectories_SMD_mach-",j,".pdf")
  pdf(nom_fichier, width = 14, height = 10)
  
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  
  x_vals = 1:length(rates_Strm_with_quantcorr$FN_rate)
  
  
  # =====================================================
  # 1. LABELS
  # =====================================================
  
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  
  x_vals = 1:nrow(Z)
  
  
  # =====================================================
  # 1. LABELS
  # =====================================================
  plot(x_vals, labels,
       type = "n",
       ylim = c(-0.05,1.05),
       yaxt = "n",
       xaxt = "n",
       xlab = "", ylab = "",
       main = "Ground truth")
  
  points(x_vals[labels==0], labels[labels==0],
         pch = 16, cex = .4, col = "blue")
  
  points(x_vals[labels==1], labels[labels==1],
         pch = 16, cex = .6, col = "red")
  
  axis(2, at = c(0, 1), las = 1, cex.axis = 1.8)
  axis(1, at = seq(1000, max(x_vals), by = 1000),
       las = 1, cex.axis = 1.8)
  box()
  
  # =====================================================
  # 2. FALSE NEGATIVE RATE
  # =====================================================
  
  plot(x_vals, rates_Strm_with_quantcorr$FN_rate*100,
       type = "l", lwd = 4, col = "red",
       xlab = "", ylab = "",
       yaxt = "n", xaxt = "n",
       main = "False Negative Rate",
       ylim = c(0,100)
  )
  
  lines(x_vals, rates_Strm_without_quantcorr$FN_rate*100,
        lty = "longdash", col = "red", lwd = 3)
  
  lines(x_vals, rates_online_with_quantcorr$FN_rate*100,
        lty = "longdash", col = "blue", lwd = 3)
  
  lines(x_vals, rates_online_without_quantcorr$FN_rate*100,
        lty = "twodash", col = "blue", lwd = 3)
  
  lines(x_vals, rates_samplecov_quantcorr$FN_rate*100,
        lty = "dotted", col = "darkgreen", lwd = 3)
  
  lines(x_vals, rates_samplecov_without_quantcorr$FN_rate*100,
        lty = "dotdash", col = "darkgreen", lwd = 3)
  
  lines(x_vals, rates_oracle$FN_rate*100,
        lty = "dashed", col = "purple", lwd = 3)
  
  axis(2, las = 1, cex.axis = 1.8)
  axis(1, at = seq(1000, max(x_vals), by = 1000),
       las = 1, cex.axis = 1.8)
  box()
  
  # =====================================================
  # 3. FALSE POSITIVE RATE
  # =====================================================
  
  
  plot(x_vals, rates_Strm_with_quantcorr$FP_rate*100,
       type = "l", lwd = 3, col = "red",
       xlab = "", ylab = "",
       yaxt = "n", xaxt = "n",
       main = "False Positive Rate",
       ylim = c(0,20)
  )
  
  lines(x_vals, rates_Strm_without_quantcorr$FP_rate*100,
        lty = "longdash", col = "red", lwd = 3)
  
  lines(x_vals, rates_online_with_quantcorr$FP_rate*100,
        lty = "longdash", col = "blue", lwd = 3)
  
  lines(x_vals, rates_online_without_quantcorr$FP_rate*100,
        lty = "dotdash", col = "blue", lwd = 3)
  
  lines(x_vals, rates_samplecov_quantcorr$FP_rate*100,
        lty = "dotted", col = "darkgreen", lwd = 3)
  
  lines(x_vals, rates_samplecov_without_quantcorr$FP_rate*100,
        lty = "dotdash", col = "darkgreen", lwd = 3)
  
  lines(x_vals, rates_oracle$FP_rate*100,
        lty = "dashed", col = "purple", lwd = 3)
  
  axis(2, las = 1, cex.axis = 1.8)
  axis(1, at = seq(1000, max(x_vals), by = 1000),
       las = 1, cex.axis = 1.8)
  box()
  
  
  dev.off()  
}


