smd_data_dir = "~/smd_data_dir"
res_SMD = "~/res_SMD"
crit_SMD = "~/criteres_SMD"

########################Boxplot Sigma errors##############################

erreursSigmasSMD = array(0, dim = c(6,length(nb_machines)))  


for(j in 1:nb_machines){

setwd(crit_SMD)
  
critFile <- paste0('CritSampleNaiveQuantonlinecorr-',"-machine-",j,"-cutoff",.95,"Ninit-",1e3,"-cm",1,".RData")

load(critFile)

erreursSigmaSMD[1,j] =   FrobeniusNormError

critFile = paste0('CritMCD-machine-',j,".RData")

load(critFile)

erreursSigmaSMD[2,j] =   FrobeniusNormError

critFile = paste0('CritOGK-machine-',j,".RData")

load(critFile)

erreursSigmaSMD[3,j] =   FrobeniusNormError

critFile = paste0('CritOffline-withonlinequantile-machine-', j,".RData")

load(critFile)

erreursSigmaSMD[4,j] =   FrobeniusNormError

critFile = paste0('CritStreamingUsonlineQuantcorr-machine-', j,"-cutoff",.95,"Ninit-",1e3,"-batch",batchStrm,"-cm",1,".RData")

load(critFile)

erreursSigmaSMD[5,j] =   FrobeniusNormError

critFile = paste0('CritOnlineUswithoutQuantonlinecorr-machine-',j,"-cutoff",.95,"Ninit-",1e3,"-batch",1,"-cm",1,".RData")

load(critFile)

erreursSigmaSMD[6,j] =   FrobeniusNormError

}


file <- paste0("boxploterreursSigma_SMD-maj-",j ,".pdf")

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


###############Trajectoires##################

####Calcul des trajectoires pour les 3 méthodes online et l'oracle######
for(j in 1:nbmachines){
  
  
  setwd(smd_data_dir)
  
  
  data_smd_mach = paste0("data_machine-",j,".RData")
  
  
  load(data_smd_mach)
  

    
  
  setwd(res_SMD)
  
  fitFile <- paste0('FitSampleNaiveQuantonlinecorr-',"-machine-",j,"-cutoff",.95,"Ninit-",1e3,"-cm",1,".RData")
  
  load(fitFile)
  
  outlmach[,1] = resSamplecov$outliers_labels   

  rates_samplecov_quantcorr = compute_rates(outlmach[,1], labels)
  
  fitFile <- paste0('FitSampleNaivewithoutQuantonlinecorr-',"-machine-",j,"-cutoff",.95,"Ninit-",1e3,"-cm",1,".RData")
  
  load(fitFile)
  
  outlmach[,2] = resSamplecov$outliers_labels   

  rates_samplecov_without_quantcorr = compute_rates(outlmach[,2], labels)
  
  fitFile <- paste0('FitOnlineUsQuantonlinecorr-machine-',j,"-cutoff",.95,"Ninit-",1e3,"-batch",1,"-cm",1,".RData")
  
  load(fitFile)
  
  outlmach[,3] = resUsOnline$outlier_labels  
  
  rates_online_with_quantcorr = compute_rates(outlmach[,3], labels)
  
  fitFile <- paste0('FitOnlineUsWithoutQuantonlinecorr-machine-',j,"-cutoff",.95,"Ninit-",1e3,"-batch",1,"-cm",1,".RData")
  
  load(fitFile)
  
  outlmach[,4] = resUsOnline$outlier_labels
  
  rates_online_without_quantcorr = compute_rates(outlmach[,4], labels)
  
  fitFile <- paste0('FitStreamingUsonlineQuantcorr-machine-', j,"-cutoff",.95,"Ninit-",1e3,"-batch",batchStrm,"-cm",1,".RData")
  
  load(fitFile)
  
  outlmach[,5] = resStrm$outlier_labels
  
  rates_Strm_with_quantcorr = compute_rates(outlmach[,5], labels)
  
  fitFile <- paste0('FitStreamingUs_without_onlineQuantcorr-machine-', j,"-cutoff",.95,"Ninit-",1e3,"-batch",batchStrm,"-cm",1,".RData")
  
  load(fitFile)
  
  outlmach[,6] = resStrm$outlier_labels
  
  rates_Strm_without_quantcorr = compute_rates(outlmach[,6], labels)
  
  
  fitFile <- paste0('FitOracle-',"-machine-",j,".RData")
  
  load(fitFile)
  
  outlmach[,7] = outloracle
  
  rates_oracle <- compute_rates(outlmach[,7], labels)
  
  ##################################Trajectoires################################"
  setwd("~/figures")
  nom_fichier = paste0("trajectories_SMD_mach-",j,"-Ninit1e3cm1",".pdf")
  pdf(nom_fichier, width = 14, height = 10)
  
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  
  x_vals = 1:length(rates_Strm_with_quantcorr$FN_rate)
  
  
  # =====================================================
  # 1. LABELS
  # =====================================================
  
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  
  x_vals = 1:length(rates_strm_corr$FN_rate)
  
  
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
       main = "False Negative Rate",
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


