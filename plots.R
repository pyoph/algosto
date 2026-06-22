#############################Final Frobenius norm error, false positives, false negatives V2##############################

methodes = c("SampleNaiveQuantonlinecorr","SampleNaivewithoutonlinequantilecorr","OnlineUsQuantonlinecorr","OnlineUswithoutQuantonlinecorr","StreamingUsonlineQuantcorr","StreamingUswithoutQuantonlinecorr","OfflinewithQuantcorr","OfflineUswithoutQuantcorr","OGK","MCD","Oracle")

for (sc in scenarios){
  k = sc$k
  l = sc$l
  rho1 = sc$rho1
  
  erreursSigmaPlot = array(0,dim = c(length(rList[1:9]),6))
  faux_positifsPlot= array(0,dim = c(length(rList[1:9]),11))
  faux_negatifsPlot= array(0,dim = c(length(rList[1:9]),11))
  ariPlot = array(0,dim = c(length(rList[1:9]),11))
  aucPlot = array(0,dim = c(length(rList[1:9]),11))
  
  for (m in seq_along(rList[1:9])){
    
    ###############################Sample covonline#################################
    
    for(j in seq_along(methodes)){
      
      r = rList[m]
      
      methode = methodes[j]
      
      setwd(criteres)
      
      critFile <- paste0(
        'Crit-',methode, "-d",d,
        '-n', n,
        '-k', k,
        '-l', l,
        '-rho', rho1,
        '-r', r,
        '-sim', sim,
        ".RData"
      )
      load(critFile)
      
      erreursSigmaPlot[m,j] = crit$erreurFrob
      faux_negatifsPlot[m,j] = crit$FN
      faux_positifsPlot[m,j] = crit$FP
      ariPlot[m,j] = crit$ARI
      aucPlot[m,j] = crit$AUC
      
    }
  
    }  
}



for(sc in scenarios){
  k = sc$k
  l = sc$l
  rho1 = sc$rho1
setwd("~/figuresRobustVariance/fig/fig")

file <- paste0("scen-k", k, "-l", l, "-rho1", rho1,".pdf")

pdf(file, width = 25, height = 4)
par(mfrow = c(1,5), mar = c(4,4,2,1))
############################################################
################ ERREUR DE FROBENIUS #######################
############################################################

plot(rList[1:9], erreursSigmaPlot[1:9,3],
     type="l", lwd=4, col="red",
     log="y",
     ylim=c(1e-1,1e2),
     xlab="", ylab="",
     xaxt="n", yaxt="n")

lines(rList[1:9], erreursSigmaPlot[1:9,2],
      lwd=4, col="blue", lty="dashed")

lines(rList[1:9], erreursSigmaPlot[1:9,1],
      lwd=4, col="darkgreen", lty="dotted")

lines(rList[1:9], erreursSigmaPlot[1:9,4],
      lwd=4, col="orange", lty="dotted")

lines(rList[1:9], erreursSigmaPlot[1:9,5],
      lwd=4, col="black", lty="dotted")

lines(rList[1:9], erreursSigmaPlot[1:9,6],
      lwd=4, col="brown", lty="dotted")

axis(1, at=rList[1:9], cex.axis = 1.8,las=1)
axis(2,cex.axis = 1.8,
     at=10^seq(-1,2),
     labels=parse(text=paste0("10^", -1:2)),
     las=1)
box()

############################################################
################### FAUX NEGATIFS ##########################
############################################################
pseudo_log <- function(x) log10(1 + x)

plot(rList[2:9],
     pseudo_log(faux_negatifsPlot[2:9,5] / ((rList[2:9]/100)*n) * 100),
     type = "l", lwd = 4, col = "red",
     ylim = range(pseudo_log(c(0, 100))),
     xlab = "", ylab = "",
     xaxt = "n", yaxt = "n")

lines(rList[2:9],
      pseudo_log(faux_negatifsPlot[2:9,6] / ((rList[2:9]/100)*n) * 100),
      lwd = 4, col = "red", lty = "longdash")

lines(rList[2:9],
      pseudo_log(faux_negatifsPlot[2:9,1] / ((rList[2:9]/100)*n) * 100),
      lwd = 4, col = "darkgreen", lty = "dotted")

lines(rList[2:9],
      pseudo_log(faux_negatifsPlot[2:9,2] / ((rList[2:9]/100)*n) * 100),
      lwd = 4, col = "darkgreen", lty = "dotdash")

lines(rList[2:9],
      pseudo_log(faux_negatifsPlot[2:9,3] / ((rList[2:9]/100)*n) * 100),
      lwd = 4, col = "blue", lty = "longdash")

lines(rList[2:9],
      pseudo_log(faux_negatifsPlot[2:9,4] / ((rList[2:9]/100)*n) * 100),
      lwd = 4, col = "blue", lty = "twodash")

axis(1, at = rList[2:9], las = 1, cex.axis = 1.8)

axis(2,
     at = log10(1 + c(0, 1, 10, 100)),
     labels = c("0", "1", "10", "100"),
     las = 1,
     cex.axis = 1.8)

box()
############################################################
################### FAUX POSITIFS ##########################
############################################################

plot(rList[1:9],
     faux_positifsPlot[,5] /
                  ((1 - rList[1:9]/100)*n)*100,
     type="l", lwd=4, col="red",
     ylim=c(0,20),
     xlab="", ylab="",
     xaxt="n", yaxt="n")

lines(rList[1:9],
      faux_positifsPlot[,6] /
              ((1 - rList[1:9]/100)*n)*100,
      lwd=4, col="red", lty="longdash")


lines(rList[1:9],
      faux_positifsPlot[1:9,1] /
                   ((1 - rList[1:9]/100)*n)*100,
      lwd=4, col="darkgreen", lty="dotted")



lines(rList[1:9],
      faux_positifsPlot[1:9,2] /
                   ((1 - rList[1:9]/100)*n)*100,
      lwd=4, col="darkgreen", lty="dotdash")


lines(rList[1:9],
      faux_positifsPlot[1:9,3] /
                   ((1 - rList[1:9]/100)*n)*100,
      lwd=4, col="blue", lty="longdash")

lines(rList[1:9],
      faux_positifsPlot[1:9,4] /
                   ((1 - rList[1:9]/100)*n)*100,
      lwd=4, col="blue", lty="twodash")

lines(rList[1:9],
      faux_positifsPlot[1:9,11] /
                   ((1 - rList[1:9]/100)*n)*100,
      lwd=4, col="purple4", lty="longdash")

lines(rList[1:9],
      faux_positifsPlot[1:9,7] /
                   ((1 - rList[1:9]/100)*n)*100,
      lwd=4, col="orange", lty="longdash")

lines(rList[1:9],
      pseudo_log(faux_positifsPlot[1:9,8] /
                  (1 - rList[1:9]/100)*n)*100,
      lwd=4, col="orange", lty="twodash")


lines(rList[1:9],
      faux_positifsPlot[1:9,9] /
                   ((1 - rList[1:9]/100)*n)*100,
      lwd=4, col="brown", lty="dotted")

lines(rList[1:9],
      faux_positifsPlot[1:9,10] /
                   ((1 - rList[1:9]/100)*n)*100,
      lwd=4, col="black", lty="twodash")



axis(1, at = rList[1:9], las = 1, cex.axis = 1.8)
axis(2,
     at = c(0, 5, 10, 15, 20),
     labels = c("0", "5", "10", "15", "20"),
     las = 1,
     cex.axis = 1.8)
box()

############################################################
######################## ARI ###############################
############################################################

plot(rList[2:9],
     ariPlot[2:9,5],
     type="l", lwd=3, col="red",
     ylim=c(-1,1),
     xlab="", ylab="",
     xaxt="n", yaxt="n")
lines(rList[2:9],
      ariPlot[2:9,6],
      lwd=4, col="red", lty="longdash")


lines(rList[2:9],
      ariPlot[2:9,1],
      lwd=4, col="darkgreen", lty="dotted")



lines(rList[2:9],
      ariPlot[2:9,2],
      lwd=4, col="darkgreen", lty="dotdash")


lines(rList[2:9],
      ariPlot[2:9,3],
      lwd=4, col="blue", lty="longdash")

lines(rList[2:9],
      ariPlot[2:9,4],
      lwd=4, col="blue", lty="twodash")

lines(rList[2:9],
      ariPlot[2:9,11],
      lwd=4, col="purple4", lty="longdash")

lines(rList[2:9],
      ariPlot[2:9,7],
      lwd=4, col="orange", lty="longdash")

lines(rList[2:9],
      ariPlot[2:9,8],
      lwd=4, col="orange", lty="twodash")


lines(rList[2:9],
      ariPlot[2:9,9],
      lwd=4, col="brown", lty="dotted")

lines(rList[2:9],
      ariPlot[2:9,10],
      lwd=4, col="black", lty="twodash")

axis(1, at = rList[2:9], las = 1, cex.axis = 1.8)
axis(2, at = seq(-1,1,0.5),
     labels = seq(-1,1,0.5),
     las = 1, cex.axis = 1.8)
box()
############################################################
######################## AUC ###############################
############################################################



plot(rList[2:9],
     aucPlot[2:9,5],
     type="l", lwd=3, col="red",
     ylim=c(0,1),
     xlab="", ylab="",
     xaxt="n", yaxt="n")
lines(rList[2:9],
      aucPlot[2:9,6],
      lwd=4, col="red", lty="longdash")


lines(rList[2:9],
      aucPlot[2:9,1],
      lwd=4, col="darkgreen", lty="dotted")



lines(rList[2:9],
      aucPlot[2:9,2],
      lwd=4, col="darkgreen", lty="dotdash")


lines(rList[2:9],
      aucPlot[2:9,3],
      lwd=4, col="blue", lty="longdash")

lines(rList[2:9],
      aucPlot[2:9,4],
      lwd=4, col="blue", lty="twodash")

lines(rList[2:9],
      aucPlot[2:9,11],
      lwd=4, col="purple4", lty="longdash")

lines(rList[2:9],
      aucPlot[2:9,7],
      lwd=4, col="orange", lty="longdash")

lines(rList[2:9],
      aucPlot[2:9,8],
      lwd=4, col="orange", lty="twodash")


lines(rList[2:9],
      aucPlot[2:9,9],
      lwd=4, col="brown", lty="dotted")

lines(rList[2:9],
      aucPlot[2:9,10],
      lwd=4, col="black", lty="twodash")

axis(1, at = rList[2:9], cex.axis = 1.8, las = 1)
axis(2, at = seq(0,1,0.2),
     labels = seq(0,1,0.2),
     cex.axis = 1.8, las = 1)

box()

dev.off()


}
#######################Trajectoires##########################

for(sc in scenarios){
  
  k = sc$k
  l = sc$l
  rho1 = sc$rho1
  

  
  for (r in rList[1:9]){
  
  outlabTraj= array(0, dim = c(n,7))
  
  setwd(SimDir)  
  
  dataFile <- paste0('SimData-d', d, '-n', n, '-k', k, '-l', l, '-rho', rho1,'-r',r ,"-sim",sim,".RData")
  
  print(dataFile)
  
  
  load(dataFile)

  labels = data$labelsVrais
  
  setwd(resAlgo)
  
  fitFile <- paste0('FitSampleNaiveQuantonlinecorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,"-cutoff",.95,"Ninit-",Ninit,"-cm",cm,'-sim', sim,".RData")
  
  load(fitFile)
  
  outlabTraj[,1] = resultats$outliers_labels
  
  
  fitFile <- paste0('FitSampleNaivewithoutonlinequantilecorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,"-Ninit",Ninit,'-sim', sim,".RData")
  
  load(fitFile)
  
  outlabTraj[,2] = resultats$outliers_labels
  
  
  fitFile <- paste0('FitOnlineUsQuantonlinecorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,"-cutoff",.95,"Ninit-",Ninit,"-batch",1,"-cm",cm,'-sim', sim,".RData")
  
  
  load(fitFile)
  
  outlabTraj[,3] = resultats$outliers_labels

  fitFile <- paste0('FitOnlineUswithoutQuantonlinecorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,"-cutoff",.95,"Ninit-",1e2,"-batch",1,'-sim', sim,".RData")
  
  load(fitFile)
  
  outlabTraj[,4] = resultats$outliers_labels
  
  fitFile <- paste0('FitStreamingUsonlineQuantcorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,"-cutoff",.95,"Ninit-",Ninit,"-batch",batch,"-cm",cm,'-sim', sim,".RData")
  
  
  load(fitFile)
  
  outlabTraj[,5] = resultats$outliers_labels
  
  
  fitFile <- paste0('FitStreamingUswithoutQuantonlinecorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,"Ninit-",1e2,"-batch",batch,'-sim', sim,".RData")
  
  
  load(fitFile)
  
  outlabTraj[,6] = resultats$outliers_labels
  
  
  
  fitFile <- paste0('FitStreamingUswithoutQuantonlinecorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,"Ninit-",1e2,"-batch",batch,'-sim', sim,".RData")
  
  load(fitFile)
  
  outlabTraj[,6] = resultats$outliers_labels
  
  
  fitFile <- paste0('FitOracle-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,"-cutoff",.95,'-sim', sim,".RData")
  
  load(fitFile)
  
  outlabTraj[,7] = resultats$outliers_labels
  
  
  #rates_samplecov = compute_rates(outlmach[,1], labels)
  
  rates_samplecov_wcorr = compute_rates(outlabTraj[,1], labels)
  
  rates_samplecov_wocorr = compute_rates(outlabTraj[,2], labels)
  
  rates_online_corr = compute_rates(outlabTraj[,3], labels)
  
  rates_online_wocorr = compute_rates(outlabTraj[,4], labels)
  
  rates_strm_corr = compute_rates(outlabTraj[,5], labels)
  
  rates_strm_wocorr = compute_rates(outlabTraj[,6], labels)
  
  rates_oracle = compute_rates(outlabTraj[,7], labels)
  
  setwd("~/figures")
  
  nom_fichier = paste0("trajectories_k-",k,"-l",l,"-rho1",rho1,"-r",r,"-sim",sim,".pdf")
  pdf(nom_fichier, width = 14, height = 10)
  
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
  
  plot(x_vals, rates_strm_corr$FN_rate*100,
       type = "l", lwd = 4, col = "red",
       xlab = "", ylab = "",
       yaxt = "n", xaxt = "n",
       main = "False Negative Rate",
       ylim = c(0,100)
  )
  
  lines(x_vals, rates_strm_wocorr$FN_rate*100,
        lty = "longdash", col = "red", lwd = 3)
  
  lines(x_vals, rates_online_corr$FN_rate*100,
        lty = "longdash", col = "blue", lwd = 3)
  
  lines(x_vals, rates_online_wocorr$FN_rate*100,
        lty = "twodash", col = "blue", lwd = 3)
  
  lines(x_vals, rates_samplecov_wcorr$FN_rate*100,
        lty = "dotted", col = "darkgreen", lwd = 3)
  
  lines(x_vals, rates_samplecov_wocorr$FN_rate*100,
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
  
  
  plot(x_vals, rates_strm_corr$FP_rate*100,
       type = "l", lwd = 3, col = "red",
       xlab = "", ylab = "",
       yaxt = "n", xaxt = "n",
       main = "False Positive Rate",
       ylim = c(0,20)
  )
  
  lines(x_vals, rates_strm_wocorr$FP_rate*100,
        lty = "longdash", col = "red", lwd = 3)
  
  lines(x_vals, rates_online_corr$FP_rate*100,
        lty = "longdash", col = "blue", lwd = 3)
  
  lines(x_vals, rates_online_wocorr$FP_rate*100,
        lty = "dotdash", col = "blue", lwd = 3)
  
  lines(x_vals, rates_samplecov_wcorr$FP_rate*100,
        lty = "dotted", col = "darkgreen", lwd = 3)
  
  lines(x_vals, rates_samplecov_wocorr$FP_rate*100,
        lty = "dotdash", col = "darkgreen", lwd = 3)
  
  lines(x_vals, rates_oracle$FP_rate*100,
        lty = "dashed", col = "purple", lwd = 3)
  
  axis(2, las = 1, cex.axis = 1.8)
  axis(1, at = seq(1000, max(x_vals), by = 1000),
       las = 1, cex.axis = 1.8)
  box()
  
  
  dev.off()

  
  }
  
  }



