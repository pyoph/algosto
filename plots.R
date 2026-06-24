###################################Boxplots inliers outliers##########

r = 30

for (sc in scenarios){
  k = sc$k
  l = sc$l
  rho1 = sc$rho1
  
  setwd(SimDir)
  
  dataFile <- paste0('SimData-d', d, '-n', n, '-k', k, '-l', l, '-rho', rho1,'-r',r ,"-sim",sim,".RData")
  
  print(dataFile)
  
  load(dataFile)
  
  setwd(resAlgo)
  
  fitFile = paste0('Fit-Oracle-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
  
  load(fitFile)
  
  distances = resultats$distances
  
  setwd("~")
  
  file <- paste0("boxplot_inliers_outliers-k",k,"-l",l,"-rho",rho1,"-r",r,".pdf")
  
  pdf(file, width = 25, height = 4)
  
  
  boxplot(
    distances[data$labelsVrais == 0],
    distances[data$labelsVrais == 1],
    names = c("0", "1"),
    col = c("blue", "red")
  )  
  dev.off()
}





#############################Final Frobenius norm error, false positives, false negatives V2##############################

methodes = c("SampleNaiveQuantonlinecorr","SampleNaivewithoutonlinequantilecorr","OnlineUsQuantonlinecorr","OnlineUswithoutQuantonlinecorr","StreamingUsonlineQuantcorr","StreamingUswithoutQuantonlinecorr","OfflinewithQuantcorr","OfflineUswithoutQuantcorr","OGK","MCD","Oracle")

for (sc in scenarios){
  k = sc$k
  l = sc$l
  rho1 = sc$rho1
  
  erreursSigmaPlot = array(0,dim = c(length(rList[1:9]),11))
  faux_positifsPlot= array(0,dim = c(length(rList[1:9]),11))
  faux_negatifsPlot= array(0,dim = c(length(rList[1:9]),11))
  ariPlot = array(0,dim = c(length(rList[1:9]),11))
  aucPlot = array(0,dim = c(length(rList[1:9]),11))
  
  for (m in seq_along(rList[1:9])){
    
    r = rList[m]
    
    
    for(j in seq_along(methodes)){
      
      
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
    setwd("~/figures2")
    
    file <- paste0("scen-k", k, "-l", l, "-rho1", rho1,".pdf")
    
    pdf(file, width = 25, height = 4)
    par(mfrow = c(1,5), mar = c(4,4,2,1))
    ############################################################
    ################ ERREUR DE FROBENIUS #######################
    ############################################################
    
    plot(rList[1:9], erreursSigmaPlot[1:9,5],
         type="l", lwd=4, col="red",
         log="y",
         ylim=c(1e-1,1e2),
         xlab="", ylab="",
         xaxt="n", yaxt="n")
    
    lines(rList[1:9], erreursSigmaPlot[1:9,3],
          lwd=4, col="blue", lty="dashed")
    
    lines(rList[1:9], erreursSigmaPlot[1:9,1],
          lwd=4, col="darkgreen", lty="dotted")
    
    lines(rList[1:9], erreursSigmaPlot[1:9,7],
          lwd=4, col="orange", lty="dotted")
    
    lines(rList[1:9], erreursSigmaPlot[1:9,10],
          lwd=4, col="black", lty="dotted")
    
    lines(rList[1:9], erreursSigmaPlot[1:9,9],
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
          pseudo_log(faux_negatifsPlot[2:9,7] / ((rList[2:9]/100) * n) * 100),
          lwd=4, col="orange", lty="longdash")
    
    lines(rList[2:9],
          pseudo_log(faux_negatifsPlot[2:9,8] / ((rList[2:9]/100) * n) * 100),
          lwd=4, col="orange", lty="twodash")
    
    
    lines(rList[2:9],
          pseudo_log(faux_negatifsPlot[2:9,11] /
                       ((rList[2:9]/100) * n) * 100),
          lwd=4, col="purple", lty="longdash")
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
}

######Plots pour k qui varie

r = 30

k_values = seq(0,100,step = 1)

erreursSigmaPlot = array(0,dim = c(length(k_values),11))
faux_positifsPlot= array(0,dim = c(length(k_values),11))
faux_negatifsPlot= array(0,dim = c(length(k_values),11))
ariPlot = array(0,dim = c(length(k_values),11))
aucPlot = array(0,dim = c(length(k_values),11))


for (k in k_values){
  
  l = 1 
  
  rho1 = 0.3
  
  
  
  for(j in seq_along(methodes)){
    
    
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
    
  
    erreursSigmaPlot[k,j] = crit$erreurFrob
    faux_negatifsPlot[k,j] = crit$FN
    faux_positifsPlot[k,j] = crit$FP
    ariPlot[k,j] = crit$ARI
    aucPlot[k,j] = crit$AUC
    
  }}  
  setwd("~/figures")
  
  k_values = k_values[1:1e2]
  erreursSigmaPlot = erreursSigmaPlot[1:1e2,]
  faux_negatifsPlot = faux_negatifsPlot[1:1e2,]
  faux_positifsPlot = faux_positifsPlot[1:1e2,]
  ariPlot = ariPlot[1:1e2,]
  aucPlot = aucPlot[1:1e2,]
  
  file <- paste0("scen-k",".pdf")
  
  pdf(file, width = 25, height = 4)
  par(mfrow = c(1,5), mar = c(4,4,2,1))
  ############################################################
  ################ ERREUR DE FROBENIUS #######################
  ############################################################
  
  plot(k_values, erreursSigmaPlot[,5],
       type="l", lwd=4, col="red",
       log="y",
       ylim=c(1e-1,1e2),
       xlab="", ylab="",
       xaxt="n", yaxt="n")
  
  lines(k_values, erreursSigmaPlot[,3],
        lwd=4, col="blue", lty="dashed")
  
  lines(k_values, erreursSigmaPlot[,1],
        lwd=4, col="darkgreen", lty="dotted")
  
  lines(k_values, erreursSigmaPlot[,7],
        lwd=4, col="orange", lty="dotted")
  
  lines(k_values, erreursSigmaPlot[,10],
        lwd=4, col="black", lty="dotted")
  
  lines(k_values, erreursSigmaPlot[,9],
        lwd=4, col="brown", lty="dotted")
  
  axis(1, at= k_values, cex.axis = 1.8,las=1)
  axis(2,cex.axis = 1.8,
       at=10^seq(-1,2),
       labels=parse(text=paste0("10^", -1:2)),
       las=1)
  box()
  
  pseudo_log <- function(x) log10(1 + x)
  
  plot(k_values,
       pseudo_log(faux_negatifsPlot[,5] / ((30/100)*n) * 100),
       type = "l", lwd = 4, col = "red",
       ylim = range(pseudo_log(c(0, 100))),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  
  lines(k_values,
        pseudo_log(faux_negatifsPlot[,6] / ((30/100)*n) * 100),
        lwd = 4, col = "red", lty = "longdash")
  
  lines(k_values,
        pseudo_log(faux_negatifsPlot[,1] / ((30/100)*n) * 100),
        lwd = 4, col = "darkgreen", lty = "dotted")
  
  lines(k_values,
        pseudo_log(faux_negatifsPlot[,2] / ((30/100)*n) * 100),
        lwd = 4, col = "darkgreen", lty = "dotdash")
  
  lines(k_values,
        pseudo_log(faux_negatifsPlot[,3] / ((30/100)*n) * 100),
        lwd = 4, col = "blue", lty = "longdash")
  
  lines(k_values,
        pseudo_log(faux_negatifsPlot[,4] / ((30/100)*n) * 100),
        lwd = 4, col = "blue", lty = "twodash")
  
  lines(k_values,
        pseudo_log(faux_negatifsPlot[,7] / ((30/100)*n) * 100),
        lwd = 4, col = "orange", lty = "longdash")
  
  lines(k_values,
        pseudo_log(faux_negatifsPlot[,8] / ((30/100)*n) * 100),
        lwd = 4, col = "orange", lty = "twodash")
  
  axis(1,
       at = k_values,
       labels = k_values,
       las = 1,
       cex.axis = 1.6)
  
  axis(2,
       at = log10(1 + c(0, 1, 10, 100)),
       labels = c("0", "1", "10", "100"),
       las = 1,
       cex.axis = 1.6)
  
  box()
  
  
  ##############Faux positifs
  
  plot(k_values,
       faux_positifsPlot[,5] /
         ((1 - 30/100)*n)*100,
       type="l", lwd=4, col="red",
       ylim=c(0,20),
       xlab="", ylab="",
       xaxt="n", yaxt="n")
  
  lines(k_values,
        faux_positifsPlot[,6] /
          ((1 - 30/100)*n)*100,
        lwd=4, col="red", lty="longdash")
  
  
  lines(k_values,
        faux_positifsPlot[,1] /
          ((1 - 30/100)*n)*100,
        lwd=4, col="darkgreen", lty="dotted")
  
  
  
  lines(k_values,
        faux_positifsPlot[,2] /
          ((1 - 30/100)*n)*100,
        lwd=4, col="darkgreen", lty="dotdash")
  
  
  lines(k_values,
        faux_positifsPlot[,3] /
          ((1 - 30/100)*n)*100,
        lwd=4, col="blue", lty="longdash")
  
  lines(k_values,
        faux_positifsPlot[,4] /
          ((1 - 30/100)*n)*100,
        lwd=4, col="blue", lty="twodash")
  
  lines(k_values,
        faux_positifsPlot[,11] /
          ((1 - 30/100)*n)*100,
        lwd=4, col="purple4", lty="longdash")
  
  lines(k_values,
        faux_positifsPlot[,7] /
          ((1 - 30/100)*n)*100,
        lwd=4, col="orange", lty="longdash")
  
  lines(k_values,
        faux_positifsPlot[,8] /
                     ((1 - 30/100)*n)*100,
        lwd=4, col="orange", lty="twodash")
  
  
  lines(k_values,
        faux_positifsPlot[,9] /
          ((1 - 30/100)*n)*100,
        lwd=4, col="brown", lty="dotted")
  
  lines(k_values,
        faux_positifsPlot[,10] /
          ((1 - 30/100)*n)*100,
        lwd=4, col="black", lty="twodash")
  
  
  
  axis(1, at = k_values, las = 1, cex.axis = 1.8)
  axis(2,
       at = c(0, 5, 10, 15, 20),
       labels = c("0", "5", "10", "15", "20"),
       las = 1,
       cex.axis = 1.8)
  box()
  
  ############################################################
  ######################## ARI ###############################
  ############################################################
  
  plot(k_values,
       ariPlot[,5],
       type="l", lwd=3, col="red",
       ylim=c(-1,1),
       xlab="", ylab="",
       xaxt="n", yaxt="n")
  lines(k_values,
        ariPlot[,6],
        lwd=4, col="red", lty="longdash")
  
  
  lines(k_values,
        ariPlot[,1],
        lwd=4, col="darkgreen", lty="dotted")
  
  
  
  lines(k_values,
        ariPlot[,2],
        lwd=4, col="darkgreen", lty="dotdash")
  
  
  lines(k_values,
        ariPlot[,3],
        lwd=4, col="blue", lty="longdash")
  
  lines(k_values,
        ariPlot[,4],
        lwd=4, col="blue", lty="twodash")
  
  lines(k_values,
        ariPlot[,11],
        lwd=4, col="purple4", lty="longdash")
  
  lines(k_values,
        ariPlot[,7],
        lwd=4, col="orange", lty="longdash")
  
  lines(k_values,
        ariPlot[,8],
        lwd=4, col="orange", lty="twodash")
  
  
  lines(k_values,
        ariPlot[,9],
        lwd=4, col="brown", lty="dotted")
  
  lines(k_values,
        ariPlot[,10],
        lwd=4, col="black", lty="twodash")
  
  axis(1, at = k_values, las = 1, cex.axis = 1.8)
  axis(2, at = seq(-1,1,0.5),
       labels = seq(-1,1,0.5),
       las = 1, cex.axis = 1.8)
  box()
  ############################################################
  ######################## AUC ###############################
  ############################################################
  
  
  
  plot(k_values,
       aucPlot[,5],
       type="l", lwd=3, col="red",
       ylim=c(0,1),
       xlab="", ylab="",
       xaxt="n", yaxt="n")
  lines(k_values,
        aucPlot[,6],
        lwd=4, col="red", lty="longdash")
  
  
  lines(k_values,
        aucPlot[,1],
        lwd=4, col="darkgreen", lty="dotted")
  
  
  
  lines(k_values,
        aucPlot[,2],
        lwd=4, col="darkgreen", lty="dotdash")
  
  
  lines(k_values,
        aucPlot[,3],
        lwd=4, col="blue", lty="longdash")
  
  lines(k_values,
        aucPlot[,4],
        lwd=4, col="blue", lty="twodash")
  
  lines(k_values,
        aucPlot[,11],
        lwd=4, col="purple4", lty="longdash")
  
  lines(k_values,
        aucPlot[,7],
        lwd=4, col="orange", lty="longdash")
  
  lines(k_values,
        aucPlot[,8],
        lwd=4, col="orange", lty="twodash")
  
  
  lines(k_values,
        aucPlot[,9],
        lwd=4, col="brown", lty="dotted")
  
  lines(k_values,
        aucPlot[,10],
        lwd=4, col="black", lty="twodash")
  
  axis(1, at = k_values, cex.axis = 1.8, las = 1)
  axis(2, at = seq(0,1,0.2),
       labels = seq(0,1,0.2),
       cex.axis = 1.8, las = 1)
  
  box()
  
  dev.off()
  
  
  ######Plots pour l < 1 qui varie
  
  r = 30
  
  l_inf1_values = seq(0.01,0.99,by = 0.01)
  
  erreursSigmaPlot = array(0,dim = c(length(l_inf1_values),11))
  faux_positifsPlot= array(0,dim = c(length(l_inf1_values),11))
  faux_negatifsPlot= array(0,dim = c(length(l_inf1_values),11))
  ariPlot = array(0,dim = c(length(l_inf1_values),11))
  aucPlot = array(0,dim = c(length(l_inf1_values),11))
  
  
  for (s in seq_along(l_inf1_values)){
    
   k = 0 
    
   l = l_inf1_values[s]
    rho1 = 0.3
    
    
    
    for(j in seq_along(methodes)){
      
      
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
      
      
      erreursSigmaPlot[s,j] = crit$erreurFrob
      faux_negatifsPlot[s,j] = crit$FN
      faux_positifsPlot[s,j] = crit$FP
      ariPlot[s,j] = crit$ARI
      aucPlot[s,j] = crit$AUC
      
    }}  
  setwd("~/figures")
  # 
  # erreursSigmaPlot = erreursSigmaPlot[1:1e2,]
  # faux_negatifsPlot = faux_negatifsPlot[1:1e2,]
  # faux_positifsPlot = faux_positifsPlot[1:1e2,]
  # ariPlot = ariPlot[1:1e2,]
  # aucPlot = aucPlot[1:1e2,]
  
  file <- paste0("scen-linf1",".pdf")
  
  pdf(file, width = 25, height = 4)
  par(mfrow = c(1,5), mar = c(4,4,2,1))
  ############################################################
  ################ ERREUR DE FROBENIUS #######################
  ############################################################
  
  plot(l_inf1_values, erreursSigmaPlot[,5],
       type="l", lwd=4, col="red",
       log="y",
       ylim=c(1e-1,1e2),
       xlab="", ylab="",
       xaxt="n", yaxt="n")
  
  lines(l_inf1_values, erreursSigmaPlot[,3],
        lwd=4, col="blue", lty="dashed")
  
  lines(l_inf1_values, erreursSigmaPlot[,1],
        lwd=4, col="darkgreen", lty="dotted")
  
  lines(l_inf1_values, erreursSigmaPlot[,7],
        lwd=4, col="orange", lty="dotted")
  
  lines(l_inf1_values, erreursSigmaPlot[,10],
        lwd=4, col="black", lty="dotted")
  
  lines(l_inf1_values, erreursSigmaPlot[,9],
        lwd=4, col="brown", lty="dotted")
  
  axis(1, at= l_inf1_values, cex.axis = 1.8,las=1)
  axis(2,cex.axis = 1.8,
       at=10^seq(-1,2),
       labels=parse(text=paste0("10^", -1:2)),
       las=1)
  box()
  
  ############################################################
  ################### FAUX NEGATIFS ##########################
  ############################################################
  pseudo_log <- function(x) log10(1 + x)
  
  plot(l_inf1_values,
       pseudo_log(faux_negatifsPlot[,5] / ((30/100)*n) * 100),
       type = "l", lwd = 4, col = "red",
       ylim = range(pseudo_log(c(0, 100))),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  
  lines(l_inf1_values,
        pseudo_log(faux_negatifsPlot[,6] / ((30/100)*n) * 100),
        lwd = 4, col = "red", lty = "longdash")
  
  
  
  
  
  lines(l_inf1_values,
        pseudo_log(faux_negatifsPlot[,1] / ((30/100)*n) * 100),
        lwd = 4, col = "darkgreen", lty = "dotted")
  
  lines(l_inf1_values,
        pseudo_log(faux_negatifsPlot[,2] / ((30/100)*n) * 100),
        lwd = 4, col = "darkgreen", lty = "dotdash")
  
  lines(l_inf1_values,
        pseudo_log(faux_negatifsPlot[,3] / ((30/100)*n) * 100),
        lwd = 4, col = "blue", lty = "longdash")
  
  lines(l_inf1_values,
        pseudo_log(faux_negatifsPlot[,4] / ((30/100)*n) * 100),
        lwd = 4, col = "blue", lty = "twodash")
  
  lines(l_inf1_values,
        pseudo_log(faux_negatifsPlot[,7] / ((30/100)*n) * 100),
        lwd = 4, col = "orange", lty = "longdash")
  
  lines(l_inf1_values,
        pseudo_log(faux_negatifsPlot[,8] / ((30/100)*n) * 100),
        lwd = 4, col = "orange", lty = "twodash")
  
  lines(l_inf1_values,
        pseudo_log(faux_negatifsPlot[,11] / ((30/100) * n) * 100),
        lwd = 4, col = "purple", lty = "longdash")
  axis(1, at = l_inf1_values, las = 1, cex.axis = 1.8)
  
  axis(2,
       at = log10(1 + c(0, 1, 10, 100)),
       labels = c("0", "1", "10", "100"),
       las = 1,
       cex.axis = 1.8)
  
  box()
  ############################################################
  ################### FAUX POSITIFS ##########################
  ############################################################
  
  plot(l_inf1_values,
       faux_positifsPlot[,5] /
         ((1 - 30/100)*n)*100,
       type="l", lwd=4, col="red",
       ylim=c(0,20),
       xlab="", ylab="",
       xaxt="n", yaxt="n")
  
  lines(l_inf1_values,
        faux_positifsPlot[,6] /
          ((1 - 30/100)*n)*100,
        lwd=4, col="red", lty="longdash")
  
  
  lines(l_inf1_values,
        faux_positifsPlot[,1] /
          ((1 - 30/100)*n)*100,
        lwd=4, col="darkgreen", lty="dotted")
  
  
  
  lines(l_inf1_values,
        faux_positifsPlot[,2] /
          ((1 - 30/100)*n)*100,
        lwd=4, col="darkgreen", lty="dotdash")
  
  
  lines(l_inf1_values,
        faux_positifsPlot[,3] /
          ((1 - 30/100)*n)*100,
        lwd=4, col="blue", lty="longdash")
  
  lines(l_inf1_values,
        faux_positifsPlot[,4] /
          ((1 - 30/100)*n)*100,
        lwd=4, col="blue", lty="twodash")
  
  lines(l_inf1_values,
        faux_positifsPlot[,11] /
          ((1 - 30/100)*n)*100,
        lwd=4, col="purple4", lty="longdash")
  
  lines(l_inf1_values,
        faux_positifsPlot[,7] /
          ((1 - 30/100)*n)*100,
        lwd=4, col="orange", lty="longdash")
  
  lines(l_inf1_values,
        pseudo_log(faux_positifsPlot[,8] /
                     (1 - 30/100)*n)*100,
        lwd=4, col="orange", lty="twodash")
  
  
  lines(l_inf1_values,
        faux_positifsPlot[,9] /
          ((1 - 30/100)*n)*100,
        lwd=4, col="brown", lty="dotted")
  
  lines(l_inf1_values,
        faux_positifsPlot[,10] /
          ((1 - 30/100)*n)*100,
        lwd=4, col="black", lty="twodash")
  
  
  
  axis(1, at = l_inf1_values, las = 1, cex.axis = 1.8)
  axis(2,
       at = c(0, 5, 10, 15, 20),
       labels = c("0", "5", "10", "15", "20"),
       las = 1,
       cex.axis = 1.8)
  box()
  
  ############################################################
  ######################## ARI ###############################
  ############################################################
  
  plot(l_inf1_values,
       ariPlot[,5],
       type="l", lwd=3, col="red",
       ylim=c(-1,1),
       xlab="", ylab="",
       xaxt="n", yaxt="n")
  lines(l_inf1_values,
        ariPlot[,6],
        lwd=4, col="red", lty="longdash")
  
  
  lines(l_inf1_values,
        ariPlot[,1],
        lwd=4, col="darkgreen", lty="dotted")
  
  
  
  lines(l_inf1_values,
        ariPlot[,2],
        lwd=4, col="darkgreen", lty="dotdash")
  
  
  lines(l_inf1_values,
        ariPlot[,3],
        lwd=4, col="blue", lty="longdash")
  
  lines(l_inf1_values,
        ariPlot[,4],
        lwd=4, col="blue", lty="twodash")
  
  lines(l_inf1_values,
        ariPlot[,11],
        lwd=4, col="purple4", lty="longdash")
  
  lines(l_inf1_values,
        ariPlot[,7],
        lwd=4, col="orange", lty="longdash")
  
  lines(l_inf1_values,
        ariPlot[,8],
        lwd=4, col="orange", lty="twodash")
  
  
  lines(l_inf1_values,
        ariPlot[,9],
        lwd=4, col="brown", lty="dotted")
  
  lines(l_inf1_values,
        ariPlot[,10],
        lwd=4, col="black", lty="twodash")
  
  axis(1, at = l_inf1_values, las = 1, cex.axis = 1.8)
  axis(2, at = seq(-1,1,0.5),
       labels = seq(-1,1,0.5),
       las = 1, cex.axis = 1.8)
  box()
  ############################################################
  ######################## AUC ###############################
  ############################################################
  
  
  
  plot(l_inf1_values,
       aucPlot[,5],
       type="l", lwd=3, col="red",
       ylim=c(0,1),
       xlab="", ylab="",
       xaxt="n", yaxt="n")
  lines(l_inf1_values,
        aucPlot[,6],
        lwd=4, col="red", lty="longdash")
  
  
  lines(l_inf1_values,
        aucPlot[,1],
        lwd=4, col="darkgreen", lty="dotted")
  
  
  
  lines(l_inf1_values,
        aucPlot[,2],
        lwd=4, col="darkgreen", lty="dotdash")
  
  
  lines(l_inf1_values,
        aucPlot[,3],
        lwd=4, col="blue", lty="longdash")
  
  lines(l_inf1_values,
        aucPlot[,4],
        lwd=4, col="blue", lty="twodash")
  
  lines(l_inf1_values,
        aucPlot[,11],
        lwd=4, col="purple4", lty="longdash")
  
  lines(l_inf1_values,
        aucPlot[,7],
        lwd=4, col="orange", lty="longdash")
  
  lines(l_inf1_values,
        aucPlot[,8],
        lwd=4, col="orange", lty="twodash")
  
  
  lines(l_inf1_values,
        aucPlot[,9],
        lwd=4, col="brown", lty="dotted")
  
  lines(l_inf1_values,
        aucPlot[,10],
        lwd=4, col="black", lty="twodash")
  
  axis(1, at = l_inf1_values, cex.axis = 1.8, las = 1)
  axis(2, at = seq(0,1,0.2),
       labels = seq(0,1,0.2),
       cex.axis = 1.8, las = 1)
  
  box()
  
  dev.off()
  
  
  ######Plots pour l > 1 qui varie
  
  r = 30
  
  l_up1_values = seq(1,1e2,by = 1)
  
  erreursSigmaPlot = array(0,dim = c(length(l_up1_values),11))
  faux_positifsPlot= array(0,dim = c(length(l_up1_values),11))
  faux_negatifsPlot= array(0,dim = c(length(l_up1_values),11))
  ariPlot = array(0,dim = c(length(l_up1_values),11))
  aucPlot = array(0,dim = c(length(l_up1_values),11))
  
  
  for (s in seq_along(l_up1_values)){
    
    k = 0 
    
    l = l_up1_values[s]
    rho1 = 0.3
    
    
    
    for(j in seq_along(methodes)){
      
      
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
      
      
      erreursSigmaPlot[s,j] = crit$erreurFrob
      faux_negatifsPlot[s,j] = crit$FN
      faux_positifsPlot[s,j] = crit$FP
      ariPlot[s,j] = crit$ARI
      aucPlot[s,j] = crit$AUC
      
    }}  
  setwd("~/figures")
  # 
  # erreursSigmaPlot = erreursSigmaPlot[1:1e2,]
  # faux_negatifsPlot = faux_negatifsPlot[1:1e2,]
  # faux_positifsPlot = faux_positifsPlot[1:1e2,]
  # ariPlot = ariPlot[1:1e2,]
  # aucPlot = aucPlot[1:1e2,]
  
  file <- paste0("scen-lup1",".pdf")
  
  pdf(file, width = 25, height = 4)
  par(mfrow = c(1,5), mar = c(4,4,2,1))
  ############################################################
  ################ ERREUR DE FROBENIUS #######################
  ############################################################
  
  plot(l_up1_values, erreursSigmaPlot[,5],
       type="l", lwd=4, col="red",
       log="y",
       ylim=c(1e-1,1e2),
       xlab="", ylab="",
       xaxt="n", yaxt="n")
  
  lines(l_up1_values, erreursSigmaPlot[,3],
        lwd=4, col="blue", lty="dashed")
  
  lines(l_up1_values, erreursSigmaPlot[,1],
        lwd=4, col="darkgreen", lty="dotted")
  
  lines(l_up1_values, erreursSigmaPlot[,7],
        lwd=4, col="orange", lty="dotted")
  
  lines(l_up1_values, erreursSigmaPlot[,10],
        lwd=4, col="black", lty="dotted")
  
  lines(l_up1_values, erreursSigmaPlot[,9],
        lwd=4, col="brown", lty="dotted")
  
  axis(1, at= l_up1_values, cex.axis = 1.8,las=1)
  axis(2,cex.axis = 1.8,
       at=10^seq(-1,2),
       labels=parse(text=paste0("10^", -1:2)),
       las=1)
  box()
  
  ############################################################
  ################### FAUX NEGATIFS ##########################
  ############################################################
  pseudo_log <- function(x) log10(1 + x)
  
  plot(l_up1_values,
       pseudo_log(faux_negatifsPlot[,5] / ((30/100)*n) * 100),
       type = "l", lwd = 4, col = "red",
       ylim = range(pseudo_log(c(0, 100))),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  
  lines(l_up1_values,
        pseudo_log(faux_negatifsPlot[,6] / ((30/100)*n) * 100),
        lwd = 4, col = "red", lty = "longdash")
  
  
  
  
  
  lines(l_up1_values,
        pseudo_log(faux_negatifsPlot[,1] / ((30/100)*n) * 100),
        lwd = 4, col = "darkgreen", lty = "dotted")
  
  lines(l_up1_values,
        pseudo_log(faux_negatifsPlot[,2] / ((30/100)*n) * 100),
        lwd = 4, col = "darkgreen", lty = "dotdash")
  
  lines(l_up1_values,
        pseudo_log(faux_negatifsPlot[,3] / ((30/100)*n) * 100),
        lwd = 4, col = "blue", lty = "longdash")
  
  lines(l_up1_values,
        pseudo_log(faux_negatifsPlot[,4] / ((30/100)*n) * 100),
        lwd = 4, col = "blue", lty = "twodash")
  
  lines(l_up1_values,
        pseudo_log(faux_negatifsPlot[,7] / ((30/100)*n) * 100),
        lwd = 4, col = "orange", lty = "longdash")
  
  lines(l_up1_values,
        pseudo_log(faux_negatifsPlot[,8] / ((30/100)*n) * 100),
        lwd = 4, col = "orange", lty = "twodash")
  
  
  axis(1, at = l_up1_values, las = 1, cex.axis = 1.8)
  
  axis(2,
       at = log10(1 + c(0, 1, 10, 100)),
       labels = c("0", "1", "10", "100"),
       las = 1,
       cex.axis = 1.8)
  
  box()
  ############################################################
  ################### FAUX POSITIFS ##########################
  ############################################################
  
  plot(l_up1_values,
       faux_positifsPlot[,5] /
         ((1 - 30/100)*n)*100,
       type="l", lwd=4, col="red",
       ylim=c(0,20),
       xlab="", ylab="",
       xaxt="n", yaxt="n")
  
  lines(l_up1_values,
        faux_positifsPlot[,6] /
          ((1 - 30/100)*n)*100,
        lwd=4, col="red", lty="longdash")
  
  
  lines(l_up1_values,
        faux_positifsPlot[,1] /
          ((1 - 30/100)*n)*100,
        lwd=4, col="darkgreen", lty="dotted")
  
  
  
  lines(l_up1_values,
        faux_positifsPlot[,2] /
          ((1 - 30/100)*n)*100,
        lwd=4, col="darkgreen", lty="dotdash")
  
  
  lines(l_up1_values,
        faux_positifsPlot[,3] /
          ((1 - 30/100)*n)*100,
        lwd=4, col="blue", lty="longdash")
  
  lines(l_up1_values,
        faux_positifsPlot[,4] /
          ((1 - 30/100)*n)*100,
        lwd=4, col="blue", lty="twodash")
  
  lines(l_up1_values,
        faux_positifsPlot[,11] /
          ((1 - 30/100)*n)*100,
        lwd=4, col="purple4", lty="longdash")
  
  lines(l_up1_values,
        faux_positifsPlot[,7] /
          ((1 - 30/100)*n)*100,
        lwd=4, col="orange", lty="longdash")
  
  lines(l_up1_values,
        pseudo_log(faux_positifsPlot[,8] /
                     (1 - 30/100)*n)*100,
        lwd=4, col="orange", lty="twodash")
  
  
  lines(l_up1_values,
        faux_positifsPlot[,9] /
          ((1 - 30/100)*n)*100,
        lwd=4, col="brown", lty="dotted")
  
  lines(l_up1_values,
        faux_positifsPlot[,10] /
          ((1 - 30/100)*n)*100,
        lwd=4, col="black", lty="twodash")
  
  
  
  axis(1, at = l_up1_values, las = 1, cex.axis = 1.8)
  axis(2,
       at = c(0, 5, 10, 15, 20),
       labels = c("0", "5", "10", "15", "20"),
       las = 1,
       cex.axis = 1.8)
  box()
  
  ############################################################
  ######################## ARI ###############################
  ############################################################
  
  plot(l_up1_values,
       ariPlot[,5],
       type="l", lwd=3, col="red",
       ylim=c(-1,1),
       xlab="", ylab="",
       xaxt="n", yaxt="n")
  lines(l_up1_values,
        ariPlot[,6],
        lwd=4, col="red", lty="longdash")
  
  
  lines(l_up1_values,
        ariPlot[,1],
        lwd=4, col="darkgreen", lty="dotted")
  
  
  
  lines(l_up1_values,
        ariPlot[,2],
        lwd=4, col="darkgreen", lty="dotdash")
  
  
  lines(l_up1_values,
        ariPlot[,3],
        lwd=4, col="blue", lty="longdash")
  
  lines(l_up1_values,
        ariPlot[,4],
        lwd=4, col="blue", lty="twodash")
  
  lines(l_up1_values,
        ariPlot[,11],
        lwd=4, col="purple4", lty="longdash")
  
  lines(l_up1_values,
        ariPlot[,7],
        lwd=4, col="orange", lty="longdash")
  
  lines(l_up1_values,
        ariPlot[,8],
        lwd=4, col="orange", lty="twodash")
  
  
  lines(l_up1_values,
        ariPlot[,9],
        lwd=4, col="brown", lty="dotted")
  
  lines(l_up1_values,
        ariPlot[,10],
        lwd=4, col="black", lty="twodash")
  
  axis(1, at = l_up1_values, las = 1, cex.axis = 1.8)
  axis(2, at = seq(-1,1,0.5),
       labels = seq(-1,1,0.5),
       las = 1, cex.axis = 1.8)
  box()
  ############################################################
  ######################## AUC ###############################
  ############################################################
  
  
  
  plot(l_up1_values,
       aucPlot[,5],
       type="l", lwd=3, col="red",
       ylim=c(0,1),
       xlab="", ylab="",
       xaxt="n", yaxt="n")
  lines(l_up1_values,
        aucPlot[,6],
        lwd=4, col="red", lty="longdash")
  
  
  lines(l_up1_values,
        aucPlot[,1],
        lwd=4, col="darkgreen", lty="dotted")
  
  
  
  lines(l_up1_values,
        aucPlot[,2],
        lwd=4, col="darkgreen", lty="dotdash")
  
  
  lines(l_up1_values,
        aucPlot[,3],
        lwd=4, col="blue", lty="longdash")
  
  lines(l_up1_values,
        aucPlot[,4],
        lwd=4, col="blue", lty="twodash")
  
  lines(l_up1_values,
        aucPlot[,11],
        lwd=4, col="purple4", lty="longdash")
  
  lines(l_up1_values,
        aucPlot[,7],
        lwd=4, col="orange", lty="longdash")
  
  lines(l_up1_values,
        aucPlot[,8],
        lwd=4, col="orange", lty="twodash")
  
  
  lines(l_up1_values,
        aucPlot[,9],
        lwd=4, col="brown", lty="dotted")
  
  lines(l_up1_values,
        aucPlot[,10],
        lwd=4, col="black", lty="twodash")
  
  axis(1, at = l_up1_values, cex.axis = 1.8, las = 1)
  axis(2, at = seq(0,1,0.2),
       labels = seq(0,1,0.2),
       cex.axis = 1.8, las = 1)
  
  box()
  
  dev.off()

  
  
###Plots pour rho1 qui varie

  
  r = 30
  
  rho1_values = seq(-0.99,0.99,by = 0.01)
  
  erreursSigmaPlot = array(0,dim = c(length(rho1_values),11))
  faux_positifsPlot= array(0,dim = c(length(rho1_values),11))
  faux_negatifsPlot= array(0,dim = c(length(rho1_values),11))
  ariPlot = array(0,dim = c(length(rho1_values),11))
  aucPlot = array(0,dim = c(length(rho1_values),11))
  
  
  for (s in seq_along(rho1_values)){
    
    k = 0 
    
    l = 1
    rho1 = rho1_values[s]
    
    
    
    for(j in seq_along(methodes)){
      
      
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
      
      
      erreursSigmaPlot[s,j] = crit$erreurFrob
      faux_negatifsPlot[s,j] = crit$FN
      faux_positifsPlot[s,j] = crit$FP
      ariPlot[s,j] = crit$ARI
      aucPlot[s,j] = crit$AUC
      
    }}  
  setwd("~/figures")
  # 
  # erreursSigmaPlot = erreursSigmaPlot[1:1e2,]
  # faux_negatifsPlot = faux_negatifsPlot[1:1e2,]
  # faux_positifsPlot = faux_positifsPlot[1:1e2,]
  # ariPlot = ariPlot[1:1e2,]
  # aucPlot = aucPlot[1:1e2,]
  
  file <- paste0("scen-rho1",".pdf")
  
  pdf(file, width = 25, height = 4)
  par(mfrow = c(1,5), mar = c(4,4,2,1))
  ############################################################
  ################ ERREUR DE FROBENIUS #######################
  ############################################################
  
  plot(rho1_values, erreursSigmaPlot[,5],
       type="l", lwd=4, col="red",
       log="y",
       ylim=c(1e-1,1e2),
       xlab="", ylab="",
       xaxt="n", yaxt="n")
  
  lines(rho1_values, erreursSigmaPlot[,3],
        lwd=4, col="blue", lty="dashed")
  
  lines(rho1_values, erreursSigmaPlot[,1],
        lwd=4, col="darkgreen", lty="dotted")
  
  lines(rho1_values, erreursSigmaPlot[,7],
        lwd=4, col="orange", lty="dotted")
  
  lines(rho1_values, erreursSigmaPlot[,10],
        lwd=4, col="black", lty="dotted")
  
  lines(rho1_values, erreursSigmaPlot[,9],
        lwd=4, col="brown", lty="dotted")
  
  axis(1, at= rho1_values, cex.axis = 1.8,las=1)
  axis(2,cex.axis = 1.8,
       at=10^seq(-1,2),
       labels=parse(text=paste0("10^", -1:2)),
       las=1)
  box()
  
  ############################################################
  ################### FAUX NEGATIFS ##########################
  ############################################################
  pseudo_log <- function(x) log10(1 + x)
  
  plot(rho1_values,
       pseudo_log(faux_negatifsPlot[,5] / ((30/100)*n) * 100),
       type = "l", lwd = 4, col = "red",
       ylim = range(pseudo_log(c(0, 100))),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  
  lines(rho1_values,
        pseudo_log(faux_negatifsPlot[,6] / ((30/100)*n) * 100),
        lwd = 4, col = "red", lty = "longdash")
  
  
  
  
  
  lines(rho1_values,
        pseudo_log(faux_negatifsPlot[,1] / ((30/100)*n) * 100),
        lwd = 4, col = "darkgreen", lty = "dotted")
  
  lines(rho1_values,
        pseudo_log(faux_negatifsPlot[,2] / ((30/100)*n) * 100),
        lwd = 4, col = "darkgreen", lty = "dotdash")
  
  lines(rho1_values,
        pseudo_log(faux_negatifsPlot[,3] / ((30/100)*n) * 100),
        lwd = 4, col = "blue", lty = "longdash")
  
  lines(rho1_values,
        pseudo_log(faux_negatifsPlot[,4] / ((30/100)*n) * 100),
        lwd = 4, col = "blue", lty = "twodash")
  
  lines(rho1_values,
        faux_negatifsPlot[,7] /
          ((30/100)*n)*100,
        lwd=4, col="orange", lty="longdash")
  
  lines(rho1_values,
        pseudo_log(faux_negatifsPlot[,8] /
                     (30/100)*n)*100,
        lwd=4, col="orange", lty="twodash")
  
  
  axis(1, at = rho1_values, las = 1, cex.axis = 1.8)
  
  axis(2,
       at = log10(1 + c(0, 1, 10, 100)),
       labels = c("0", "1", "10", "100"),
       las = 1,
       cex.axis = 1.8)
  
  box()
  ############################################################
  ################### FAUX POSITIFS ##########################
  ############################################################
  
  plot(rho1_values,
       faux_positifsPlot[,5] /
         ((1 - 30/100)*n)*100,
       type="l", lwd=4, col="red",
       ylim=c(0,20),
       xlab="", ylab="",
       xaxt="n", yaxt="n")
  
  lines(rho1_values,
        faux_positifsPlot[,6] /
          ((1 - 30/100)*n)*100,
        lwd=4, col="red", lty="longdash")
  
  
  lines(rho1_values,
        faux_positifsPlot[,1] /
          ((1 - 30/100)*n)*100,
        lwd=4, col="darkgreen", lty="dotted")
  
  
  
  lines(rho1_values,
        faux_positifsPlot[,2] /
          ((1 - 30/100)*n)*100,
        lwd=4, col="darkgreen", lty="dotdash")
  
  
  lines(rho1_values,
        faux_positifsPlot[,3] /
          ((1 - 30/100)*n)*100,
        lwd=4, col="blue", lty="longdash")
  
  lines(rho1_values,
        faux_positifsPlot[,4] /
          ((1 - 30/100)*n)*100,
        lwd=4, col="blue", lty="twodash")
  
  lines(rho1_values,
        faux_positifsPlot[,11] /
          ((1 - 30/100)*n)*100,
        lwd=4, col="purple4", lty="longdash")
  
  lines(rho1_values,
        faux_positifsPlot[,7] /
          ((1 - 30/100)*n)*100,
        lwd=4, col="orange", lty="longdash")
  
  lines(rho1_values,
        pseudo_log(faux_positifsPlot[,8] /
                     (1 - 30/100)*n)*100,
        lwd=4, col="orange", lty="twodash")
  
  
  lines(rho1_values,
        faux_positifsPlot[,9] /
          ((1 - 30/100)*n)*100,
        lwd=4, col="brown", lty="dotted")
  
  lines(rho1_values,
        faux_positifsPlot[,10] /
          ((1 - 30/100)*n)*100,
        lwd=4, col="black", lty="twodash")
  
  
  
  axis(1, at = rho1_values, las = 1, cex.axis = 1.8)
  axis(2,
       at = c(0, 5, 10, 15, 20),
       labels = c("0", "5", "10", "15", "20"),
       las = 1,
       cex.axis = 1.8)
  box()
  
  ############################################################
  ######################## ARI ###############################
  ############################################################
  
  plot(rho1_values,
       ariPlot[,5],
       type="l", lwd=3, col="red",
       ylim=c(-1,1),
       xlab="", ylab="",
       xaxt="n", yaxt="n")
  lines(rho1_values,
        ariPlot[,6],
        lwd=4, col="red", lty="longdash")
  
  
  lines(rho1_values,
        ariPlot[,1],
        lwd=4, col="darkgreen", lty="dotted")
  
  
  
  lines(rho1_values,
        ariPlot[,2],
        lwd=4, col="darkgreen", lty="dotdash")
  
  
  lines(rho1_values,
        ariPlot[,3],
        lwd=4, col="blue", lty="longdash")
  
  lines(rho1_values,
        ariPlot[,4],
        lwd=4, col="blue", lty="twodash")
  
  lines(rho1_values,
        ariPlot[,11],
        lwd=4, col="purple4", lty="longdash")
  
  lines(rho1_values,
        ariPlot[,7],
        lwd=4, col="orange", lty="longdash")
  
  lines(rho1_values,
        ariPlot[,8],
        lwd=4, col="orange", lty="twodash")
  
  
  lines(rho1_values,
        ariPlot[,9],
        lwd=4, col="brown", lty="dotted")
  
  lines(rho1_values,
        ariPlot[,10],
        lwd=4, col="black", lty="twodash")
  
  axis(1, at = l_up1_values, las = 1, cex.axis = 1.8)
  axis(2, at = seq(-1,1,0.5),
       labels = seq(-1,1,0.5),
       las = 1, cex.axis = 1.8)
  box()
  ############################################################
  ######################## AUC ###############################
  ############################################################
  
  
  
  plot(rho1_values,
       aucPlot[,5],
       type="l", lwd=3, col="red",
       ylim=c(0,1),
       xlab="", ylab="",
       xaxt="n", yaxt="n")
  lines(rho1_values,
        aucPlot[,6],
        lwd=4, col="red", lty="longdash")
  
  
  lines(rho1_values,
        aucPlot[,1],
        lwd=4, col="darkgreen", lty="dotted")
  
  
  
  lines(rho1_values,
        aucPlot[,2],
        lwd=4, col="darkgreen", lty="dotdash")
  
  
  lines(rho1_values,
        aucPlot[,3],
        lwd=4, col="blue", lty="longdash")
  
  lines(rho1_values,
        aucPlot[,4],
        lwd=4, col="blue", lty="twodash")
  
  lines(rho1_values,
        aucPlot[,11],
        lwd=4, col="purple4", lty="longdash")
  
  lines(rho1_values,
        aucPlot[,7],
        lwd=4, col="orange", lty="longdash")
  
  lines(rho1_values,
        aucPlot[,8],
        lwd=4, col="orange", lty="twodash")
  
  
  lines(rho1_values,
        aucPlot[,9],
        lwd=4, col="brown", lty="dotted")
  
  lines(rho1_values,
        aucPlot[,10],
        lwd=4, col="black", lty="twodash")
  
  axis(1, at = rho1_values, cex.axis = 1.8, las = 1)
  axis(2, at = seq(0,1,0.2),
       labels = seq(0,1,0.2),
       cex.axis = 1.8, las = 1)
  
  box()
  
  dev.off()
  
  

#######################Trajectoires##########################

methodes_online = c("SampleNaiveQuantonlinecorr","SampleNaivewithoutonlinequantilecorr","OnlineUsQuantonlinecorr","OnlineUswithoutQuantonlinecorr","StreamingUsonlineQuantcorr","StreamingUswithoutQuantonlinecorr","Oracle")


for(sc in scenarios){
  
  k = sc$k
  l = sc$l
  rho1 = sc$rho1
  

  
  for (j in seq_along(rList[1:9])){
  
    
  r = rList[j]
    
  outlabTraj= array(0, dim = c(n,7))
  
  
  setwd(SimDir)  
  
  dataFile <- paste0('SimData-d', d, '-n', n, '-k', k, '-l', l, '-rho', rho1,'-r',r ,"-sim",sim,".RData")
  
  print(dataFile)
  
  load(dataFile)
  labels = data$labelsVrais
  
  
  for(s in seq_along(methodes_online)){
  
  methode = methodes_online[s]

  
  setwd(resAlgo)
  
  fitFile <- paste0('Fit-',methode,'-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
  
  load(fitFile)
  
  outlabTraj[,s] = resultats$outliers_labels
  
  }
  
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
  



