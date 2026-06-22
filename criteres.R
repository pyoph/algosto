################Calcul criteres###################

methodes = c("SampleNaiveQuantonlinecorr","SampleNaivewithoutonlinequantilecorr","OnlineUsQuantonlinecorr","OnlineUswithoutQuantonlinecorr","StreamingUsonlineQuantcorr","StreamingUswithoutQuantonlinecorr","OfflinewithQuantcorr","OfflineUswithoutQuantcorr","OGK","MCD","Oracle")

for (sc in scenarios){
  k = sc$k
  l = sc$l
  rho1 = sc$rho1
  
  for (m in seq_along(rList[1:9])){
    

    for(methode in methodes){
    
    r = rList[m]
    
    setwd(SimDir)
    
    dataFile <- paste0('SimData-d', d, '-n', n, '-k', k, '-l', l, '-rho', rho1,'-r',r ,".RData")
    
    load(dataFile)
    
    setwd(resAlgo)
    
    
    fitFile <- paste0(
      "Fit-", methode,
      "-d", d,
      "-n", n,
      "-k", k,
      "-l", l,
      "-rho", rho1,
      "-r", r,
      "-sim", sim,
      ".RData"
    )
    print(fitFile)
    
    
    load(fitFile)
    
    crit = compute_criteres(variance = resultats$variance, outlab = resultats$outliers_labels, distances = resultats$distances, labels = as.numeric(data$labelsVrais), SigmaTrue = Sigma0, r = r)
    
    
    
    setwd(criteres)
    
    
    
    critFile <- paste0(
      'Crit-',methode,"-d" ,d,
      '-n', n,
      '-k', k,
      '-l', l,
      '-rho', rho1,
      '-r', r,
      '-sim', sim,
      ".RData"
    )
    save(  crit,    file = critFile)
    }
  }  
}



