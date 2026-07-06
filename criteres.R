################Calcul criteres###################
methodes = c("SampleNaiveQuantonlinecorr","SampleNaivewithoutonlinequantilecorr","OnlineUsQuantonlinecorr","OnlineUswithoutQuantonlinecorr","StreamingUsonlineQuantcorr","StreamingUswithoutQuantonlinecorr","OfflinewithQuantcorr","OfflineUswithoutQuantcorr","OGK","MCD","Oracle")
methodes_add  = c("oracleRD","OracleQC","SampleRaw","OnlRaw","StrmRaw","OfflRaw","OGKRD","OGKQC","MCDRD","MCDQC")
for(sim in 1:1){
for (sc in scen_conc){
  k = sc$k
  l = sc$l
  rho1 = sc$rho1
  
  
  
  for (m in seq_along(rList[1:9])){
    
    r = rList[m]
     
    for(methode in methodes){
  
    setwd(SimDir)
    
    dataFile <- paste0('SimData-d', d, '-n', n, '-k', k, '-l', l, '-rho', rho1,'-r',r ,"-sim",sim,".RData")
    
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

}

#######################Calcul moyenne erreur norme de Frobenius faux positifs faux négatifs, ari, auc
#################### Moyenne des critères sur les simulations ####################

setwd(criteres)

for (sc in scenarios){
  
  k <- sc$k
  l <- sc$l
  rho1 <- sc$rho1
  
  for (r in rList[1:9]){
    
    for (methode in methodes){
      
      erreurFrob <- 0
      FP <- 0
      FN <- 0
      ARI <- 0
      AUC <- 0
      
      for (sim in 1:simNb){
        
        critFile <- paste0(
          "Crit-", methode,
          "-d", d,
          "-n", n,
          "-k", k,
          "-l", l,
          "-rho", rho1,
          "-r", r,
          "-sim", sim,
          ".RData"
        )
        
        
        load(critFile)
        
        erreurFrob <- erreurFrob + crit$erreurFrob
        FP <- FP + crit$FP
        FN <- FN + crit$FN
        ARI <- ARI + crit$ARI
        AUC <- AUC + crit$AUC
      }
      
      crit_mean <- list(
        erreurFrob = erreurFrob / simNb,
        FP = FP / simNb,
        FN = FN / simNb,
        ARI = ARI / simNb,
        AUC = AUC / simNb
      )
      
      save(
        crit_mean,
        file = paste0(
          "Crit-", methode,
          "-d", d,
          "-n", n,
          "-k", k,
          "-l", l,
          "-rho", rho1,
          "-r", r,
          "-mean.RData"
        )
      )
    }
  }
}