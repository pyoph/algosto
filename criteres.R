################Calcul criteres###################
methodes = c("SampleNaiveQuantonlinecorr","OnlineUsQuantonlinecorr","StreamingUsonlineQuantcorr","OfflinewithQuantcorr","OGK","MCD")
methodes_add  = c("SampleNaivewithoutonlinequantilecorr","OnlineUswithoutQuantonlinecorr","StreamingUswithoutQuantonlinecorr","OfflineUswithoutQuantcorr","OracleRD","OracleQC","SampleRaw","OnlRaw","StrmRaw","OfflRaw","OGKRD","OGKQC","MCDRD","MCDQC")
for(sim in 1:10){
  for (sc in scen_strong_conc){
    k = sc$k
    l = sc$l
    rho1 = sc$rho1
    
    
    
    for (m in seq_along(rList[1:13])){
      
      r = rList[m]
      
      for(methode in c(methodes,methodes_add)){
        
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
        if(methode %in% methodes)
        {
          
          crit <- compute_criteres(
            variance  = resultats$variance,
            outlab    = resultats$outliers_labels,
            distances = resultats$distances,
            labels    = as.numeric(data$labelsVrais),
            SigmaTrue = Sigma0,
            r         = r
          )
          setwd(criteres)
          
          
          
          save(  crit,    file = critFile)
          
        }
        
        if(methode %in% c("Oracle"))
        {
          
          crit <- compute_criteres(
            variance  = Sigma0,
            outlab    = resultats$outliers_labels,
            distances = resultats$distances,
            labels    = as.numeric(data$labelsVrais),
            SigmaTrue = Sigma0,
            r         = r
          )
          
          
          setwd(criteres)
          
          
          
          save(  crit,    file = critFile)
          
        }
        if(methode %in% methodes_add)
        {
          
          crit <- compute_criteres(
            variance  = Sigma0,
            outlab    = resultats$outliers_labels,
            distances = rep(0,n),
            labels    = as.numeric(data$labelsVrais),
            SigmaTrue = Sigma0,
            r         = r
          )
          
          setwd(criteres)
          
          
          
          save(  crit,    file = critFile)
          
          
        }
        
      }
    }  
  }
  
}

#######################Calcul moyenne erreur norme de Frobenius faux positifs faux négatifs, ari, auc
#################### Moyenne des critères sur les simulations ####################

setwd(criteres)

for (sc in scen_strong_conc){
  
  k <- sc$k
  l <- sc$l
  rho1 <- sc$rho1
  
  for (r in rList[1:13]){
    
    for (methode in c(methodes,methodes_add)){
      
      erreurFrob <- 0
      FP <- 0
      FN <- 0
      ARI <- 0
      AUC <- 0
      prop_hors_diag = 0    
      for (sim in 1:10){
        
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
        if(methode %in% methodes)
        {   
          erreurFrob <- erreurFrob + crit$erreurFrob
          FP <- FP + crit$FP
          FN <- FN + crit$FN
          ARI <- ARI + crit$ARI
          AUC <- AUC + crit$AUC
          prop_hors_diag = prop_hors_diag + crit$prop_hors_diag
          
        }
        if(methode %in% c("Oracle"))
        {   
          #erreurFrob <- erreurFrob + crit$erreurFrob
          FP <- FP + crit$FP
          FN <- FN + crit$FN
          ARI <- ARI + crit$ARI
          AUC <- AUC + crit$AUC
          prop_hors_diag = prop_hors_diag + crit$prop_hors_diag
          
          
        }
        
        
        else if (methode %in% methodes_add){
          #erreurFrob <- erreurFrob + crit$erreurFrob
          FP <- FP + crit$FP
          FN <- FN + crit$FN
          ARI <- ARI + crit$ARI
          #AUC <- AUC + crit$AUC
          prop_hors_diag = prop_hors_diag + crit$prop_hors_diag
          
          
        }
        
      }
      
      
      
      if(methode %in% c("SampleNaiveQuantonlinecorr","OnlineUsQuantonlinecorr","StreamingUsonlineQuantcorr","OfflinewithQuantcorr","OGK","MCD"))
      {
        crit_mean <- list(
          erreurFrob = erreurFrob / simNb,
          FP = FP / simNb,
          FN = FN / simNb,
          ARI = ARI / simNb,
          AUC = AUC / simNb,
          prop_hors_diag = prop_hors_diag/simNb
          
          
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
      
      
      if(methode %in% c("Oracle"))
      {
        crit_mean <- list(
          #erreurFrob = erreurFrob / simNb,
          FP = FP / simNb,
          FN = FN / simNb,
          ARI = ARI / simNb,
          AUC = AUC / simNb,
          prop_hors_diag = prop_hors_diag/simNb
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
      
      
      if(!methode %in% c("SampleNaiveQuantonlinecorr","OnlineUsQuantonlinecorr","StreamingUsonlineQuantcorr","OfflinewithQuantcorr","OGK","MCD","Oracle"))
      {crit_mean <- list(
        # erreurFrob = erreurFrob / simNb,
        FP = FP / simNb,
        FN = FN / simNb,
        ARI = ARI / simNb,
        prop_hors_diag = prop_hors_diag/simNb
        #AUC = AUC / simNb
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
      )} 
    }
  }}