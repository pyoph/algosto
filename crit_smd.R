
methodes = c("SampleNaiveQuantonlinecorr","OnlineUsQuantonlinecorr","StreamingUsonlineQuantcorr","Offline-withonlinequantile","OGK","MCD")
methodes_add  = c("SampleNaivewithoutonlinequantilecorr","OnlineUswithoutQuantonlinecorr","StreamingUswithoutQuantonlinecorr","OfflineUswithoutQuantcorr","OracleRD","OracleQC","SampleRaw","OnlRaw","StrmRaw","OfflRaw","OGKRD","OGKQC","MCDRD","MCDQC")
methode_oracle = c("Oracle")

for(j in 1:56){
  
  
    
  setwd(smd_data_dir)
  
  
  data_smd_mach = paste0("data_machine-",j,".RData")
  
  
  load(data_smd_mach)  
  
  setwd(res_SMD)  

  
    
      for(methode in c(methodes,methodes_add,methode_oracle)){
        
        setwd(res_SMD)
        
        
        fitFile = paste0("Fit-",methode,"-machine-",j,".RData")
        
        load(fitFile)
        
        critFile <- paste0("Crit-",methode,"-machine-",j,".RData")
        
        print(critFile)
        
        if(methode %in% methodes){
        crit = compute_criteres(variance = resultats$variance, outlab = resultats$outliers_labels, distances = resultats$distances, labels = as.numeric(labels), SigmaTrue = resultats$variance_ref, r = 5)
      
        setwd(crit_SMD)
        
        save(  crit,file = critFile)
        }
        
        if(methode %in% c("Oracle"))
        {
          
          crit <- compute_criteres(
            variance  = Sigma0,
            outlab    = resultats$outliers_labels,
            distances = resultats$distances,
            labels    = as.numeric(labels),
            SigmaTrue = Sigma0,
            r         = 5
          )
          
          
          setwd(crit_SMD)
          
          
          
          
          save(  crit,    file = critFile)
          
        }
        if(methode %in% methodes_add)
        {
          
          crit <- compute_criteres(
            variance  = Sigma0,
            outlab    = resultats$outliers_labels,
            distances = rep(0,length(labels)),
            labels    = as.numeric(labels),
            SigmaTrue = Sigma0,
            r         = 5
          )
          
          setwd(crit_SMD)
          
          
          
          save(  crit,    file = critFile)
          
          
        }
        
      }
     
  
  
}

