smd_data_dir = "~/smd_data_dir"
res_SMD = "~/res_SMD"
crit_SMD = "~/criteres_SMD"

methodes = c("Oracle","SampleNaiveQuantonlinecorr","SampleNaivewithoutQuantonlinecorr","MCD","Offline-without_onlinequantile","Offline-without_onlinequantile","OnlineUsQuantonlinecorr","OnlineUsWithoutQuantonlinecorr","StreamingUsonlineQuantcorr","StreamingUs_without_onlineQuantcorr","OGK")

for(j in 11:11){
  
  
    
  setwd(smd_data_dir)
  
  
  data_smd_mach = paste0("data_machine-",j,".RData")
  
  
  load(data_smd_mach)  
  
  setwd(res_SMD)  

  
    
      for(methode in methodes){
        
        setwd(res_SMD)
        
        fitFile = paste0("Fit-",methode,"-machine-",j,".RData")
        
        load(fitFile)
        
        crit = compute_criteres(variance = resultats$variance, outlab = resultats$outliers_labels, distances = resultats$distances, labels = as.numeric(labels), SigmaTrue = resultats$variance_ref, r = 5)
      
        setwd(crit_SMD)
        
        critFile <- paste0("Crit-",methode,"-machine-",j,".RData")
        save(  crit,file = critFile)
      }
     
  
  
}

