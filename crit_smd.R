smd_data_dir = "~/smd_data_dir"
res_SMD = "~/res_SMD"
crit_SMD = "~/criteres_SMD"

methodes = c("Oracle","SampleNaiveQuantonlinecorr","SampleNaivewithoutQuantonlinecorr","MCD","Offline-without_onlinequantile","Offline-without_onlinequantile","OnlineUsQuantonlinecorr","OnlineUsWithoutQuantonlinecorr","StreamingUsonlineQuantcorr","StreamingUs_without_onlineQuantcorr","OGK")

for(j in 1:28){
  
  
    
  setwd(smd_data_dir)
  
  
  data_smd_mach = paste0("data_machine-",j,".RData")
  
  
  load(data_smd_mach)  
  
  setwd(res_SMD)  

  
  fitFile <- paste0('Fit-Oracle-',"machine-",j,".RData")
  
  load(fitFile)
  
  Sigma_trueCov = resultats$variance
    
  
    
      for(methode in methodes){
        
        
        
        for(methode in methodes){
        
        setwd(res_SMD)
        
        fitFile = paste0("Fit-",methode,"-machine-",j,".RData")
        
        crit = compute_criteres(variance = resultats$variance, outlab = resultats$outliers_labels, distances = resultats$distances, labels = as.numeric(labels), SigmaTrue = Sigma0, r = r,Sigma_ref = Sigma_trueCov,smd = TRUE)
      
        
        }
        setwd(crit_SMD)
        
        critFile <- paste0("Crit-",methode,"-machine-",j,".RData")
        save(  crit,file = critFile)
      }
     
  
  
}

