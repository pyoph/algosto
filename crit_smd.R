smd_data_dir = "~/smd_data_dir"
res_SMD = "~/res_SMD"
crit_SMD = "~/criteres_SMD"

for (j in 1:28){
  
  setwd(smd_data_dir)
  
  
  data_smd_mach = paste0("data_machine-",j,".RData")
  
  
  load(data_smd_mach)  

  Sigma_trueCov = cov(Z_clean)
  
  
  setwd(res_SMD)

  
  
      
  fitFile <- paste0('FitSampleNaiveQuantonlinecorr-',"-machine-",j,"-cutoff",.95,"Ninit-",1e3,"-cm",1,".RData")
  
  load(fitFile)
  
  setwd(crit_SMD)
  
  critFile <- paste0('CritSampleNaiveQuantonlinecorr-',"-machine-",j,"-cutoff",.95,"Ninit-",1e3,"-cm",1,".RData")
  
  
  FrobeniusNormError = norm(resSamplecov$Sigma - Sigma_trueCov,"F")
  FP <- sum(resSamplecov$outliers_labels == 1 & labels == 0)
  FN <- sum(resSamplecov$outliers_labels == 0 & labels == 1)
  
  save(FrobeniusNormError,FP,FN,file = critFile)
  
  setwd(res_SMD)
  

  fitFile <- paste0('FitSampleNaiveQuantonlinecorr-',"-machine-",j,"-cutoff",.95,"Ninit-",1e3,"-cm",1,".RData")
  
  load(fitFile)
  
  ##############Sample cov online################
  
  setwd(crit_SMD)
  
  critFile <- paste0('CritSampleNaiveQuantonlinecorr-',"-machine-",j,"-cutoff",.95,"Ninit-",1e3,"-cm",1,".RData")
  
  
  FrobeniusNormError = norm(resSamplecov$Sigma - Sigma_trueCov,"F")
  FP <- sum(resSamplecov$outliers_labels == 1 & labels == 0)
  FN <- sum(resSamplecov$outliers_labels == 0 & labels == 1)
  
  save(FrobeniusNormError,FP,FN,file = critFile)
  
  setwd(res_SMD)
  
  fitFile <- paste0('FitSampleNaivewithoutQuantonlinecorr-',"-machine-",j,"-cutoff",.95,"Ninit-",1e3,"-cm",1,".RData")
    
  load(fitFile)
  
  CritFile <- paste0('CritSampleNaivewithoutQuantonlinecorr-',"machine-",j,"-cutoff",.95,"Ninit-",1e3,"-cm",1,".RData")
  
  
  FrobeniusNormError = norm(resSamplecov$Sigma - Sigma_trueCov,"F")
  FP <- sum(resSamplecov$outliers_labels == 1 & labels == 0)
  FN <- sum(resSamplecov$outliers_labels == 0 & labels == 1)
  
  critFile = paste0('CritSampleNaiveQuantonlinecorr-',"machine-",j,"-cutoff",.95,"Ninit-",1e3,"-cm",1,".RData")
  
  setwd(crit_SMD)
  
  save(FrobeniusNormError,FP,FN,file = critFile)
  
  ######################MCD#######################
  
  setwd(res_SMD)
  
  fitFile <- paste0('FitMCD-machine-',j,".RData")
  
  load(fitFile)
  
  FrobeniusNormError = norm(resmcd$cov - Sigma_trueCov,"F")
  FP <- sum(outlmcd == 1 & labels == 0)
  FN <- sum(outlmcd == 0 & labels == 1)
  
  critFile = paste0('CritMCD-machine-',j,".RData")
  
  setwd(crit_SMD)
  
  save(FrobeniusNormError,FP,FN,file = critFile)
  
  
  ######################OGK#######################
  
  setwd(res_SMD)
  
  fitFile <- paste0('FitOGK-machine-',j,".RData")
  
  load(fitFile)
  
  critFile = paste0('CritOGK-machine-',j,".RData")
  
  setwd(crit_SMD)
  
  
  FrobeniusNormError = norm(resogk$cov - Sigma_trueCov,"F")
  FP <- sum(outlogk == 1 & labels == 0)
  FN <- sum(outlogk == 0 & labels == 1)
  
  
  save(FrobeniusNormError,FP,FN,file = critFile)
  
  ##########################Online###############################
  
  setwd(res_SMD)
  
  fitFile <- paste0('FitOnlineUsQuantonlinecorr-machine-',j,"-cutoff",.95,"Ninit-",1e3,"-batch",1,"-cm",1,".RData")
  
  load(fitFile)
  
  setwd(crit_SMD)
  
  critFile = paste0('CritOnlineUsQuantonlinecorr-machine-',j,"-cutoff",.95,"Ninit-",1e3,"-batch",1,"-cm",1,".RData")
  
  FrobeniusNormError = norm(resUsOnline$variance - Sigma_trueCov,"F")
  FP <- sum(resUsOnline$outlier_labels == 1 & labels == 0)
  FN <- sum(resUsOnline$outlier_labels == 0 & labels == 1)
  
  save(FrobeniusNormError,FP,FN,file = critFile)

  setwd(res_SMD)
  
  fitFile <- paste0('FitOnlineUsWithoutQuantonlinecorr-machine-',j,"-cutoff",.95,"Ninit-",1e3,"-batch",1,"-cm",1,".RData")
  
  load(fitFile)
  
  setwd(crit_SMD)
  
  critFile = paste0('CritOnlineUswithoutQuantonlinecorr-machine-',j,"-cutoff",.95,"Ninit-",1e3,"-batch",1,"-cm",1,".RData")
  
  FrobeniusNormError = norm(resUsOnline$variance - Sigma_trueCov,"F")
  FP <- sum(resUsOnline$outlier_labels == 1 & labels == 0)
  FN <- sum(resUsOnline$outlier_labels == 0 & labels == 1)
  
  save(FrobeniusNormError,FP,FN,file = critFile)
  
  
  ##########################Streaming###############################
  
  setwd(res_SMD)
  fitFile <- paste0('FitStreamingUsonlineQuantcorr-machine-', j,"-cutoff",.95,"Ninit-",1e3,"-batch",batchStrm,"-cm",1,".RData")
  
  load(fitFile)
  
  setwd(crit_SMD)
  
  critFile = paste0('CritStreamingUsonlineQuantcorr-machine-', j,"-cutoff",.95,"Ninit-",1e3,"-batch",batchStrm,"-cm",1,".RData")
  
  FrobeniusNormError = norm(resUsStreaming$variance - Sigma_trueCov,"F")
  FP <- sum(resStrm$outlier_labels == 1 & labels == 0)
  FN <- sum(resStrm$outlier_labels == 0 & labels == 1)
  
  save(FrobeniusNormError,FP,FN,file = critFile)
  
  setwd(res_SMD)
  
  fitFile <- paste0('FitStreamingUs_without_onlineQuantcorr-machine-', j,"-cutoff",.95,"Ninit-",1e3,"-cm",1,".RData")
  
  load(fitFile)
  
  setwd(crit_SMD)
  
  critFile = paste0('CritStreamingUs_without_onlineQuantcorr-machine-', j,"-cutoff",.95,"Ninit-",1e3,"-cm",1,".RData")
  
  FrobeniusNormError = norm(resStrm$variance - Sigma_trueCov,"F")
  FP <- sum(resStrm$outlier_labels == 1 & labels == 0)
  FN <- sum(resStrm$outlier_labels == 0 & labels == 1)
  save(FrobeniusNormError,FP,FN,file = critFile)
  
  
  ##########################Offline################################
  
  
  setwd(res_SMD)
  
  fitFile <- paste0('FitOffline-withonlinequantile-machine-',j,".RData")
  
  load(fitFile)
  
  setwd(crit_SMD)
  
  critFile = paste0('CritOffline-withonlinequantile-machine-', j,".RData")
  
  FrobeniusNormError = norm(res0$variance - Sigma_trueCov,"F")
  FP <- sum(res0$outlier_labels == 1 & labels == 0)
  FN <- sum(resOffline$outlier_labels == 0 & labels == 1)
  
  save(FrobeniusNormError,FP,FN,file = critFile)
  
  setwd(res_SMD)
  
  fitFile <- paste0('FitOffline-without_onlinequantile-machine-',j,".RData")
  
  load(fitFile)
  
  setwd(crit_SMD)
  
  critFile = paste0('CritOffline-without_onlinequantile-machine-', j,".RData")
  
  FrobeniusNormError = norm(resOffline$variance - Sigma_trueCov,"F")
  FP <- sum(res0$outlier_labels == 1 & labels == 0)
  FN <- sum(res0$outlier_labels == 0 & labels == 1)
  
  save(FrobeniusNormError,FP,FN,file = critFile)
  
  
  
  
}