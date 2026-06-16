smd_data_dir = "~/smd_data_dir"
res_SMD = "~/res_SMD"

for(j in 1:28){
  setwd(smd_data_dir)
  
  
  data_smd_mach = paste0("data_machine-",j,".RData")
  

  load(data_smd_mach)  
  
  print(data_smd_mach)
  
  setwd(res_SMD)
  
  
  #Calcul taille de batch
  
  possible_div <- which(nrow(Z) %% 1:nrow(Z) == 0)
  
  # ne garder que ceux <= ncol(Z)
  possible_div <- possible_div[possible_div <= ncol(Z)]
  
  # prendre le plus grand
  batchStrm <- max(possible_div)
  
  ###########Extraction moyenne et covariance sans outliers
  
  meanTrueCov = colMeans(Z_clean)
  
  Sigma_trueCov = cov(Z_clean)
  
  invSigmaTrueCov = inv_safe(Sigma_trueCov)
  
  Sigma_trueMCD = covMcd(Z_clean)$cov
  
  Sigma_trueOGK = covOGK(Z_clean,sigmamu = scaleTau2)$cov
  
  Sigma_trueOffl = offlineRobustVariance(Z_clean)$variance
  
  Sigma_trueOnl = onlineRobustVariance(Z_clean,computeOutliers = TRUE,cutoff=0.95,cutinit=0.6,nDataInit = 1e3,c_m=1,batch = 1)$variance
  
  Sigma_trueStrm = onlineRobustVariance(Z_clean,computeOutliers = TRUE,cutoff=0.95,cutinit=0.6,nDataInit = 1e3,c_m=1,batch = batchStrm)$variance

  ################Distances without outliers and with sample mean and covariance without ourliers####################
  
  distrue = rep(0,nrow(Z))
  
  for(i in 1:nrow(Z))
  {
    distrue[i] = t(Z[i,] - meanTrueCov)%*%invSigmaTrueCov%*%(Z[i,] - meanTrueCov)
    
  }
  distrueClean = rep(0,nrow(Z_clean))
  for (s in 1:nrow(Z_clean))
  {
    distrueClean[s] = t(Z_clean[s,] - meanTrueCov)%*%invSigmaTrueCov%*%(Z_clean[s,] - meanTrueCov)
  }
  
#####################Cutoffs##################################
  
  quantoracle = quantile(distrueClean,.95)
  
  cut_off_ogk = qchisq(.95,df = ncol(Z))
  cut_off_mcd = qchisq(.95,df = ncol(Z))
  
  
######################################################algorithms#######################  
  
  
  fitFile <- paste0('FitOracle-',"-machine-",j,".RData")
  
  
  oultoracle = rep(0,nrow(Z))
  
  for(s in 1:nrow(Z)){
    if(distrue[s] > quantoracle){
      outloracle[s] = 1
    }
  }
  
save(distrue,distrueClean,outloracle,file = fitFile)    
  temps_covonline = system.time(
    {
      resSamplecov= SampleCovOnline(Z,quantcutoff = TRUE,nDataInit = 1e3,cutoffquant =  .95,c_m = 1)
    })

  fitFile <- paste0('FitSampleNaiveQuantonlinecorr-',"-machine-",j,"-cutoff",.95,"Ninit-",1e3,"-cm",1,".RData")
  
  save(temps_covonline,resSamplecov,file = fitFile)
  
  
  temps_covonline = system.time(
    {
      resSamplecov= SampleCovOnline(Z,quantcutoff = FALSE,nDataInit = 1e3,cutoffquant =  .95,c_m = 1)
    })

  fitFile <- paste0('FitSampleNaivewithoutQuantonlinecorr-',"-machine-",j,"-cutoff",.95,"Ninit-",1e3,"-cm",1,".RData")
  
  save(temps_covonline,resSamplecov,file = fitFile)
  
  
  fitFile <- paste0('FitMCD-machine-',j,".RData")
  
  distmcd = rep(0,nrow(Z))
  
  outlmcd = rep(0,nrow(Z))
  
  
  temps_mcd = system.time(
    
    {
      
      resmcd = covMcd(Z)
      invSigmaMCD = inv_safe(resmcd$cov)
      
      for(i in 1:nrow(Z))
      {
        distmcd[i] = t(Z[i,] - (resmcd$center))%*%invSigmaMCD%*%(Z[i,] - (resmcd$center))
        
        if (distmcd[s] > cut_off_mcd){
          outlmcd[i] = 1
          }
      }
    }
  )

  

  save(temps_mcd,resmcd,distmcd,outlmcd,file = fitFile)
  
  temps_offline = system.time(
    {
      res0 <- offlineRobustVariance(Z,computeOutliers = TRUE,cutinit=0.6,c_m=1)
    }
  )
  
   fitFile <- paste0('FitOffline-withonlinequantile-machine-',j,".RData")
   
   save(temps_offline,res0,file = fitFile)
  # 
   temps_offline = system.time(
     {
       res0 <- offlineRobustVariance_old(Z)
   
     }
   )
   
   fitFile <- paste0('FitOffline-without_onlinequantile-machine-',j,".RData")
   
   save(temps_offline,res0,file = fitFile)
  
  
  ###############################################Online us#########################################
  
  fitFile <- paste0('FitOnlineUsQuantonlinecorr-machine-',j,"-cutoff",.95,"Ninit-",1e3,"-batch",1,"-cm",1,".RData")
  
  
  temps_online = system.time(
    {
      
      resUsOnline= onlineRobustVariance(Z,computeOutliers = TRUE,batch = 1,cutoff=cut,cutinit=0.6,nDataInit = 1e3,c_m= 1)

      
    })
  

  save(resUsOnline,temps_online,file = fitFile)
  
  fitFile <- paste0('FitOnlineUsWithoutQuantonlinecorr-machine-',j,"-cutoff",.95,"Ninit-",1e3,"-batch",1,"-cm",1,".RData")
  
   temps_online = system.time(
     {
       
       resUsOnline= onlineRobustVariance_old(Z,batch = 1,nDataInit = 1e3,computeOutliers = TRUE)
       
     })
  # 
  # 
  save(resUsOnline,temps_online,file = fitFile)
  
  
  
  ###############################################Streaming us#########################################
  
  fitFile <- paste0('FitStreamingUsonlineQuantcorr-machine-', j,"-cutoff",.95,"Ninit-",1e3,"-cm",1,".RData")
  
  
  temps_streaming = system.time(
    {
      
      resStrm <- onlineRobustVariance(Z,computeOutliers = TRUE,cutoff=cut,cutinit=0.6,nDataInit = 1e3,c_m= 1,batch = batchStrm)
      
    })
  
  
  save(resStrm,temps_streaming,file = fitFile)
  
  
  
  # 
  fitFile <- paste0('FitStreamingUs_without_onlineQuantcorr-machine-', j,"-cutoff",.95,"Ninit-",1e3,"-cm",1,".RData")
  # 
   temps_streaming = system.time(
     {
       
       resStrm <- onlineRobustVariance_old(Z,computeOutliers = TRUE,nDataInit = 1e3,batch = batchStrm)
       
       
     })
   
   
   save(resStrm,temps_streaming,file = fitFile)
  
  outlogk = rep(0,nrow(Z))
  
  temps_ogk = system.time(
    {
      resogk = covOGK(Z, sigmamu = scaleTau2)
      for (s in 1:nrow(Z))
      if(resogk$distances[s] > cut_off_ogk){outlogk[s] = 1}    
    }
  )

  fitFile <- paste0('FitOGK-machine-',j,".RData")
  
  
  save(resogk,outlogk,temps_ogk,file = fitFile)
  
  

  
  
}