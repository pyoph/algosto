smd_data_dir = "~/smd_data_dir"
res_SMD = "~/res_SMD"

for(j in 1:56){
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
  
  
  fitFile <- paste0('Fit-Oracle-',"machine-",j,".RData")
  
  
  outloracle = rep(0,nrow(Z))
  
  for(s in 1:nrow(Z)){
    if(distrue[s] > quantoracle){
      outloracle[s] = 1
    }
  }

  resultats <- list(
    variance = Sigma_trueCov,
    variance_ref = Sigma_trueCov,
    outliers_labels = outloracle,
    distances = distrue,
    temps = 0
  )
  
  save(resultats,file = fitFile)
  
  check_fit(fitFile = fitFile,variance = Sigma_trueCov,outliers_labels = outloracle,distances = distrue)
  
  
    


  temps_covonline = system.time(
    {
      resSamplecov= SampleCovOnline(Z,quantcutoff = TRUE,nDataInit = 1e3,cutoffquant =  .95,c_m = 1)
    })

  fitFile <- paste0('Fit-SampleNaiveQuantonlinecorr-',"machine-",j,".RData")
  
  #########Calcul sans les outliers#####################

  
  
  resSamplecov_ref = SampleCovOnline(Z[labels == 0,],quantcutoff = TRUE,nDataInit = 1e3,cutoffquant =  .95,c_m = 1)

  
  
  fitFile <- paste0('Fit-SampleNaiveQuantonlinecorr-',"machine-",j,".RData")
  
  
  
  resultats <- list(
    variance = resSamplecov$Sigma,
    variance_ref =  resSamplecov_ref$Sigma,
    outliers_labels = resSamplecov$outliers_labels,
    distances = resSamplecov$distances,
    temps = temps_covonline
  )
  
  save(resultats,file = fitFile)
  
  check_fit(fitFile = fitFile,variance = resSamplecov$Sigma,outliers_labels = resSamplecov$outliers_labels,distances = resSamplecov$distances)
  
  
  
  temps_covonline = system.time(
    {
      resSamplecov= SampleCovOnline(Z,quantcutoff = FALSE,nDataInit = 1e3,cutoffquant =  .95,c_m = 1)
    })

  fitFile <- paste0('Fit-SampleNaivewithoutQuantonlinecorr',"-machine-",j,".RData")
  
  resSamplecov_ref = SampleCovOnline(Z[labels == 0,],quantcutoff = TRUE,nDataInit = 1e3,cutoffquant =  .95,c_m = 1)
  
  
  
  resultats <- list(
    variance = resSamplecov$Sigma,
    variance_ref = resSamplecov_ref$Sigma,
    resSamplecov_ref$Sigma,
    outliers_labels = resSamplecov$outliers_labels,
    distances = resSamplecov$distances,
    temps = temps_covonline
  )
  
  save(resultats,file = fitFile)
  
  check_fit(fitFile = fitFile,variance = resSamplecov$Sigma,outliers_labels = resSamplecov$outliers_labels,distances = resSamplecov$distances)
  
  fitFile <- paste0('Fit-MCD-machine-',j,".RData")
  
  distmcd = rep(0,nrow(Z))
  
  outlmcd = rep(0,nrow(Z))

  resmcd_ref = covMcd(Z[labels == 0,])     
  
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
  
  resultats <- list(
    variance = resmcd$cov,
    variance_ref = resmcd_ref$cov,
    outliers_labels = outlmcd,
    distances = distmcd,
    temps = temps_mcd
  )
  
  save(resultats,file = fitFile)
  
  check_fit(fitFile = fitFile,variance = resmcd$cov,outliers_labels = outlmcd,distances = distmcd)
  
  

  #save(temps_mcd,resmcd,distmcd,outlmcd,file = fitFile)
  
  temps_offline = system.time(
    {
      res0 <- offlineRobustVariance(Z,computeOutliers = TRUE,cutinit=0.6,c_m=1)
    }
  )
  
   fitFile <- paste0('Fit-Offline-withonlinequantile-machine-',j,".RData")
   
   res0ref <- offlineRobustVariance(Z[labels == 0,],computeOutliers = TRUE,cutinit=0.6,c_m=1)
   
   resultats <- list(
     variance = res0$variance,
     variance_ref = res0ref$variance,
     outliers_labels = res0$outlier_labels,
     distances = res0$distances,
     temps = temps_offline
   )
   
   save(resultats,file = fitFile)

   check_fit(fitFile = fitFile,variance = res0$variance,outliers_labels = res0$outliers_labels,distances = res0$distances)
   
   
  # 
   temps_offline = system.time(
     {
       res0 <- offlineRobustVariance_old(Z,computeOutliers = TRUE)
   
     }
   )
   
   fitFile <- paste0('Fit-Offline-without_onlinequantile-machine-',j,".RData")
  
   res0ref <- offlineRobustVariance_old(Z[labels == 0,],computeOutliers = TRUE)
   
   
   resultats <- list(
     variance = res0$variance,
     variance_ref = res0ref$variance,
     outliers_labels = res0$outlier_labels,
     distances = res0$distances,
     temps = temps_offline
   )
   
   save(resultats,file = fitFile)
   
   check_fit(fitFile = fitFile,variance = res0$variance,outliers_labels = res0$outliers_labels,distances = res0$distances)
   
  
  ###############################################Online us#########################################
  
  fitFile <- paste0('Fit-OnlineUsQuantonlinecorr-machine-',j,".RData")
  
  
  temps_online = system.time(
    {
      
      resUsOnline= onlineRobustVariance(Z,computeOutliers = TRUE,batch = 1,cutoff=.95,cutinit=0.6,nDataInit = 1e3,c_m= 1)

      
    })

  resUsOnlineref= onlineRobustVariance(Z[labels == 0,],computeOutliers = TRUE,batch = 1,cutoff=.95,cutinit=0.6,nDataInit = 1e3,c_m= 1)
  
  
  resultats <- list(
    variance = resUsOnline$variance,
    variance_ref = resUsOnlineref$variance,
    outliers_labels = resUsOnline$outlier_labels,
    distances = resUsOnline$distances,
    temps = temps_online
  )
  
  save(resultats,file = fitFile)
  
  check_fit(fitFile = fitFile,variance = resUsOnline$variance,outliers_labels = resUsOnline$outliers_labels,distances = resUsOnline$distances)
  

  
  fitFile <- paste0('Fit-OnlineUsWithoutQuantonlinecorr-machine-',j,".RData")
  
   temps_online = system.time(
     {
       
       resUsOnline= onlineRobustVariance_old(Z,batch = 1,nDataInit = 1e3,computeOutliers = TRUE)
       
     })
  # 
  # 
   
   resUsOnlineref= onlineRobustVariance_old(Z[labels == 0,],nDataInit = 1e3,computeOutliers = TRUE,batch = 1)
   
   resultats <- list(
     variance = resUsOnline$variance,
     variance_ref = resUsOnlineref$variance,
     outliers_labels = resUsOnline$outlier_labels,
     distances = resUsOnline$distances,
     temps = temps_online
   )
   
   save(resultats,file = fitFile)
   
   check_fit(fitFile = fitFile,variance = resUsOnline$variance,outliers_labels = resUsOnline$outliers_labels,distances = resUsOnline$distances)
   
  
  
  ###############################################Streaming us#########################################
  
  fitFile <- paste0('Fit-StreamingUsonlineQuantcorr-machine-', j,".RData")
  
  
  temps_streaming = system.time(
    {
      
      resStrm <- onlineRobustVariance(Z,computeOutliers = TRUE,cutoff=.95,cutinit=0.6,nDataInit = 1e3,c_m= 1,batch = batchStrm)
      
    })
  
  resStrmref <- onlineRobustVariance(Z[labels == 0,],computeOutliers = TRUE,cutoff=.95,cutinit=0.6,nDataInit = 1e3,c_m= 1,batch = batchStrm)
  
  resultats <- list(
    variance = resStrm$variance,
    variance_ref = resStrmref$variance,
    outliers_labels = resStrm$outlier_labels,
    distances = resStrm$distances,
    temps = temps_streaming
  )
  
  save(resultats,file = fitFile)
  
  check_fit(fitFile = fitFile,variance = resStrm$variance,outliers_labels = resStrm$outliers_labels,distances = resStrm$distances)
  
  
  
  # 
  fitFile <- paste0('Fit-StreamingUs_without_onlineQuantcorr-machine-', j,".RData")
  # 
   temps_streaming = system.time(
     {
       
       resStrm <- onlineRobustVariance_old(Z,computeOutliers = TRUE,nDataInit = 1e3,batch = batchStrm)
       
       
     })
   
   resStrmref <- onlineRobustVariance_old(Z[labels == 0,],nDataInit = 1e3,computeOutliers = TRUE,batch = batchStrm)
   
   
   resultats <- list(
     variance = resStrm$variance,
     variance_ref = resStrmref$variance,
     outliers_labels = resStrm$outlier_labels,
     distances = resStrm$distances,
     temps = temps_streaming
   )
   
   save(resultats,file = fitFile)
   
   check_fit(fitFile = fitFile,variance = resStrm$variance,outliers_labels = resStrm$outliers_labels,distances = resStrm$distances)

   fitFile <- paste0('Fit-OGK-machine-',j,".RData")
   
      
  outlogk = rep(0,nrow(Z))
  
  temps_ogk = system.time(
    {
      resogk = covOGK(Z, sigmamu = scaleTau2)
      for (s in 1:nrow(Z)){
      if(resogk$distances[s] > cut_off_ogk){outlogk[s] = 1}    }
    }
  )

  resogkref = covOGK(Z[labels == 0,], sigmamu = scaleTau2)
  
  resultats <- list(
    variance = resogk$cov,
    variance_ref = resogkref$cov,
    outliers_labels = outlogk,
    distances = resogk$distances,
    temps = temps_ogk
  )
  
  save(resultats,file = fitFile)
  
  check_fit(fitFile = fitFile,variance = resogk$cov,outliers_labels = outlogk,distances = resogk$distances)
  
  

  
  
}