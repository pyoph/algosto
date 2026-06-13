##############Simulation parameters for our parameters##################
batch = d
cm = 2
Ninit = 1e2
########################################################################


################Computations of the algorithms#############################
for(sim in 1:100){
for(sc in scenarios){
  
  k = sc$k
  l = sc$l
  rho1 = sc$rho1
  

  
  

  
  for (m in seq_along(rList[1:9])){
  setwd(SimDir)
    r = rList[m]
    dataFile <- paste0('SimData-d', d, '-n', n, '-k', k, '-l', l, '-rho', rho1,'-r',r ,".RData")
    
    
    print(dataFile)
    
    
    load(dataFile)
        
    
    setwd(resAlgo)
    
    
    #######################Sample cov with online quantile correction########################################
    
    fitFile <- paste0('FitSampleNaiveQuantonlinecorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,"-cutoff",.95,"Ninit-",Ninit,"-cm",cm,'-sim', sim,".RData")
    temps_naif <- system.time(
      resNaif <- tryCatch(
        {
          SampleCovOnline(data$Z, quantcutoff = TRUE,nDataInit = Ninit,cutoffquant = .95)
        },
        error = function(e) {
          message("SampleCovOnline failed: ", e$message)
          NULL
        }
      )
      
      
      
    )
    
    fitNaif = resNaif
    
    save(resNaif,temps_naif,file = fitFile)
    
    
    #######################Sample cov without online quantile correction########################################
    
    
    fitFile <- paste0('FitSampleNaivewithoutonlinequantilecorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,"-Ninit",Ninit,'-sim', sim,".RData")
    temps_naif <- system.time(
      resNaif <- tryCatch(
        {
          SampleCovOnline(data$Z, quantcutoff = FALSE,nDataInit = Ninit)
        },
        error = function(e) {
          message("SampleCovOnline failed: ", e$message)
          NULL
        }
      )
      
      
      
    )
    
    fitNaif = resNaif
    
    save(resNaif,temps_naif,file = fitFile)
    
    ###############################################Online us#########################################
    
     fitFile <- paste0('FitOnlineUsQuantonlinecorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,"-cutoff",.95,"Ninit-",Ninit,"-batch",1,"-cm",cm,'-sim', sim,".RData")
     
     
     temps_online = system.time(
       {
     
         resUsOnline= onlineRobustVariance(data$Z,batch = 1,computeOutliers = TRUE,cutoff = .95)
         
       })
     
     fitUsOnline = resUsOnline
     
     save(resUsOnline,temps_online,file = fitFile)
     
     
     
     
  fitFile <- paste0('FitOnlineUswithoutQuantonlinecorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,"-cutoff",.95,"Ninit-",1e2,"-batch",1,'-sim', sim,".RData")
  
  temps_online = system.time(
    {
      
      resUsOnline= onlineRobustVariance_old(data$Z,batch = 1,computeOutliers = TRUE)
      
    })
  

  save(resUsOnline,temps_online,file = fitFile)

  
  ###############################################Streaming us#########################################

  fitFile <- paste0('FitStreamingUsonlineQuantcorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,"-cutoff",.95,"Ninit-",Ninit,"-batch",batch,"-cm",cm,'-sim', sim,".RData")
  
  
  temps_streaming = system.time(
    {
      
      resUsStreaming= onlineRobustVariance(data$Z,computeOutliers = TRUE,cutoff = .95)
      
    })
  

  save(resUsStreaming,temps_streaming,file = fitFile)
  

  
  
  fitFile <- paste0('FitStreamingUswithoutQuantonlinecorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,"Ninit-",1e2,"-batch",batch,'-sim', sim,".RData")
  
  temps_streaming = system.time(
    {
      
      resUsStreaming= onlineRobustVariance_old(data$Z,computeOutliers = TRUE)
      
    })
  

  save(resUsStreaming,temps_streaming,file = fitFile)
  
  
  #####################################################Offline Us########################################################
  
  fitFile <- paste0('FitOfflinewithQuantcorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,"-cutoff",.95,"-cm",cm,'-sim', sim,".RData")
  
  
  temps_offline = system.time(
    {
      
      resOffline = offlineRobustVariance(data$Z,computeOutliers = TRUE,cutoff = .95)
      
    })
  

  save(resOffline,temps_offline,file = fitFile)
  
  
  
  
  fitFile <- paste0('FitOfflineUswithoutQuantcorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,"-cutoff",.95,"-cm",cm,'-sim', sim,".RData")
  
  temps_offline = system.time(
    {
      
      resOffline= offlineRobustVariance_old(data$Z,computeOutliers = TRUE)
      
    })
  
  fitUsOffline = resUsOnline
  
  save(resOffline,temps_offline,file = fitFile)
  
  #############################OGK#########################################################################
  
  fitFile <- paste0('FitOGK-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,"-cutoff",.95,'-sim', sim,".RData")
  
  outlogk = rep(0,n)

  temps_ogk = system.time({
    resOGK = covOGK(data$Z,sigmamu = scaleTau2)})
  
  for(s in 1:n)
  {
    if(resOGK$distances[s] >qchisq(.95,df = d)) {outlogk[s] = 1}
  }
  
  save(resOGK,outlogk,temps_ogk,file = fitFile)
  
  #############################MCD#########################################################################
  
  fitFile <- paste0('FitMCD-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,"-cutoff",.95,'-sim', sim,".RData")
  
  distmcd = rep(0,n)
  outlmcd = rep(0,n)
  temps_mcd = system.time({
    resMcd = covMcd(data$Z)
    invSigmaMCD = solve(resMcd$cov) 
    for(s in (1:n))
    {
      distmcd[s] = t(data$Z[s,] - resMcd$center)%*%invSigmaMCD%*%(data$Z[s,] - resMcd$center)
      if (distmcd[s] > qchisq(.95,df = d)){outlmcd[s] = 1}
    }
  }
  )
  save(resMcd,temps_mcd,distmcd,outlmcd,file = fitFile)
  
  ################################Oracle###############################################################
  
  fitFile <- paste0('FitOracle-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,"-cutoff",.95,'-sim', sim,".RData")
  
    
  distoracle = rep(0,n)
  outloracle = rep(0,n)
  invSigma0 = solve(Sigma0)
  
  for (s in 1:n){
    distoracle[s] = t(data$Z[s,] - mu0)%*%invSigma0%*%(data$Z[s,] - mu0)
    if(distoracle[s] > qchisq(.95,df = d)){
      outloracle[s] = 1
      }
  }
  
  save(distoracle,outloracle,file = fitFile)
  
    
}
  
}}




