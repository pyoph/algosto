##############Simulation parameters for our parameters##################
batch = d
cm = 2
Ninit = 1e2
k_values = seq(0,100,step = 1)

########################################################################

################Computations of the algorithms#############################


for(sim in 1:1){
  
for(sc in sc_test)
  {
  
  k = sc$k
  l = sc$l
  rho1 = sc$rho1
  

  
  

  
  for (m in seq_along(rList[1:9])){
  setwd(SimDir)
    r = rList[m]
    dataFile <- paste0('SimData-d', d, '-n', n, '-k', k, '-l', l, '-rho', rho1,'-r',r ,"-sim",sim,".RData")
    
    
    print(dataFile)
    
    
    load(dataFile)
        
    
    setwd(resAlgo)
    
    
    #######################Sample cov with online quantile correction########################################
    
    fitFile <- paste0('Fit-SampleNaiveQuantonlinecorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
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
    
    
    resultats <- list(
      variance = resNaif$Sigma,
      outliers_labels = resNaif$outliers_labels,
      distances = resNaif$distances,
      temps = temps_naif
    )
    
    save(resultats, file = fitFile)    
    
    
    load(fitFile)
    
    
    check_fit(fitFile = fitFile,variance = resNaif$Sigma,outliers_labels = resNaif$outliers_labels,distances = resNaif$distances)
    
    
    #######################Sample cov without online quantile correction########################################
    
    
    fitFile <- paste0('Fit-SampleNaivewithoutonlinequantilecorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
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
    
    resultats <- list(
      variance = resNaif$Sigma,
      outliers_labels = resNaif$outliers_labels,
      distances = resNaif$distances,
      temps = temps_naif
    )
    
    save(resultats,file = fitFile)
    
    
    load(fitFile)
    
    
    check_fit(fitFile = fitFile,variance = resNaif$Sigma,outliers_labels = resNaif$outliers_labels,distances = resNaif$distances)
    
    
    
    ###############################################Online us#########################################
    
     fitFile <- paste0('Fit-OnlineUsQuantonlinecorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
     
     
     temps_online = system.time(
       {
     
         resUsOnline= onlineRobustVariance(data$Z,batch = 1,computeOutliers = TRUE,cutoff = .95)
         
       })
     
     
     resultats <- list(
       variance = resUsOnline$variance,
       outliers_labels = resUsOnline$outlier_labels,
       distances = resUsOnline$distances,
       temps = temps_online
     )
     
     save(resultats,file = fitFile)
     
     

     
     check_fit(fitFile = fitFile,variance = resUsOnline$variance,outliers_labels = resUsOnline$outliers_labels,distances = resUsOnline$distances)
     
     
     
  fitFile <- paste0('Fit-OnlineUswithoutQuantonlinecorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
  
  temps_online = system.time(
    {
      
      resUsOnline= onlineRobustVariance_old(data$Z,batch = 1,computeOutliers = TRUE)
      
    })
  
  
  resultats <- list(
    variance = resUsOnline$variance,
    outliers_labels = resUsOnline$outlier_labels,
    distances = resUsOnline$distances,
    temps = temps_online
  )
  
  save(resultats,file = fitFile)
  
  
  load(fitFile)
  
  check_fit(fitFile = fitFile,variance = resUsOnline$variance,outliers_labels = resUsOnline$outliers_labels,distances = resUsOnline$distances)
  
  ###############################################Streaming us#########################################

  setwd(resAlgo)
  
  fitFile <- paste0('Fit-StreamingUsonlineQuantcorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
  
  
  temps_streaming = system.time(
    {
      
      resUsStreaming= onlineRobustVariance(data$Z,computeOutliers = TRUE,cutoff = .95)
      
    })
  

  
  resultats <- list(
    variance = resUsStreaming$variance,
    outliers_labels = resUsStreaming$outlier_labels,
    distances = resUsStreaming$distances,
    temps = temps_streaming
  )
  
  save(resultats,file = fitFile)
  
  load(fitFile)
  
  check_fit(fitFile = fitFile,variance = resUsStreaming$variance,outliers_labels = resUsStreaming$outliers_labels,distances = resUsStreaming$distances)
  
  

  
  
  fitFile <- paste0('Fit-StreamingUswithoutQuantonlinecorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
  
  temps_streaming = system.time(
    {
      
      resUsStreaming= onlineRobustVariance_old(data$Z,computeOutliers = TRUE)
      
    })
  
  
  
  resultats <- list(
    variance = resUsStreaming$variance,
    outliers_labels = resUsStreaming$outlier_labels,
    distances = resUsStreaming$distances,
    temps = temps_streaming
  )
  
  save(resultats,file = fitFile)
  
  
  load(fitFile)
  
  check_fit(fitFile = fitFile,variance = resUsStreaming$variance,outliers_labels = resUsStreaming$outliers_labels,distances = resUsStreaming$distances)
  
  
  #####################################################Offline Us########################################################
  
  fitFile <- paste0('Fit-OfflinewithQuantcorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
  
  
  temps_offline = system.time(
    {
      
      resOffline = offlineRobustVariance(data$Z,computeOutliers = TRUE,cutoff = .95)
      
    })
  
  
  
  resultats <- list(
    variance = resOffline$variance,
    outliers_labels = resOffline$outlier_labels,
    distances = resOffline$distances,
    temps = temps_offline
  )
  
  save(resultats,file = fitFile)
  
  load(fitFile)
  
  check_fit(fitFile = fitFile,variance = resOffline$variance,outliers_labels = resOffline$outliers_labels,distances = resOffline$distances)
  

  fitFile <- paste0('Fit-OfflineUswithoutQuantcorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
  
  temps_offline = system.time(
    {
      
      resOffline= offlineRobustVariance_old(data$Z,computeOutliers = TRUE)
      
    })
  
  resultats <- list(
    variance = resOffline$variance,
    outliers_labels = resOffline$outlier_labels,
    distances = resOffline$distances,
    temps = temps_offline
  )
  
  save(resultats,file = fitFile)
  
  
  
  check_fit(fitFile = fitFile,variance = resOffline$variance,outliers_labels = resOffline$outliers_labels,distances = resOffline$distances)
  
  
  #############################OGK#########################################################################
  
  fitFile <- paste0('Fit-OGK-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
  
  outlogk = rep(0,n)

  temps_ogk = system.time({
    resOGK = covOGK(data$Z,sigmamu = scaleTau2)})
  
  for(s in 1:n)
  {
    if(resOGK$distances[s] >qchisq(.95,df = d)) {outlogk[s] = 1}
  }
  
  
  resultats <- list(
    variance = resOGK$cov,
    outliers_labels = outlogk,
    distances = resOGK$distances,
    temps = temps_ogk
  )
  
  save(resultats,file = fitFile)
  
  
  check_fit(fitFile = fitFile,variance = resOGK$cov,outliers_labels = outlogk,distances = resOGK$distances)
  
  
  
  #############################MCD#########################################################################
  
  fitFile <- paste0('Fit-MCD-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
  
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
  
  resultats <- list(
    variance = resMcd$cov,
    outliers_labels = outlmcd,
    distances = distmcd,
    temps = temps_mcd
  )
  
  save(resultats,file = fitFile)
  
  
  
  check_fit(fitFile = fitFile,variance = resMcd$cov,outliers_labels = outlmcd,distances = distmcd)
  
  
  
  ################################Oracle###############################################################
  
  fitFile <- paste0('Fit-Oracle-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
  
    
  distoracle = rep(0,n)
  outloracle = rep(0,n)
  invSigma0 = solve(Sigma0)
  
  for (s in 1:n){
    distoracle[s] = t(data$Z[s,] - mu0)%*%invSigma0%*%(data$Z[s,] - mu0)
    if(distoracle[s] > qchisq(.95,df = d)){
      outloracle[s] = 1
      }
  }
  
  
  resultats <- list(
    variance = Sigma0,
    outliers_labels = outloracle,
    distances = distoracle,
    temps = 0
  )
  
  save(resultats,file = fitFile)

  
  check_fit(fitFile = fitFile,variance = Sigma0,outliers_labels = outloracle,distances = distoracle)
  
  
    
}}}
  




