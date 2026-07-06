##############Simulation parameters for our parameters##################
batch = d
cm = 2
Ninit = 1e2

################################Chargement des scénarios#################################

setwd("~/algosto/")

load("scenarios.RData")

scenarios = c(scenarios_1_param,scenarios_2_param)

################Computations of the algorithms#############################


for(sim in 1:20){
  
for(sc in scen_strong_conc)
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
    
    if(!file.exists(fitFile)){
    
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
    }
    
    else{load(fitFile)}
    
    

    
    #######################Sample cov without online quantile correction########################################
    
    
    fitFile <- paste0('Fit-SampleNaivewithoutonlinequantilecorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
    
    if(!file.exists(fitFile)){
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
    }
    else{
    print("File exists")
    load(fitFile)
    }
    
    
    
    ###############################################Online us#########################################
    
     fitFile <- paste0('Fit-OnlineUsQuantonlinecorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
     
      if(!file.exists(fitFile)){  
     temps_online = system.time(
       {
     
         resUsOnline= onlineRobustVariance(data$Z,batch = 1,computeOutliers = TRUE,cutoff = .95,eps_lambda = .01)
         
       })
     
     
     resultats <- list(
       variance = resUsOnline$variance,
       outliers_labels = resUsOnline$outlier_labels,
       distances = resUsOnline$distances,
       temps = temps_online
     )
     
     save(resultats,file = fitFile)
      } else {
        print("File exists")
        
        load(fitFile)}
     

     
     
     
  fitFile <- paste0('Fit-OnlineUswithoutQuantonlinecorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
  
  
  if(!file.exists(fitFile)){
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
  }
  
  else {
    print("File exists")
    
    load(fitFile)}
  
  ###############################################Streaming us#########################################

  setwd(resAlgo)
  
  fitFile <- paste0('Fit-StreamingUsonlineQuantcorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
  
  
  if(!file.exists(fitFile)){
  temps_streaming = system.time(
    {
      
      resUsStreaming= onlineRobustVariance(data$Z,computeOutliers = TRUE,cutoff = .95,eps_lambda = .01)
      
    })
  

  
  resultats <- list(
    variance = resUsStreaming$variance,
    outliers_labels = resUsStreaming$outlier_labels,
    distances = resUsStreaming$distances,
    temps = temps_streaming
  )
  
  save(resultats,file = fitFile)
  }
  else {
    print("File exists")
    
    load(fitFile)}
  
  

  
  
  fitFile <- paste0('Fit-StreamingUswithoutQuantonlinecorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
  
  if(!file.exists(fitFile)){
  
  temps_streaming = system.time(
    {
      
      resUsStreaming= onlineRobustVariance_old(data$Z,computeOutliers = TRUE,eps_lambda = .01)
      
    })
  
  
  
  resultats <- list(
    variance = resUsStreaming$variance,
    outliers_labels = resUsStreaming$outlier_labels,
    distances = resUsStreaming$distances,
    temps = temps_streaming
  )
  
  save(resultats,file = fitFile)
  }
  
  else{
    print("File exists")
    
    load(fitFile)}
  
  
  #####################################################Offline Us########################################################
  
  fitFile <- paste0('Fit-OfflinewithQuantcorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
  
  if(!file.exists(fitFile)){
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
  
  save(resultats,file = fitFile)}
  else {load(fitFile)}
  
  fitFile <- paste0('Fit-OfflineUswithoutQuantcorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
  
  
  if(!file.exists(fitFile)){
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
  }
  else{
    print("File exists")
    
    load(fitFile)}
  
  
  
  #############################OGK#########################################################################
  
  fitFile <- paste0('Fit-OGK-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
  
  
  if(!file.exists(fitFile)){
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
  }
  else{
    print("File exists")
    
    load(fitFile)}
  
  #############################MCD#########################################################################
  
  fitFile <- paste0('Fit-MCD-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
  
  if(!file.exists(fitFile)){
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
  } else{
    print("File exists")
    
    load(fitFile)}
  
  
  
  ################################Oracle###############################################################
  
  fitFile <- paste0('Fit-Oracle-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
  
  if(!file.exists(fitFile)){  
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

  } else{
    print("File exists")
    
    load(fitFile)}
  
    
}}}
  