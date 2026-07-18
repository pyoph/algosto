##############Simulation parameters for our parameters##################
batch = d
cm = 2
Ninit = 500

###################################################

save_add_method <- function(dist, fitFile, type) {
  
  cutoff <- calcule_cutoff(dist, type = type)
  
  
  resultats <- list(
    outliers_labels = as.integer(dist > cutoff)
  )
  
  save(resultats, file = fitFile)
}



################################Chargement des scénarios#################################

setwd("~/algosto/")

load("scenarios.RData")

scenarios = c(scenarios_1_param,scenarios_2_param)

################Computations of the algorithms#############################


for(sim in 1:2){
  
  for(sc in scen_strong_conc)
  {
    
    k = sc$k
    l = sc$l
    rho1 = sc$rho1
    
    
    
    for (m in seq_along(rList)){
      setwd(SimDir)
      r = rList[m]
      dataFile <- paste0('SimData-d', d, '-n', n, '-k', k, '-l', l, '-rho', rho1,'-r',r ,"-sim",sim,".RData")
      
      
      print(dataFile)
      
      load(dataFile)
      
      
      setwd(resAlgo)
      
      ################################Sample naive########################################################################################
      
      fitFile <- paste0(
        "Fit-SampleNaiveQuantonlinecorr-d", d,
        "-n", n,
        "-k", k,
        "-l", l,
        "-rho", rho1,
        "-r", r,
        "-sim", sim,
        ".RData"
      )
      
      
      temps_naif <- system.time(
        resNaif <- tryCatch(
          {
            SampleCovOnline(
              data$Z,
              quantcutoff = TRUE,
              nDataInit = Ninit,
              cutoffquant = .95
            )
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
      
      dist = resultats$distances
      
      fitFile <- paste0(
        "Fit-SampleNaivewithoutonlinequantilecorr-d", d,
        "-n", n,
        "-k", k,
        "-l", l,
        "-rho", rho1,
        "-r", r,
        "-sim", sim,
        ".RData"
      )
      
      
      # if(!file.exists(fitFile)){
      #   temps_naif <- system.time(
      
      save(dist,file = fitFile)
      # }
      # else{
      # print("File exists")
      # load(fitFile)
      # }
      # 
      
      
      # if(!file.exists(fitFile)){
      #   temps_naif <- system.time(
      #   resNaif <- tryCatch(
      #     {
      #       SampleCovOnline(data$Z, quantcutoff = FALSE,nDataInit = Ninit)
      #     },
      #     error = function(e) {
      #       message("SampleCovOnline failed: ", e$message)
      #       NULL
      #     }
      #   )
      #   
      #   
      #   
      # )
      # 
      # resultats <- list(
      #   #variance = resNaif$Sigma,
      #   outliers_labels = resNaif$outliers_labels,
      #   distances = resNaif$distances,
      #   temps = temps_naif
      # )
      # 
      # save(resultats,file = fitFile)
      # }
      # else{
      # print("File exists")
      # load(fitFile)
      # }
      # 
      
      fitFile <- paste0('Fit-SampleRaw-d', d,
                        '-n', n,
                        '-k', k,
                        '-l', l,
                        '-rho', rho1,
                        '-r', r,
                        '-sim', sim,
                        ".RData")
      
      save_add_method(dist, fitFile, "raw")
      
      # if (!file.exists(fitFile)) {
      #   save_add_method(resultats, fitFile, "raw")
      #       } 
      # else {
      #   
      #   print("File exists")
      #   load(fitFile)
      #   
      # }
      
      ###############################################Online us#########################################
      
      fitFile <- paste0('Fit-OnlineUsQuantonlinecorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
      
      temps_online = system.time(
        {
          
          resUsOnline= onlineRobustVariance(data$Z,batch = 1,computeOutliers = TRUE,cutoff = .95,cutoff_method = "quantile",nDataInit =  Ninit)
          
        })
      
      
      resultats <- list(
        variance = resUsOnline$variance,
        outliers_labels = resUsOnline$outlier_labels,
        distances = resUsOnline$distances,
        temps = temps_online
      )
      
      save(resultats,file = fitFile)
      
      load(fitFile)
      
      dist = resultats$distances
      
      
      fitFile <- paste0(
        'Fit-OnlineUswithoutQuantonlinecorr-d', d,
        '-n', n,
        '-k', k,
        '-l', l,
        '-rho', rho1,
        '-r', r,
        '-sim', sim,
        ".RData"
      )
      
      
      
      temps_online <- system.time({
        
        resUsOnline <- onlineRobustVariance(
          data$Z,
          batch = 1,
          computeOutliers = TRUE,cutoff_method="Chi-square",
          ,nDataInit =  Ninit
        )
        
      })
      
      resultats <- list(
        #variance = resUsOnline$variance,
        outliers_labels = resUsOnline$outlier_labels
        #distances = resUsOnline$distances,
        #temps = temps_online
      )
      
      save(dist, file = fitFile)
      
      
      
      fitFile <- paste0('Fit-OnlRaw-d', d,
                        '-n', n,
                        '-k', k,
                        '-l', l,
                        '-rho', rho1,
                        '-r', r,
                        '-sim', sim,
                        ".RData")
      
      
      save_add_method(dist, fitFile, "raw")
      # if (!file.exists(fitFile)) {
      #   
      #   save_add_method(resultats, fitFile, "raw")
      #   
      # } else {
      #   
      #   print("File exists")
      #   load(fitFile)
      #   
      # }
      # 
      
      
      ###############################################Streaming us#########################################
      
      setwd(resAlgo)
      
      fitFile <- paste0('Fit-StreamingUsonlineQuantcorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
      
      
      temps_streaming = system.time(
        {
          
          resUsStreaming= onlineRobustVariance(data$Z,computeOutliers = TRUE,cutoff = .95,cutoff_method = "quantile",nDataInit =  Ninit)
          
        })
      
      
      
      resultats <- list(
        variance = resUsStreaming$variance,
        outliers_labels = resUsStreaming$outlier_labels,
        distances = resUsStreaming$distances,
        temps = temps_streaming
      )
      
      save(resultats,file = fitFile)
      
      
      dist = resultats$distances
      
      fitFile <- paste0('Fit-StreamingUswithoutQuantonlinecorr-d', d,
                        '-n', n,
                        '-k', k,
                        '-l', l,
                        '-rho', rho1,
                        '-r', r,
                        '-sim', sim,
                        ".RData")
      
      
      temps_streaming <- system.time({
        
        resUsStreaming <- onlineRobustVariance(
          data$Z,
          computeOutliers = TRUE,
          cutoff_method="Chi-square",nDataInit =  Ninit
        )
        
      })
      
      resultats <- list(
        #variance = resUsStreaming$variance,
        outliers_labels = resUsStreaming$outlier_labels
        #distances = resUsStreaming$distances,
        #temps = temps_streaming
      )
      
      save(dist, file = fitFile)
      
      
      #print("File exists")
      #load(fitFile)
      
      
      fitFile <- paste0('Fit-StrmRaw-d', d,
                        '-n', n,
                        '-k', k,
                        '-l', l,
                        '-rho', rho1,
                        '-r', r,
                        '-sim', sim,
                        ".RData")
      
      
      save_add_method(dist, fitFile, "raw")
      
      # if (!file.exists(fitFile)) {
      #   
      #   save_add_method(resultats, fitFile, "raw")
      #   
      # } else {
      #   
      #   print("File exists")
      #   load(fitFile)
      #   
      # }
      # 
      
      #####################################################Offline Us########################################################
      
      fitFile <- paste0('Fit-OfflinewithQuantcorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
      
      temps_offline = system.time(
        {
          
          resOffline = offlineRobustVariance(data$Z,computeOutliers = TRUE,cutoff = .95,cutoff_method = "quantile")
          
        })
      
      
      
      resultats <- list(
        variance = resOffline$variance,
        outliers_labels = resOffline$outlier_labels,
        distances = resOffline$distances,
        temps = temps_offline
      )
      
      save(resultats,file = fitFile)
      
      
      dist = resultats$distances
      
      fitFile <- paste0(
        'Fit-OfflineUswithoutQuantcorr-d', d,
        '-n', n,
        '-k', k,
        '-l', l,
        '-rho', rho1,
        '-r', r,
        '-sim', sim,
        ".RData"
      )
      
      temps_offline <- system.time({
        
        resOffline <- offlineRobustVariance(
          data$Z,
          computeOutliers = TRUE,cutoff_method = "Chi-square"
        )
        
      })
      
      resultats <- list(
        #variance = resOffline$variance,
        outliers_labels = resOffline$outlier_labels
        #distances = resOffline$distances,
        #temps = temps_offline
      )
      
      save(resultats, file = fitFile)
      
      
      
      
      fitFile <- paste0(
        'Fit-OfflRaw-d', d,
        '-n', n,
        '-k', k,
        '-l', l,
        '-rho', rho1,
        '-r', r,
        '-sim', sim,
        ".RData"
      )
      
      
      save_add_method(dist, fitFile, "raw")
      
      # 
      # if (!file.exists(fitFile)) {
      #   
      #   save_add_method(resultats, fitFile, "raw")
      #   
      # } else {
      #   
      #   print("File exists")
      #   load(fitFile)
      #   
      # }
      # 
      # 
      
      #############################OGK#########################################################################
      fitFile <- paste0(
        'Fit-OGK-d', d,
        '-n', n,
        '-k', k,
        '-l', l,
        '-rho', rho1,
        '-r', r,
        '-sim', sim,
        ".RData"
      )
      
      temps_ogk <- system.time({
        
        resOGK <- covOGK(data$Z, sigmamu = scaleTau2)
        
      })
      
      cutoff <- qchisq(0.95, df = d)
      
      outlogk <- as.integer(resOGK$distances > cutoff)
      
      resultats <- list(
        variance = resOGK$cov,
        outliers_labels = outlogk,
        distances = resOGK$distances,
        temps = temps_ogk
      )
      
      save(resultats, file = fitFile)
      
      
      dist = resultats$distances
      
      
      fitFile <- paste0(
        'Fit-OGKRD-d', d,
        '-n', n,
        '-k', k,
        '-l', l,
        '-rho', rho1,
        '-r', r,
        '-sim', sim,
        ".RData"
      )
      
      save_add_method(
        dist , fitFile, "rescale_dist")
      # if (!file.exists(fitFile)) {
      #   
      #   save_add_method(resultats, fitFile, "rescale_dist")
      #   
      # } else {
      #   
      #   print("File exists")
      #   load(fitFile)
      #   
      # }
      # 
      
      fitFile <- paste0(
        'Fit-OGKQC-d', d,
        '-n', n,
        '-k', k,
        '-l', l,
        '-rho', rho1,
        '-r', r,
        '-sim', sim,
        ".RData"
      )
      
      
      save_add_method(
        dist , fitFile, "quantcorr")
      # if (!file.exists(fitFile)) {
      #   
      #   save_add_method(resultats, fitFile, "quantcorr")
      #   
      # } else {
      #   
      #   print("File exists")
      #   load(fitFile)
      #   
      # }
      # 
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
      dist = resultats$distances
      
      fitFile <- paste0(
        'Fit-MCDRD-d', d,
        '-n', n,
        '-k', k,
        '-l', l,
        '-rho', rho1,
        '-r', r,
        '-sim', sim,
        ".RData"
      )
      
      save_add_method(dist, fitFile, "rescale_dist")
      
      # 
      # if (!file.exists(fitFile)) {
      #   
      #   save_add_method(resultats, fitFile, "rescale_dist")
      #   
      # } else {
      #   
      #   print("File exists")
      #   load(fitFile)
      #   
      # }
      # 
      # 
      
      
      
      
      fitFile <- paste0(
        'Fit-MCDQC-d', d,
        '-n', n,
        '-k', k,
        '-l', l,
        '-rho', rho1,
        '-r', r,
        '-sim', sim,
        ".RData"
      )
      
      save_add_method(dist, fitFile, "quantcorr")
      
      # 
      # if (!file.exists(fitFile)) {
      #   
      #   save_add_method(resultats, fitFile, "quantcorr")
      #   
      # } else {
      #   
      #   print("File exists")
      #   load(fitFile)
      #   
      # }
      # 
      
      
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
        #variance = Sigma0,
        outliers_labels = outloracle,
        distances = distoracle
        #temps = 0
      )
      
      save(resultats,file = fitFile)
      
      
      fitFile <- paste0(
        'Fit-OracleRD-d', d,
        '-n', n,
        '-k', k,
        '-l', l,
        '-rho', rho1,
        '-r', r,
        '-sim', sim,
        ".RData"
      )
      
      save_add_method(distoracle, fitFile, "rescale_dist")
      
      # 
      # if (!file.exists(fitFile)) {
      #   
      #   save_add_method(resultats, fitFile, "rescale_dist")
      #   
      # } else {
      #   
      #   print("File exists")
      #   load(fitFile)
      #   
      # }
      
      
      
      
      
      
      fitFile <- paste0(
        'Fit-OracleQC-d', d,
        '-n', n,
        '-k', k,
        '-l', l,
        '-rho', rho1,
        '-r', r,
        '-sim', sim,
        ".RData"
      )
      
      
      save_add_method(distoracle, fitFile, "quantcorr")
      
      #   
      # if (!file.exists(fitFile)) {
      #   
      #   save_add_method(resultats, fitFile, "quantcorr")
      #   
      # } else {
      #   
      #   print("File exists")
      #   load(fitFile)
      #   
      # }
      # 
      # 
      # 
      
    }}
}


