###################Offline cumulative computation times
Ninit = 500
inv_safe <- function(S, eps = 1e-4) {
  
  S <- as.matrix(S)
  invS = diag(ncol(S))
  if (any(!is.finite(S))) {
    stop("Matrix contains NA/Inf")
  }
  
  tryCatch({
    
  invS =   solve(S)
    
  }, error = function(e) {
    
    cat("Singular matrix -> ridge correction applied\n")
    
    invS = solve(S + eps * diag(ncol(S)))
    
  })
  
  return (invS)
}

####Naive baseline method
SampleCovOnline = function(Z,nDataInit=100)
{
  nblignes = nrow(Z)
  n0 = nDataInit 
  
  # Initialization
  mean = colMeans(Z[1:n0, , drop = FALSE])
  meanOld = mean
  Sigma = cov(Z[1:n0, , drop = FALSE])
  
  M = n0*mean
  
  A = matrix(0,ncol(Z),ncol(Z))
  
  for (j in (1:n0))
  {
    A = A + (Z[j,] - mean)%*%t(Z[j,] - mean)
  }
  
  B <- NULL
  
  for (lambda in c(0, 1e-10, 1e-8, 1e-6, 1e-4)) {
    
    B <- tryCatch(
      solve(A + lambda * diag(ncol(A))),
      error = function(e) NULL
    )
    
    if (!is.null(B)) break
  }
  
  if (is.null(B))
    stop("Unable to invert covariance matrix.")
  
  # B = solve(A)
  
  invSigma = n0*B
  
  meanIter = matrix(0, nrow(Z), ncol(Z))
  SigmaIter = array(0, dim = c(nrow(Z), ncol(Z), ncol(Z)))
  
  meanIter[1:n0,] = mean
  
  distances = rep(0, nrow(Z))
  outliers_labels = rep(0, nrow(Z))
  
  meanOld = mean
  
  distances = rep(0, nrow(Z))
  outliers_labels = rep(0,nrow(Z))
  
  cutoff = qchisq(.95,df = ncol(Z))
  
  for(j in 1:n0)
  {
    distances[j] = t(Z[j,] - meanIter[j,])%*%invSigma%*%(Z[j,] - meanIter[j,])
    if(distances[j] > cutoff)
    {
      outliers_labels[j] = 1
    }
  }
  
  
  
  for (i in ((n0 + 1):(nblignes-1)))
  {
    
    M = M + Z[i,]
    Xbar = M/i
    
    U = Z[i,] - Xbar
    
    scal = 1 + t(U)%*%B%*%U
    
    if(scal !=0)
    {
      B = B - 1/scal[1]*(B%*%U)%*%(t(B%*%U))
    }
    
    invSigma = i*B
    
    mean = mean + (1.0 / (i + 1)) * (Z[i+1,] - mean)
    
    Sigma = (i - 1)/i*Sigma +
      1/(i + 1)*((Z[i+1,] - meanOld)%*%t(Z[i+1,] - meanOld))
    
    meanOld = mean
    
    meanIter[i,] = mean
    SigmaIter[i,,] = Sigma
    
    distances[i] = t(Z[i,] - meanIter[i,])%*%invSigma%*%(Z[i,] - meanIter[i,])
    
    S = distances[i]
    
    
    cutoff = qchisq(.95,df = ncol(Z))
    
    if (S > cutoff)
    {
      outliers_labels[i] = 1
    }
    
    
  }
  
  
  SigmaIter[nrow(Z),,] = Sigma
  meanIter[nrow(Z),] = mean
  
  return(list(
    mean = mean,
    Sigma = Sigma,
    meanIter = meanIter,
    SigmaIter = SigmaIter,
    distances = distances,
    outliers_labels = outliers_labels
  ))
  
}


offlineComputationTimes <- function(Z,
                                    Ninit = 5*1e2,
                                    batch = 10,
                                    machine,
                                    repRes = res_SMD) {
  
  n <- nrow(Z)
  
  # Clock times :
  # Ninit, Ninit + batch, Ninit + 2*batch, ..., n
  sizes <- seq(
    from = Ninit,
    to = n,
    by = batch
  )
  
  # Ajouter n si la dernière étape ne tombe pas exactement sur n
  if (tail(sizes, 1) != n) {
    sizes <- c(sizes, n)
  }
  
  # Temps cumulés, enregistrés observation par observation
  cum_times_MCD <- numeric(n)
  cum_times_OGK <- numeric(n)
  
  # Temps cumulés au niveau des batches
  cumulative_MCD <- 0
  cumulative_OGK <- 0
  
  #outliers labels
  
  outl_labels_mcd <- rep(0,n)
  outl_labels_ogk <- rep(0,n)
  
  # Première observation du batch courant
  start_idx <- 1
  
  for (i in seq_along(sizes)) {
    
    ni <- sizes[i]
    
    # On prend les lignes 1,...,ni
    Zi <- Z[1:ni, , drop = FALSE]
    
    # -----------------
    # MCD
    # -----------------
    t0 <- Sys.time()
    
    tryCatch({ 
      mcd <- covMcd(Zi)
      dist_mcd = mcd$mah
      cutoff = qchisq(p = .95, df = ncol(Z))
      outl_labels_mcd[1:ni] = as.numeric(dist_mcd > cutoff)
      
      
      
    }, error = function(e) {
      
      message(
        "Erreur MCD pour i = ", i,
        " (n = ", ni, ") : ",
        conditionMessage(e)
      )
      print(
        paste0(
          "ERREUR MCD - i = ", i,
          " - n = ", ni,
          " - message : ", conditionMessage(e)
        )
      )
      # En cas d'erreur, on met NA pour ce batch
      outl_labels_mcd[1:ni] <- NA_real_
    })
    current_time_MCD <- as.numeric(Sys.time() - t0, units = "secs")
    # -----------------
    # OGK
    # -----------------
    t0 <- Sys.time()
    
    tryCatch({
      
      ogk <- covOGK(
        Zi,
        sigmamu = scaleTau2
      )
      
      dist_ogk <- ogk$distances
      
      cutoff <- qchisq(
        p = .95,
        df = ncol(Z)
      )
      
      outl_labels_ogk[1:ni] <- as.numeric(
        dist_ogk > cutoff
      )
      
      
      
    }, error = function(e) {
      
      message(
        "Erreur OGK pour i = ", i,
        " (n = ", ni, ") : ",
        conditionMessage(e)
      )
      
      print(
        paste0(
          "ERREUR OGK - i = ", i,
          " - n = ", ni,
          " - message : ", conditionMessage(e)
        )
      )
    })
    current_time_OGK <- as.numeric(
      Sys.time() - t0,
      units = "secs"
    )
    
    # -----------------
    # Temps cumulés
    # -----------------
    cumulative_MCD <- cumulative_MCD + current_time_MCD
    cumulative_OGK <- cumulative_OGK + current_time_OGK
    
    # Toutes les observations correspondant à ce batch
    # reçoivent le même temps cumulé
    cum_times_MCD[start_idx:ni] <- cumulative_MCD
    cum_times_OGK[start_idx:ni] <- cumulative_OGK
    
    # Batch suivant
    start_idx <- ni + 1
  }
  
  setwd(repRes)
  
  fitFile = paste0("Fit-MCDonlineQC-","machine-",machine,".RData")
  
  resultats = list(outliers_labels = outl_labels_mcd,cum_times = cum_times_MCD)
  
  save(resultats,file = fitFile)
  
  fitFile = paste0("Fit-OGKonlineQC-","machine-",machine,".RData")
  
  resultats = list(outliers_labels = outl_labels_ogk,cum_times = cum_times_OGK)
  
  save(resultats,file = fitFile)
  
  
  
  result <- data.frame( N = 1:n, time_MCD = cum_times_MCD, time_OGK = cum_times_OGK, outlMCD = outl_labels_mcd, outlOGK = outl_labels_ogk ) 
  return(result)
}





#####################Execution of the algorithms on 3 machines of the SMD dataset###################


for(j in c(5,35,56)){
  
  setwd(smd_data_dir)
  
  data_smd_mach = paste0("data_machine-",j,".RData")
  
  load(data_smd_mach)  
  
  print(data_smd_mach)
  
  setwd(res_SMD)
  
  
  Z_clean = Z[labels == 0,]
  
  ###########Extraction moyenne et covariance sans outliers
  
  meanTrueCov = colMeans(Z_clean)
  
  Sigma_trueCov = cov(Z_clean)
  
  invSigmaTrueCov = inv_safe(Sigma_trueCov)
  
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
  

  ##########################Cov online####################
  
  temps_covonline = system.time(
    {
      resSamplecov= SampleCovOnline(Z,nDataInit = Ninit)
    })
  
  
  fitFile <- paste0(
    "Fit-SampleNaivewithoutonlinequantilecorr-", "machine-",j,".RData"
    
  )
  
  resSamplecov_ref = SampleCovOnline(Z[labels == 0,],nDataInit = Ninit)
  
  
  
  
  resultats <- list(
    variance = resSamplecov$Sigma,
    variance_ref =  resSamplecov_ref$Sigma,
    outliers_labels = resSamplecov$outliers_labels,
    distances = resSamplecov$distances,
    temps = temps_covonline
  )
  
  save(resultats,file = fitFile)
  
  
  
  
  ###########################Offline#####################
  
  
  temps_offline = system.time(
    {
      res0 <- offlineRobustVariance(Z,computeOutliers = TRUE,cutoff_method = "Chi-square")
    }
  )
  
  fitFile <- paste0(
    'Fit-OfflineUswithoutQuantcorr-',"machine-",j,".RData"
  )
  
  res0ref <- offlineRobustVariance(Z[labels == 0,],computeOutliers = TRUE,cutoff_method = "Chi-square")
  
  resultats <- list(
    variance = res0$variance,
    variance_ref = res0ref$variance,
    outliers_labels = res0$outlier_labels,
    distances = res0$distances,
    temps = temps_offline
  )
  
  save(resultats,file = fitFile)
  
  
  
  
  ###############################################Streaming us#########################################
  
  
  fitFile <- paste0('Fit-StreamingUswithoutQuantonlinecorr-',"machine-",j,".RData"
  )
  
  temps_streaming = system.time(
    {
      
      resUsStreaming <- onlineRobustVariance(
        Z,
        computeOutliers = TRUE,
        cutoff_method="Chi-square",nDataInit =  Ninit,batch = batch_size
      )      
    })
  
  resStrmref <- onlineRobustVariance(Z[labels == 0,],computeOutliers = TRUE,
                                     cutoff_method="Chi-square",nDataInit =  Ninit
                                     ,batch = batch_size)
  
  resultats <- list(
    variance = resUsStreaming$variance,
    variance_ref = resStrmref$variance,
    outliers_labels = resUsStreaming$outlier_labels,
    distances = resUsStreaming$distances,
    temps = temps_streaming,
    cum_times = resUsStreaming$time_cumulative
  )
  
  save(resultats,file = fitFile)
  
  
  

  
  ###############################################Online us#########################################
  
  
  fitFile <- paste0(
    'Fit-OnlineUswithoutQuantonlinecorr-',"machine-",j,".RData"
    
  )
  
  temps_online = system.time(
    {
      
      
      resUsOnline <- onlineRobustVariance(
        Z,
        batch = 1,
        computeOutliers = TRUE,cutoff_method="Chi-square",
        ,nDataInit =  Ninit
      )    
    })
  
  resUsOnlineref= onlineRobustVariance(Z[labels == 0,],computeOutliers = TRUE,batch = 1,
                                      cutoff_method="Chi-square",
                                       ,nDataInit =  Ninit)
  
  resultats <- list(
    variance = resUsOnline$variance,
    variance_ref = resUsOnlineref$variance,
    outliers_labels = resUsOnline$outlier_labels,
    distances = resUsOnline$distances,
    temps = temps_online,
    cum_times = resUsOnline$time_cumulative
  )
  
  save(resultats,file = fitFile)
  
  
  load(fitFile)
  

  
  
  #######OGK#####################
  
  fitFile <- paste0('Fit-OGK-machine-',j,".RData")
  
  
  outlogk = rep(0,nrow(Z))
  
  temps_ogk = system.time(
    {
      resogk = covOGK(Z, sigmamu = scaleTau2)
      for (s in 1:nrow(Z)){
        if(resogk$distances[s] > qchisq(.95,df = ncol(Z))){outlogk[s] = 1}    }
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
  
  ####################################MCD####################
  
  
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
        
        if (distmcd[s] > qchisq(.95,df = ncol(Z))){
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
  
  #############Offline cumulative computation times###############
  
   offlineComputationTimes(
    Z,
    Ninit = Ninit,
    batch = 10, machine = j
  )
  
  
}
