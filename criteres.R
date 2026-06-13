########################Computation of the evaluation criteria#######################################




for(sc in scenarios){
  
  k = sc$k
  l = sc$l
  rho1 = sc$rho1
  
  for (m in seq_along(rList[1:9])){
  
  ###############################Sample covonline#################################
  
    r = rList[m]
  
  setwd(SimDir)
  
  dataFile <- paste0('SimData-d', d, '-n', n, '-k', k, '-l', l, '-rho', rho1,'-r',r ,".RData")
  
  load(dataFile)
  
  setwd(resAlgo)
  
  
  fitFile <- paste0('FitSampleNaiveQuantonlinecorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,"-cutoff",.95,"Ninit-",Ninit,"-cm",cm,'-sim', sim,".RData")
  print(fitFile)
  
  
  load(fitFile)
  
  setwd(criteres)
  
  erreurFrob <- norm(resNaif$Sigma - Sigma0, "F")
  
  print(paste0("erreur norme de Frobenius naive ", erreurFrob))
  
  t <- table(data$labelsVrais, resNaif$outliers_labels)
  
  faux_positifs <- t[1, 2]
  
  if (r != 0) {
    faux_negatifs <- t[2, 1]
  } else {
    faux_negatifs <- 0
  }
  
  print(paste0("faux positifs naive ", faux_positifs))
  print(paste0("faux négatifs naive ", faux_negatifs))
  
  if(r != 0){
  
  # ARI
  ari <- adjustedRandIndex(
    data$labelsVrais,
    resNaif$outliers_labels
  )
  
  print(paste0("ARI naive ", round(ari, 4)))
  
  # AUC
  y_true <- as.numeric(data$labelsVrais == unique(sort(data$labelsVrais))[2])
  
  auc_val <- as.numeric(
    auc(
      response = y_true,
      predictor = as.numeric(resNaif$distances)
    )
  )
  
  print(paste0("AUC naive ", round(auc_val, 4)))
  }
  
  
  critFile <- paste0(
    'CritSampleNaiveQuantonlinecorr-d', d,
    '-n', n,
    '-k', k,
    '-l', l,
    '-rho', rho1,
    '-r', r,
    "-cutoff", .95,
    "Ninit-", Ninit,
    "-cm", cm,
    '-sim', sim,
    ".RData"
  )
  if(r!= 0){
  save(
    erreurFrob,
    faux_negatifs,
    faux_positifs,
    
    file = critFile
  )
  
  }
  else {save(
    erreurFrob,
    faux_negatifs,
    faux_positifs,
    ari,
    auc_val,
    file = critFile
  )
  }
  setwd(resAlgo)
  
  fitFile <- paste0('FitSampleNaivewithoutonlinequantilecorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,"-Ninit",Ninit,'-sim', sim,".RData")
  
  load(file = fitFile)
 
  
  erreurFrob <- norm(resNaif$Sigma - Sigma0, "F")
  
  print(paste0("erreur norme de Frobenius naive ", erreurFrob))
  
  t <- table(data$labelsVrais, resNaif$outliers_labels)
  
  faux_positifs <- t[1, 2]
  
  if (r != 0) {
    faux_negatifs <- t[2, 1]
  } else {
    faux_negatifs <- 0
  }
  
  print(paste0("faux positifs naive ", faux_positifs))
  print(paste0("faux négatifs naive ", faux_negatifs))
  if( r!= 0){
  # ARI
  ari <- adjustedRandIndex(
    data$labelsVrais,
    resNaif$outliers_labels
  )
  
  print(paste0("ARI naive ", round(ari, 4)))
  
  # AUC
  y_true <- as.numeric(data$labelsVrais == unique(sort(data$labelsVrais))[2])
  
  auc_val <- as.numeric(
    auc(
      response = y_true,
      predictor = as.numeric(resNaif$distances)
    )
  )
  
  print(paste0("AUC naive ", round(auc_val, 4)))
  }
  setwd(criteres)
  critFile <- paste0(
    'CritSampleNaivewithoutonlinequantilecorr-d', d,
    '-n', n,
    '-k', k,
    '-l', l,
    '-rho', rho1,
    '-r', r,
    "-cutoff", .95,
    "Ninit-", Ninit,
    ".RData"
  )
  
  if(r != 0){ 
  
  save(
    erreurFrob,
    faux_negatifs,
    faux_positifs,
    ari,
    auc_val,
    file = critFile
  )
  }
  else{save(
    erreurFrob,
    faux_negatifs,
    faux_positifs,
    ari,
    auc_val,
    file = critFile
  )}
  
  ############################Online us#########################################
  
  
  setwd(resAlgo)
  
  
 fitFile =  paste0('FitOnlineUsQuantonlinecorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,"-cutoff",.95,"Ninit-",Ninit,"-batch",1,"-cm",cm,'-sim', sim,".RData")
  
  print(fitFile)
  
  
  load(fitFile)
  
  
  erreurFrob <- norm(resUsOnline$variance - Sigma0, "F")
  
  print(paste0("erreur norme de Frobenius online Us ", erreurFrob))
  
  t <- table(data$labelsVrais, resUsOnline$outlier_labels)
  
  faux_positifs <- t[1, 2]
  
  if (r != 0) {
    faux_negatifs <- t[2, 1]
  } else {
    faux_negatifs <- 0
  }
  
  print(paste0("faux positifs online with quantile corr ", faux_positifs))
  print(paste0("faux négatifs online with quantile corr ", faux_negatifs))
  if(r != 0){
  # ARI
  ari <- adjustedRandIndex(
    data$labelsVrais,
    resUsOnline$outlier_labels
  )
  
  print(paste0("ARI online us ", round(ari, 4)))
  
  # AUC
  y_true <- as.numeric(data$labelsVrais == unique(sort(data$labelsVrais))[2])
  
  auc_val <- as.numeric(
    auc(
      response = y_true,
      predictor = as.numeric(resUsOnline$distances)
    )
  )
  
  print(paste0("AUC online ", round(auc_val, 4)))
  }
  setwd(criteres)
  critFile <- paste0(
    'CritOnlineUsonlinecorr-d', d,
    '-n', n,
    '-k', k,
    '-l', l,
    '-rho', rho1,
    '-r', r,
    "-cutoff", .95,
    "Ninit-", Ninit,
    "-batch",1,
    "-cm", cm,
    '-sim', sim,
    ".RData"
  )
  if( r!= 0){
  save(
    erreurFrob,
    faux_negatifs,
    faux_positifs,
    ari,
    auc_val,
    file = critFile
  )
  
  }else {
    save(
      erreurFrob,
      faux_negatifs,
      faux_positifs,
        file = critFile
    )
    
  }
  
  
  setwd(resAlgo)
  
  fitFile <- paste0('FitOnlineUswithoutQuantonlinecorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,"-cutoff",.95,"Ninit-",1e2,"-batch",1,'-sim', sim,".RData")
  
  load(fitFile)
  
  print(fitFile)
  
  
  
  erreurFrob <- norm(resUsOnline$variance - Sigma0, "F")
  
  print(paste0("erreur norme de Frobenius online Us ", erreurFrob))
  
  t <- table(data$labelsVrais, resUsOnline$outlier_labels)
  
  faux_positifs <- t[1, 2]
  
  if (r != 0) {
    faux_negatifs <- t[2, 1]
  } else {
    faux_negatifs <- 0
  }
  
  print(paste0("faux positifs online with quantile corr ", faux_positifs))
  print(paste0("faux négatifs online with quantile corr ", faux_negatifs))
  if(r != 0){
    # ARI
    ari <- adjustedRandIndex(
      data$labelsVrais,
      resUsOnline$outlier_labels
    )
    
    print(paste0("ARI online us ", round(ari, 4)))
    
    # AUC
    y_true <- as.numeric(data$labelsVrais == unique(sort(data$labelsVrais))[2])
    
    auc_val <- as.numeric(
      auc(
        response = y_true,
        predictor = as.numeric(resUsOnline$distances)
      )
    )
    
    print(paste0("AUC online ", round(auc_val, 4)))
  }
  
  setwd(criteres)
  
  critFile <- paste0(
    'CritOnlineUswithoutQuantonlinecorr-d', d,
    '-n', n,
    '-k', k,
    '-l', l,
    '-rho', rho1,
    '-r', r,
    "-cutoff", .95,
    "Ninit-", Ninit,
    "-batch",1,
    '-sim', sim,
    ".RData"
  )  

  if( r!= 0){
    save(
      erreurFrob,
      faux_negatifs,
      faux_positifs,
      ari,
      auc_val,
      file = critFile
    )
    
  } else {
    save(
      erreurFrob,
      faux_negatifs,
      faux_positifs,
      file = critFile
    )
    
  }
  
  
  
  
  ############################Streaming us#########################################
  
  fitFile <- paste0('FitStreamingUsonlineQuantcorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,"-cutoff",.95,"Ninit-",Ninit,"-batch",batch,"-cm",cm,'-sim', sim,".RData")
  
  setwd(resAlgo)
  
  print(fitFile)
  
  load(fitFile)
  
  erreurFrob <- norm(resUsStreaming$variance - Sigma0, "F")
  
  print(paste0("erreur norme de Frobenius streaming Us ", erreurFrob))
  
  t <- table(data$labelsVrais, resUsStreaming$outlier_labels)
  
  faux_positifs <- t[1, 2]
  
  if (r != 0) {
    faux_negatifs <- t[2, 1]
  } else {
    faux_negatifs <- 0
  }
  
  print(paste0("faux positifs streaming with quantile corr ", faux_positifs))
  print(paste0("faux négatifs streaming with quantile corr ", faux_negatifs))
  if(r != 0){
    # ARI
    ari <- adjustedRandIndex(
      data$labelsVrais,
      resUsStreaming$outlier_labels
    )
    
    print(paste0("ARI online us ", round(ari, 4)))
    
    # AUC
    y_true <- as.numeric(data$labelsVrais == unique(sort(data$labelsVrais))[2])
    
    auc_val <- as.numeric(
      auc(
        response = y_true,
        predictor = as.numeric(resUsStreaming$distances)
      )
    )
    
    print(paste0("AUC online ", round(auc_val, 4)))
  }
  
  setwd(criteres)
  critFile <- paste0(
    'CritStreamingUsonlinecorr-d', d,
    '-n', n,
    '-k', k,
    '-l', l,
    '-rho', rho1,
    '-r', r,
    "-cutoff", .95,
    "Ninit-", Ninit,
    "-batch",batch,
    "-cm", cm,
    '-sim', sim,
    ".RData"
  )
  if( r!= 0){
    save(
      erreurFrob,
      faux_negatifs,
      faux_positifs,
      ari,
      auc_val,
      file = critFile
    )
    
  } else {
    save(
      erreurFrob,
      faux_negatifs,
      faux_positifs,
      file = critFile
    )
    
  }
  
  fitFile <- paste0('FitStreamingUswithoutQuantonlinecorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,"Ninit-",1e2,"-batch",batch,'-sim', sim,".RData")
  
  setwd(resAlgo)
  
  print(fitFile)
  
  load(fitFile)
  
  erreurFrob <- norm(resUsStreaming$variance - Sigma0, "F")
  
  print(paste0("erreur norme de Frobenius streaming Us ", erreurFrob))
  
  t <- table(data$labelsVrais, resUsStreaming$outlier_labels)
  
  faux_positifs <- t[1, 2]
  
  if (r != 0) {
    faux_negatifs <- t[2, 1]
  } else {
    faux_negatifs <- 0
  }
  
  print(paste0("faux positifs streaming with quantile corr ", faux_positifs))
  print(paste0("faux négatifs streaming with quantile corr ", faux_negatifs))
  if(r != 0){
    # ARI
    ari <- adjustedRandIndex(
      data$labelsVrais,
      resUsStreaming$outlier_labels
    )
    
    print(paste0("ARI online us ", round(ari, 4)))
    
    # AUC
    y_true <- as.numeric(data$labelsVrais == unique(sort(data$labelsVrais))[2])
    
    auc_val <- as.numeric(
      auc(
        response = y_true,
        predictor = as.numeric(resUsStreaming$distances)
      )
    )
    
    print(paste0("AUC online ", round(auc_val, 4)))
  }
  setwd(criteres)
  critFile <- paste0(
    'CritStreamingUswithoutQuantonlinecorr-d', d,
    '-n', n,
    '-k', k,
    '-l', l,
    '-rho', rho1,
    '-r', r,
    "-cutoff", .95,
    "Ninit-", Ninit,
    "-batch",batch,
    '-sim', sim,
    ".RData"
  )
  if( r!= 0){
    save(
      erreurFrob,
      faux_negatifs,
      faux_positifs,
      ari,
      auc_val,
      file = critFile
    )
    
  } else {
    save(
      erreurFrob,
      faux_negatifs,
      faux_positifs,
      file = critFile
    )
    
  }
  
  #########################Offline######
  
  setwd(resAlgo)
  
  fitFile <- paste0('FitOfflinewithQuantcorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,"-cutoff",.95,"-cm",cm,'-sim', sim,".RData")
  
  print(fitFile)
  
  load(fitFile)
  
  erreurFrob <- norm(resOffline$variance - Sigma0, "F")
  
  print(paste0("erreur norme de Frobenius offline Us ", erreurFrob))
  
  t <- table(data$labelsVrais, resOffline$outlier_labels)
  
  faux_positifs <- t[1, 2]
  
  if (r != 0) {
    faux_negatifs <- t[2, 1]
  } else {
    faux_negatifs <- 0
  }
  
  print(paste0("faux positifs offline with quantile corr ", faux_positifs))
  print(paste0("faux négatifs offline with quantile corr ", faux_negatifs))
  if(r != 0){
    # ARI
    ari <- adjustedRandIndex(
      data$labelsVrais,
      resOffline$outlier_labels
    )
    
    print(paste0("ARI offline us ", round(ari, 4)))
    
    # AUC
    y_true <- as.numeric(data$labelsVrais == unique(sort(data$labelsVrais))[2])
    
    auc_val <- as.numeric(
      auc(
        response = y_true,
        predictor = as.numeric(resOffline$distances)
      )
    )
    
    print(paste0("AUC offline ", round(auc_val, 4)))
  }
  setwd(criteres)
  critFile <- paste0(
    'CritOfflinewithQuantcorr-d', d,
    '-n', n,
    '-k', k,
    '-l', l,
    '-rho', rho1,
    '-r', r,
    "-cutoff", .95,
    "-cm", cm,
    '-sim', sim,
    ".RData"
  )
  if( r!= 0){
    save(
      erreurFrob,
      faux_negatifs,
      faux_positifs,
      ari,
      auc_val,
      file = critFile
    )
    
  } else {
    save(
      erreurFrob,
      faux_negatifs,
      faux_positifs,
      file = critFile
    )
    
  }
  
  setwd(resAlgo)
  
  fitFile <- paste0('FitOfflineUswithoutQuantcorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,"-cutoff",.95,"-cm",cm,'-sim', sim,".RData")
  
  print(fitFile)
  
  load(fitFile)
  
  erreurFrob <- norm(resOffline$variance - Sigma0, "F")
  
  print(paste0("erreur norme de Frobenius offline Us ", erreurFrob))
  
  t <- table(data$labelsVrais, resOffline$outlier_labels)
  
  faux_positifs <- t[1, 2]
  
  if (r != 0) {
    faux_negatifs <- t[2, 1]
  } else {
    faux_negatifs <- 0
  }
  
  print(paste0("faux positifs offline without quantile corr ", faux_positifs))
  print(paste0("faux négatifs offline without quantile corr ", faux_negatifs))
  if(r != 0){
    # ARI
    ari <- adjustedRandIndex(
      data$labelsVrais,
      resOffline$outlier_labels
    )
    
    print(paste0("ARI offline us ", round(ari, 4)))
    
    # AUC
    y_true <- as.numeric(data$labelsVrais == unique(sort(data$labelsVrais))[2])
    
    auc_val <- as.numeric(
      auc(
        response = y_true,
        predictor = as.numeric(resOffline$distances)
      )
    )
    
    print(paste0("AUC online ", round(auc_val, 4)))
  }
  
    setwd(criteres)  
  
    critFile <- paste0(
    'CritOfflineUswithoutQuantcorr-d', d,
    '-n', n,
    '-k', k,
    '-l', l,
    '-rho', rho1,
    '-r', r,
    "-cutoff", .95,
  '-sim', sim,
    ".RData"
  )
  if( r!= 0){
    save(
      erreurFrob,
      faux_negatifs,
      faux_positifs,
      ari,
      auc_val,
      file = critFile
    )
    
  } else {
    save(
      erreurFrob,
      faux_negatifs,
      faux_positifs,
      file = critFile
    )
    
  }
  
  ############################OGK###########################
  
  setwd(resAlgo)
  
  
  fitFile <- paste0('FitOGK-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,"-cutoff",.95,'-sim', sim,".RData")
  
  print(fitFile)
  
  load(fitFile)
  
  
  erreurFrob <- norm(resOGK$cov - Sigma0, "F")
  
  print(paste0("erreur norme de Frobenius OGK Us ", erreurFrob))
  
  t <- table(data$labelsVrais, outlogk)
  
  faux_positifs <- t[1, 2]
  
  if (r != 0) {
    faux_negatifs <- t[2, 1]
  } else {
    faux_negatifs <- 0
  }
  
  print(paste0("faux positifs OGK without quantile corr ", faux_positifs))
  print(paste0("faux négatifs OGK without quantile corr ", faux_negatifs))
  if(r != 0){
    # ARI
    ari <- adjustedRandIndex(
      data$labelsVrais,
      outlogk
    )
    
    print(paste0("ARI OGK ", round(ari, 4)))
    
    # AUC
    y_true <- as.numeric(data$labelsVrais == unique(sort(data$labelsVrais))[2])
    
    auc_val <- as.numeric(
      auc(
        response = y_true,
        predictor = as.numeric(resOGK$distances)
      )
    )
    
    print(paste0("AUC OGK ", round(auc_val, 4)))
  }
  
  setwd(criteres)
  critFile <- paste0(
    'CritOGK-d', d,
    '-n', n,
    '-k', k,
    '-l', l,
    '-rho', rho1,
    '-r', r,
    "-cutoff", .95,
    '-sim', sim,
    ".RData"
  )
  if( r!= 0){
    save(
      erreurFrob,
      faux_negatifs,
      faux_positifs,
      ari,
      auc_val,
      file = critFile
    )
    
  } else {
    save(
      erreurFrob,
      faux_negatifs,
      faux_positifs,
      file = critFile
    )
    
  }
  
  
  
  setwd(resAlgo)
  
  
  fitFile <- paste0('FitMCD-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,"-cutoff",.95,'-sim', sim,".RData")
  
  print(fitFile)
  
  load(fitFile)
  
  
  erreurFrob <- norm(resMcd$cov - Sigma0, "F")
  
  print(paste0("erreur norme de Frobenius MCD ", erreurFrob))
  
  t <- table(data$labelsVrais, outlmcd)
  
  faux_positifs <- t[1, 2]
  
  if (r != 0) {
    faux_negatifs <- t[2, 1]
  } else {
    faux_negatifs <- 0
  }
  
  print(paste0("faux positifs MCD without quantile corr ", faux_positifs))
  print(paste0("faux négatifs MCD without quantile corr ", faux_negatifs))
  if(r != 0){
    # ARI
    ari <- adjustedRandIndex(
      data$labelsVrais,
      resOffline$outlier_labels
    )
    
    print(paste0("ARI MCD ", round(ari, 4)))
    
    # AUC
    y_true <- as.numeric(data$labelsVrais == unique(sort(data$labelsVrais))[2])
    
    auc_val <- as.numeric(
      auc(
        response = y_true,
        predictor = as.numeric(distmcd)
      )
    )
    
    print(paste0("AUC MCD ", round(auc_val, 4)))
  }
  
  setwd(criteres)
  critFile <- paste0(
    'CritMCD-d', d,
    '-n', n,
    '-k', k,
    '-l', l,
    '-rho', rho1,
    '-r', r,
    "-cutoff", .95,
    '-sim', sim,
    ".RData"
  )
  if( r!= 0){
    save(
      erreurFrob,
      faux_negatifs,
      faux_positifs,
      ari,
      auc_val,
      file = critFile
    )
    
  } else {
    save(
      erreurFrob,
      faux_negatifs,
      faux_positifs,
      file = critFile
    )
    
  }
  
  #############################Oracle#########
  
  setwd(resAlgo)
  
  
  fitFile <- paste0('FitOracle-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,"-cutoff",.95,'-sim', sim,".RData")
  
  print(fitFile)
  
  load(fitFile)
  
  t <- table(data$labelsVrais, outloracle)
  
  faux_positifs <- t[1, 2]
  
  if (r != 0) {
    faux_negatifs <- t[2, 1]
  } else {
    faux_negatifs <- 0
  }
  
  print(paste0("faux positifs Oracle ", faux_positifs))
  print(paste0("faux négatifs Oracle ", faux_negatifs))
  if(r != 0){
    # ARI
    ari <- adjustedRandIndex(
      data$labelsVrais,
      resOffline$outlier_labels
    )
    
    print(paste0("ARI Oracle ", round(ari, 4)))
    
    # AUC
    y_true <- as.numeric(data$labelsVrais == unique(sort(data$labelsVrais))[2])
    
    auc_val <- as.numeric(
      auc(
        response = y_true,
        predictor = as.numeric(distoracle)
      )
    )
    
    print(paste0("AUC Oracle ", round(auc_val, 4)))
  }
  
  setwd(criteres)
  critFile <- paste0(
    'CritOracle-d', d,
    '-n', n,
    '-k', k,
    '-l', l,
    '-rho', rho1,
    '-r', r,
    "-cutoff", .95,
    '-sim', sim,
    ".RData"
  )
  if( r!= 0){
    save(
      faux_negatifs,
      faux_positifs,
      ari,
      auc_val,
      file = critFile
    )
    
  } else {
    save(
        faux_negatifs,
      faux_positifs,
      file = critFile
    )
    
  }
  
  
  
  }
}
