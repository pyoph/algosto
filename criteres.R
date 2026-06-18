
###########################
# Helpers
###########################

load_fit <- function(base_dir, file) {
  full <- file.path(base_dir, file)
  if (!file.exists(full)) stop(paste("Missing file:", full))
  load(full)
  return(resultats)
}

save_crit <- function(base_dir, file, crit) {
  save(crit, file = file.path(base_dir, file))
}

compute_one <- function(resAlgo, criteres, data, k, l, rho1, r, method, suffix_fit, suffix_crit, extra = "") {
  
  fitFile <- paste0(
    "Fit", method, "-d", d,
    "-n", n,
    "-k", k,
    "-l", l,
    "-rho", rho1,
    "-r", r,
    suffix_fit,
    extra,
    ".RData"
  )
  
  print(fitFile)
  
  resultats <- load_fit(resAlgo, fitFile)
  
  crit <- compute_criteres(
    variance   = resultats$variance,
    outlab     = resultats$outliers_labels,
    distances  = resultats$distances,
    labels     = as.numeric(data$labelsVrais),
    SigmaTrue  = Sigma0,
    r          = r
  )
  
  critFile <- paste0(
    "Crit", method, "-d", d,
    "-n", n,
    "-k", k,
    "-l", l,
    "-rho", rho1,
    "-r", r,
    suffix_crit,
    extra,
    ".RData"
  )
  
  save_crit(criteres, critFile, crit)
}
#######################Computation of the evaluation criteria#######################################




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
  
  crit = compute_criteres(variance = resultats$variance, outlab = resultats$outliers_labels, distances = resultats$distances, labels = as.numeric(data$labelsVrais), SigmaTrue = Sigma0, r = r)
  
  
  
  setwd(criteres)
  
  
  
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
  save(  crit,    file = critFile)
  
  

  setwd(resAlgo)
  
  fitFile <- paste0('FitSampleNaivewithoutonlinequantilecorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,"-Ninit",Ninit,'-sim', sim,".RData")
  
  load(file = fitFile)
 
  
  crit = compute_criteres(variance = resultats$variance, outlab = resultats$outliers_labels, distances = resultats$distances, labels = as.numeric(data$labelsVrais), SigmaTrue = Sigma0, r = r)
  
  
  
  setwd(criteres)
  critFile <- paste0(
    'CritSampleNaivewithoutonlinequantilecorr-d', d,
    '-n', n,
    '-k', k,
    '-l', l,
    '-rho', rho1,
    '-r', r,
    "-cutoff", 0.95,
    "Ninit-", Ninit,
    ".RData"
  )
  
  save(crit,file = critFile)
  
  
  ############################Online us#########################################
  
  
  setwd(resAlgo)
  
  
 fitFile =  paste0('FitOnlineUsQuantonlinecorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,"-cutoff",.95,"Ninit-",Ninit,"-batch",1,"-cm",cm,'-sim', sim,".RData")
  
  print(fitFile)
  
  
  load(fitFile)
  
  crit = compute_criteres(variance = resultats$variance, outlab = resultats$outliers_labels, distances = resultats$distances, labels = as.numeric(data$labelsVrais), SigmaTrue = Sigma0, r = r)
  
  
  
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
  
  save(crit,    file = critFile)
  
  
  setwd(resAlgo)
  
  fitFile <- paste0('FitOnlineUswithoutQuantonlinecorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,"-cutoff",.95,"Ninit-",1e2,"-batch",1,'-sim', sim,".RData")
  
  load(fitFile)
  
  print(fitFile)
  
  
  crit = compute_criteres(variance = resultats$variance, outlab = resultats$outliers_labels, distances = resultats$distances, labels = as.numeric(data$labelsVrais), SigmaTrue = Sigma0, r = r)
  
  
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

  save(crit,    file = critFile)
  
  
  
  
  
  ############################Streaming us#########################################
  
  fitFile <- paste0('FitStreamingUsonlineQuantcorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,"-cutoff",.95,"Ninit-",Ninit,"-batch",batch,"-cm",cm,'-sim', sim,".RData")
  
  setwd(resAlgo)
  
  print(fitFile)
  
  load(fitFile)
  
  setwd(criteres)
  
  crit = compute_criteres(variance = resultats$variance, outlab = resultats$outliers_labels, distances = resultats$distances, labels = as.numeric(data$labelsVrais), SigmaTrue = Sigma0, r = r)

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
  
  save(crit,    file = critFile)
  
  
  fitFile <- paste0('FitStreamingUswithoutQuantonlinecorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,"Ninit-",1e2,"-batch",batch,'-sim', sim,".RData")
  
  setwd(resAlgo)
  
  print(fitFile)
  
  load(fitFile)
  crit = compute_criteres(variance = resultats$variance, outlab = resultats$outliers_labels, distances = resultats$distances, labels = as.numeric(data$labelsVrais), SigmaTrue = Sigma0, r = r)
  
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
  save(crit,file = critFile)
  #########################Offline######
  
  setwd(resAlgo)
  
  fitFile <- paste0('FitOfflinewithQuantcorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,"-cutoff",.95,"-cm",cm,'-sim', sim,".RData")
  
  print(fitFile)
  
  load(fitFile)
  
  setwd(criteres)
  
  crit = compute_criteres(variance = resultats$variance, outlab = resultats$outliers_labels, distances = resultats$distances, labels = as.numeric(data$labelsVrais), SigmaTrue = Sigma0, r = r)
  
  
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
  save(crit,file = critFile)
  
  
  setwd(resAlgo)
  
  fitFile <- paste0('FitOfflineUswithoutQuantcorr-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,"-cutoff",.95,"-cm",cm,'-sim', sim,".RData")
  
  print(fitFile)
  
  load(fitFile)
  
  
  crit = compute_criteres(variance = resultats$variance, outlab = resultats$outliers_labels, distances = resultats$distances, labels = as.numeric(data$labelsVrais), SigmaTrue = Sigma0, r = r)
  
  
  
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
  
    save(crit,file = critFile)
    
  ############################OGK###########################
  
  setwd(resAlgo)
  
  
  fitFile <- paste0('FitOGK-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,"-cutoff",.95,'-sim', sim,".RData")
  
  print(fitFile)
  
  load(fitFile)
  crit = compute_criteres(variance = resultats$variance, outlab = resultats$outliers_labels, distances = resultats$distances, labels = as.numeric(data$labelsVrais), SigmaTrue = Sigma0, r = r)
  
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
  save(crit,file = critFile) 
  
  
  setwd(resAlgo)
  
  
  fitFile <- paste0('FitMCD-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,"-cutoff",.95,'-sim', sim,".RData")
  
  print(fitFile)
  
  load(fitFile)
  crit = compute_criteres(variance = resultats$variance, outlab = resultats$outliers_labels, distances = resultats$distances, labels = as.numeric(data$labelsVrais), SigmaTrue = Sigma0, r = r)
  
  
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
  save(crit,file = critFile) 
  
  
  #############################Oracle#########
  
  setwd(resAlgo)
  
  
  fitFile <- paste0('FitOracle-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,"-cutoff",.95,'-sim', sim,".RData")
  
  print(fitFile)
  
  load(fitFile)
  crit = compute_criteres(variance = resultats$variance, outlab = resultats$outliers_labels, distances = resultats$distances, labels = as.numeric(data$labelsVrais), SigmaTrue = Sigma0, r = r)
  
  
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
  save(crit,file = critFile) 
  
}}
