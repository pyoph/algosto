simDir <- 'DataSim/'
resDir <- 'FitSim/'




d <- 10
rList <- 5*(0:10)
#load(paste0('SimParmsGrid-d', d, '.Rdata'))
simNb <- 2
n <- 1e4

# Fit (k for mu1) Other parameters are fixed
for(r in rList){for(k in kList){
  for(sim in 1:simNb){
    dataFile <- paste0('SimData-d', d, '-n', n, '-k', k, '-l', 1, '-rho', rho0, '-sim', sim,".RData")
    setwd("C:/Users/Paul/Documents/Simus/DataSim")
    load(file = dataFile)
    fitFile <- paste0('FitParms-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-sim', sim,".RData")
    if(!file.exists(fitFile)){
      temps_naif = system.time(
      {fitNaif <- SampleCovOnline(data$Z)})
      temps_online  = system.time({fitUsOnline <- StreamingOutlierDetection(data$Z,batch = 1)})
      temps_streaming = system.time({fitUSStreaming =StreamingOutlierDetection(data$Z,batch = ncol(data$Z))})
      setwd("C:/Users/Paul/Documents/Simus/Fitsim")
      save(fitNaif, fitUsOnline, fitUSStreaming,temps_naif,temps_online,temps_streaming,file=fitFile)
    }else{load(fitFile)}
  }
}}


# Fit (l for Sigma1) Other parameters are fixed
for(r in rList){for(l in lList){
  for(sim in 1:simNb){
    dataFile <- paste0('SimData-d', d, '-n', n, '-k', k, '-l', 1, '-rho', rho0, '-sim', sim,".RData")
    setwd("C:/Users/Paul/Documents/Simus/DataSim")
    load(file = dataFile)
    fitFile <- paste0('FitParms-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-sim', sim,".RData")
    if(!file.exists(fitFile)){
      temps_naif = system.time(
        {fitNaif <- SampleCovOnline(data$Z)})
      temps_online  = system.time({fitUsOnline <- StreamingOutlierDetection(data$Z,batch = 1)})
      temps_streaming = system.time({fitUSStreaming =StreamingOutlierDetection(data$Z,batch = ncol(data$Z))})
      setwd("C:/Users/Paul/Documents/Simus/Fitsim")
      save(fitNaif, fitUsOnline, fitUSStreaming,temps_naif,temps_online,temps_streaming,file=fitFile)
    }else{load(fitFile)}
  }
}}


# Fit (rho for Sigma1) Other parameters are fixed
for(r in rList){for(rho in rho1List){
  for(sim in 1:simNb){
    dataFile <- paste0('SimData-d', d, '-n', n, '-k', k, '-l', 1, '-rho', rho0, '-sim', sim,".RData")
    setwd("C:/Users/Paul/Documents/Simus/DataSim")
    load(file = dataFile)
    fitFile <- paste0('FitParms-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-sim', sim,".RData")
    if(!file.exists(fitFile)){
      temps_naif = system.time(
        {fitNaif <- SampleCovOnline(data$Z)})
      temps_online  = system.time({fitUsOnline <- StreamingOutlierDetection(data$Z,batch = 1)})
      temps_streaming = system.time({fitUSStreaming =StreamingOutlierDetection(data$Z,batch = ncol(data$Z))})
      setwd("C:/Users/Paul/Documents/Simus/Fitsim")
      save(fitNaif, fitUsOnline, fitUSStreaming,temps_naif,temps_online,temps_streaming,file=fitFile)
    }else{load(fitFile)}
  }
  }
}



# Fit (all parms)
for(r in rList){for(k in kList){for(l in lList){for(rho1 in rho1List){
  for(sim in 1:simNb){
    dataFile <- paste0(simDir, 'SimDate-d', d, '-n', n,  '-k', k, '-l', l, '-rho', rho1, '-sim', sim)
    fitFile <- paste0(resDir, 'FitParms-d', d, '-n', n,  '-k', k, '-l', l, '-rho', rho1, '-sim', sim)
    load(dataFile)
    fit <- Algo(data)
    save(fit, file=fitFile)
  }
}}}}

