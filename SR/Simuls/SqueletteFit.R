simDir <- 'DataSim/'
resDir <- 'FitSim/'




d <- 10
rList <- 5*(0:10)
#load(paste0('SimParmsGrid-d', d, '.Rdata'))
rList <- 5*(0:10)
simNb <- 1
n <- 1e4

# Fit (k for mu1) Other parameters are fixed
for(r in rList){for(k in kList){
  for(sim in 1:simNb){
    dataFile <- paste0('SimData-d', d, '-n', n, '-k', k, '-l', 1, '-rho', rho0, '-sim', sim)
    setwd("C:/Users/Paul GUILLOT/Documents/Simus/DataSim")
    load(file = dataFile)
    if(!file.exists(fitFile)){
      fitNaif <- SampleCovOnline(data$Z)
      fitUsOnline <- StreamingOutlierDetection(data$Z,batch = 1)
      fitUSStreaming = StreamingOutlierDetection(data$Z,batch = ncol(Z))
      setwd("C:/Users/Paul GUILLOT/Documents/Simus/Fitsim")
      save(fitNaif, fitUsOnline, fitUSStreaming,file=fitFile)
    }else{load(fitFile)}
  }
}}


# Fit (l for Sigma1) Other parameters are fixed
for(r in rList){for(l in lList){
  for(sim in 1:simNb){
    dataFile <- paste0('SimData-d', d, '-n', n, '-k', k, '-l', 1, '-rho', rho0, '-sim', sim)
    fitFile <- paste0('FitParms-d', d, '-n', n,  '-k', k, '-l', 1, '-rho', rho0, '-sim', sim)
    setwd("C:/Users/Paul GUILLOT/Documents/Simus/DataSim")
    load(file = dataFile)
    if(!file.exists(fitFile)){
      fitNaif <- SampleCovOnline(dataFile$Z)
      fitUsOnline <- StreamingOutlierDetection(dataFile$Z,batch = 1)
      fitUSStreaming = StreamingOutlierDetection(dataFile$Z,batch = ncol(Z))
      setwd("C:/Users/Paul GUILLOT/Documents/Simus/Fitsim")
      save(fitNaif, fitUsOnline, fitUSStreaming,file=fitFile)
    }else{load(fitFile)}
  }
}}


# Fit (rho for Sigma1) Other parameters are fixed
for(r in rList){for(rho in rho1List){
  for(sim in 1:simNb){
    dataFile <- paste0('SimData-d', d, '-n', n, '-k', k, '-l', 1, '-rho', rho, '-sim', sim,".RData")
    fitFile <- paste0('FitParms-d', d, '-n', n,  '-k', k, '-l', 1, '-rho', rho, '-sim', sim)
    setwd("C:/Users/Paul GUILLOT/Documents/Simus/DataSim")
    load(dataFile)
    if(!file.exists(fitFile)){
      fitNaif <- SampleCovOnline(data$Z)
      fitUsOnline <- StreamingOutlierDetection(data$Z,batch = 1)
      fitUSStreaming = StreamingOutlierDetection(data$Z,batch = ncol(Z))
      setwd("C:/Users/Paul GUILLOT/Documents/Simus/Fitsim")
      save(fitNaif, fitUsOnline, fitUSStreaming,file=fitFile)
    }else{load(fitFile)}
  }
}}



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

