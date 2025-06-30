simDir <- 'DataSim/'
resDir <- 'FitSim/'




d <- 10
rList <- 5*(0:10)
load(paste0('SimParmsGrid-d', d, '.Rdata'))
rList <- 5*(0:10)
simNb <- 10
n <- 1e4

# Fit (k for mu1)
for(r in rList){for(k in kList){
  for(sim in 1:simNb){
    dataFile <- paste0(simDir, 'SimDate-d', d, '-n', n, '-k', k, '-l', 1, '-rho', rho0, '-sim', sim)
    fitFile <- paste0(resDir, 'FitParms-d', d, '-n', n,  '-k', k, '-l', 1, '-rho', rho0, '-sim', sim)
    load(dataFile)
    if(!file.exists(fitFile)){
      fitNaif <- Naif(data)
      fitUs <- Algo(data)
      save(fitNaif, fitUs, file=fitFile)
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

