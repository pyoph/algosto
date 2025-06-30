simDir <- 'DataSim/'
resDir <- 'FitSim/'

d <- 10
rList <- 5*(0:10)
load(paste0('SimParmsGrid-d', d, '.Rdata'))
rList <- 5*(0:10)
simNb <- 10

# Fit
for(r in rList){for(k in kList){for(rho1 in rho1List){
  load(parmsFile)
  for(sim in 1:simNb){
    fitFile <- paste0(resDir, 'FitParms-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-sim', sim)
    load(dataFile)
    load(fitFile)
  }
}}}

