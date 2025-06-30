simDir <- 'DataSim/'

d <- 10
load(paste0('SimParmsGrid-d', d, '.Rdata'))
rList <- 5*(0:10)
simNb <- 10

# Simul
for(r in rList){for(k in kList){for(l in lList){for(rho1 in rho1List){
  parmsFile <- paste0('Parms-d', d, '-n', n,  '-k', k, '-l', l, '-rho1', rho1, '-sim', sim)
  # Sigma0 <- ; mu0 <- ; .....
  save(mu0, Sigma0, mu1, Sigma1, file=parmsFile)
  for(sim in 1:simNb){
    set.seed(sim)
    dataFile <- paste0('SimData-d', d, '-n', n,  '-k', k, '-l', l, '-rho', rho1, '-sim', sim)
    # data <- Simul(d=d, r=r, k=k, rho1=rho1)
    # save(data, file=simFile)
  }
}}}}
