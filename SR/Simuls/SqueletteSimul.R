paramDir = 'ParamsSim/'
simDir <- 'DataSim/'

d <- 10
#load(paste0('SimParmsGrid-d', d, '.Rdata'))
rList <- 5*(0:10)
simNb <- 10

#Extraction des paramètres
kList = k1val
lList = l1val
rho1List = rho1val

Simus = "Simus"
setwd(Simus)


# Simul
for(r in rList){for(k in kList){for(l in lList){for(rho1 in rho1List){
  for(sim in 1:simNb){
    set.seed(sim)
    setwd("C:/Users/Paul GUILLOT/Documents/Simus/ParamsSim")
    parmsFile <- paste0('Parms-d', d, '-n', n,  '-k', k, '-l', l, '-rho1', rho1, '-sim', sim)
    # Sigma0 <- ; mu0 <- ; .....
    save(mu0, Sigma0, mu1, Sigma1, file=parmsFile)
    setwd("C:/Users/Paul GUILLOT/Documents/Simus/DataSim")               # Se déplace      setwd("DataSim")
    m1 <- (-1)^(1:d)/sqrt(d)
    sigmaSq0 <- (1:d); sigmaSq0 <- sigmaSq0 / mean(sigmaSq0)
    SigmaContamin <- diag(sqrt(sigmaSq0)) %*% toeplitz(rho^(0:(d-1))) %*% diag(sqrt(sigmaSq0))
    data = genererEchantillon(n,d,mu1 = mu0, mu2 = k*m1,Sigma1 = Sigma0,Sigma2 = l*SigmaContamin,p1 =1 - r/100,p2 = r/100,contamin = "moyenne_variance")
    dataFile <- paste0('SimData-d', d, '-n', n,  '-k', k, '-l', l, '-rho', rho1, '-sim', sim,".RData")
    
    # data <- Simul(d=d, r=r, k=k, rho1=rho1)
     save(data, file=dataFile)
  }
}}}}
