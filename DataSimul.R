##############
#Library needed
################

library(MASS)


#######################################################
#Directories paths à adapter selon votre configuration
#########################################################
paramDir = "D:/Simus/ParamsSim"
dataDir =  "D:/Simus/DataSim"


########################################################
#Loading of the parameters
#######################################################
load(paste0('SimParmsGrid-n',n,'-d', d, '.Rdata'))




####################
#Simulate data with different contamination rates
####################

#Extraction des paramètres
kList = k1val
lList = l1val
rho1List = rho1val


#######################################################################
#Simulation of a dataset with different parameters####################
######################################################################

genererEchantillon <- function(n, d, mu1, mu2,Sigma1, Sigma2,r) {
  # InitialSisation
  n1 <- floor((1 -r/100) * n)  # Taille du groupe non contaminé
  n2 <- n - n1         # Taille du groupe contaminé
  
  labels_mu1 <- rep(0, n1)  # Labels pour les vecteurs avec moyenne mu1
  labels_mu2 <- rep(1, n2)  # Labels pour les vecteurs avec moyenne mu2
  
  if (r > 0) {
    # Générer les vecteurs selon le type de contamination
      vecteurs_mu1 <- mvrnorm(n1, mu1, Sigma1)
      vecteurs_mu2 <- mvrnorm(n2, mu2, Sigma2)
    
      # Combinaison des vecteurs
      Z <- rbind(vecteurs_mu1, vecteurs_mu2)
      labelsVrais <- c(labels_mu1, labels_mu2)  
    }
    
    
    
  else {
    # Pas de contamination
    Z <- mvrnorm(n, mu1, Sigma1)
    labelsVrais <- rep(0, n)
  } # Mélanger aléatoirement les données
    set.seed(123)  # Pour garantir la reproductibilité
    indices <- sample(nrow(Z))
    Z <- Z[indices, ]
    labelsVrais <- labelsVrais[indices]
  
  
  return(list(Z = Z,labelsVrais = labelsVrais))
}


simNb = 1


# Simulation of datas
for(r in rList){
  #print(paste0("r =",r))
  
  for(k in kList){for(l in lList){for(rho1 in rho1List){
  for(sim in 1:simNb){
    set.seed(sim)
    
    parmsFile <- paste0('Parms-d', d, '-n', n,  '-k', k, '-l', l, '-rho1', rho1, '-r', r,'-sim', sim,".RData")
    setwd(paramDir)
    load(parmsFile)
    setwd(dataDir)               
    m1 <- (-1)^(1:d)/sqrt(d)
    sigmaSq0 <- (1:d); sigmaSq0 <- sigmaSq0 / mean(sigmaSq0)
    SigmaContamin <- diag(sqrt(sigmaSq0)) %*% toeplitz(rho1^(0:(d-1))) %*% diag(sqrt(sigmaSq0))
    dataFile <- paste0('SimData-d', d, '-n', n,  '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
    #print(paste0("nb outliers ",sum(data$labelsVrais == 1)))
    if(!file.exists(dataFile)){data = genererEchantillon(n,d,mu1 = mu0, mu2 = k*m1,Sigma1 = Sigma0,Sigma2 = l*SigmaContamin,r)
    }
    else {load(dataFile)}  
    # data <- Simul(d=d, r=r, k=k, rho1=rho1)
    save(data, file=dataFile)
  }
}}}}


