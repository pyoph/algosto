simDir <- 'C:/Users/Paul/Documents/Simus/DataSim'
resDir <- 'C:/Users/Paul/Documents/Simus/FitSim/'

d <- 10
rList <- 5*(0:10)
load(paste0('SimParmsGrid-d', d, '.Rdata'))
rList <- 5*(0:10)
simNb <- 10

erreursSigma = array(0,dim = c(length(rList),length(kList),3))
outliers_labelsRec = array(0,dim = c(n,length(rList),length(kList),3))
labels_vraisRec = matrix(0,n,length(rList))


l = 1
rho1 = 0.3

sim = 1
setwd("C:/Users/Paul/Documents/Simus/FitSim")
for (i in seq_along(rList) ){
  setwd(simDir)
  
  dataFile <- paste0('SimData-d', d, '-n', n, '-k', k, '-l', 1, '-rho', rho0,'-r',r , '-sim', sim,".RData")
  load(dataFile)
  
  labels_vraisRec[,i] = data$labelsVrais
  r = rList[i]
  
  for (j in seq_along(kList)){
  k = kList[j]

  setwd(resDir)
  
  fitFile <- paste0('FitParms-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
  load(fitFile)
  
  erreursSigma[i,j,1] = norm(Sigma0 - fitNaif$Sigma,"F")
  
  erreursSigma[i,j,2] = norm(Sigma0 - fitUsOnline$Sigma[n,,],"F")
  
  erreursSigma[i,j,3] = norm(Sigma0 - fitUSStreaming$Sigma[n,,],"F")
  
  outliers_labelsRec[,i,j,1] = fitNaif$outliers_labels
  outliers_labelsRec[,i,j,2] = fitUsOnline$outlier_labels
  outliers_labelsRec[,i,j,3] = fitUSStreaming$outlier_labels
  
  
  
  }
  
}


erreursSigma = array(0,dim = c(length(rList),length(kList),3))
outliers_labelsRec = array(0,dim = c(n,length(rList),length(kList),3))
labels_vraisRec = matrix(0,n,length(rList))


k = 0
rho1 = 0.3



erreursSigma = array(0,dim = c(length(rList),length(kList),3))
outliers_labelsRec = array(0,dim = c(n,length(rList),length(kList),3))
labels_vraisRec = matrix(0,n,length(rList))


k = 0
rho1 = 0.3

sim = 1
setwd("C:/Users/Paul/Documents/Simus/FitSim")
for (i in seq_along(rList) ){
  setwd(simDir)
  
  dataFile <- paste0('SimData-d', d, '-n', n, '-k', k, '-l', 1, '-rho', rho0,'-r',r , '-sim', sim,".RData")
  load(dataFile)
  
  labels_vraisRec[,i] = data$labelsVrais
  r = rList[i]
  
  for (j in seq_along(lList)[4:length(lList)]){
    l = lList[j]
    
    setwd(resDir)
    
    fitFile <- paste0('FitParms-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
    load(fitFile)
    
    erreursSigma[i,(j-3),1] = norm(Sigma0 - fitNaif$Sigma,"F")
    
    erreursSigma[i,(j-3),2] = norm(Sigma0 - fitUsOnline$Sigma[n,,],"F")
    
    erreursSigma[i,(j-3),3] = norm(Sigma0 - fitUSStreaming$Sigma[n,,],"F")
    
    outliers_labelsRec[,i,(j-3),1] = fitNaif$outliers_labels
    outliers_labelsRec[,i,(j-3),2] = fitUsOnline$outlier_labels
    outliers_labelsRec[,i,(j-3),3] = fitUSStreaming$outlier_labels
    
    
    
  }
  
}



erreursSigma = array(0,dim = c(length(rList),length(kList),3))
outliers_labelsRec = array(0,dim = c(n,length(rList),length(kList),3))
labels_vraisRec = matrix(0,n,length(rList))


k = 0
rho1 = 0.3

sim = 1
setwd("C:/Users/Paul/Documents/Simus/FitSim")
for (i in seq_along(rList) ){
  setwd(simDir)
  
  dataFile <- paste0('SimData-d', d, '-n', n, '-k', k, '-l', 1, '-rho', rho0,'-r',r , '-sim', sim,".RData")
  load(dataFile)
  
  labels_vraisRec[,i] = data$labelsVrais
  r = rList[i]
  
  for (j in seq_along(lList)[1:4]){
    l = lList[j]
    
    setwd(resDir)
    
    fitFile <- paste0('FitParms-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
    load(fitFile)
    
    erreursSigma[i,j,1] = norm(Sigma0 - fitNaif$Sigma,"F")
    
    erreursSigma[i,j,2] = norm(Sigma0 - fitUsOnline$Sigma[n,,],"F")
    
    erreursSigma[i,j,3] = norm(Sigma0 - fitUSStreaming$Sigma[n,,],"F")
    
    outliers_labelsRec[,i,j,1] = fitNaif$outliers_labels
    outliers_labelsRec[,i,j,2] = fitUsOnline$outlier_labels
    outliers_labelsRec[,i,j,3] = fitUSStreaming$outlier_labels
    
    
    
  }
  
}

k = 0

l  = 1

erreursSigma = array(0,dim = c(length(rList),length(kList),3))
outliers_labelsRec = array(0,dim = c(n,length(rList),length(kList),3))
labels_vraisRec = matrix(0,n,length(rList))

sim = 1
setwd("C:/Users/Paul/Documents/Simus/FitSim")
for (i in seq_along(rList) ){
  setwd(simDir)
  
  dataFile <- paste0('SimData-d', d, '-n', n, '-k', k, '-l', 1, '-rho', rho0,'-r',r , '-sim', sim,".RData")
  load(dataFile)
  
  labels_vraisRec[,i] = data$labelsVrais
  r = rList[i]
  
  for (j in seq_along(rho1List)){
    l = lList[j]
    
    setwd(resDir)
    
    fitFile <- paste0('FitParms-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
    load(fitFile)
    
    erreursSigma[i,j,1] = norm(Sigma0 - fitNaif$Sigma,"F")
    
    erreursSigma[i,j,2] = norm(Sigma0 - fitUsOnline$Sigma[n,,],"F")
    
    erreursSigma[i,j,3] = norm(Sigma0 - fitUSStreaming$Sigma[n,,],"F")
    
    outliers_labelsRec[,i,j,1] = fitNaif$outliers_labels
    outliers_labelsRec[,i,j,2] = fitUsOnline$outlier_labels
    outliers_labelsRec[,i,j,3] = fitUSStreaming$outlier_labels
    
    
    
  }
  
}



##################################
#Plot
#################################


# Créer le plot vide avec la première courbe (KL = 0)
plot(rList, erreursSigma[,1,1], type = "l", col = "blue", lwd = 2,
     xlab = "r", ylab = "Sigma estimation error",
     ylim = range(erreursSigma[,,1], na.rm = TRUE),
     main = "Sigma estimation error for different KL (l varying streaming)")

# Ajouter les autres courbes
lines(rList, erreursSigma[,2,2], col = "red", lwd = 2)
lines(rList, erreursSigma[,3,2], col = "green", lwd = 2)
lines(rList, erreursSigma[,4,2], col = "purple", lwd = 2)

# Ajouter une légende
legend("topright", legend = c("KL = 0", "KL = 1", "KL = 10", "KL = 100"),
       col = c("blue", "red", "green", "purple"), lwd = 2, bty = "n")



##########################################






# Fit
for(r in rList){for(k in kList){for(rho1 in rho1List){
  load(parmsFile)
  for(sim in 1:simNb){
    
    fitFile <- paste0(resDir, 'FitParms-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-sim', sim)
    load(dataFile)
    load(fitFile)
    
    
  }
}}}


