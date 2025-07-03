##################################################
############Fit parameters file
##################################################



####################################################################
###########Repositories à adapter en fonction de votre configuration
####################################################################

simDir <- "D:/Simus/DataSim"
resDir <- 'D:/Simus/FitSim/'

###################################################################
#####load sim parameters
###################################################################

load("SimParmsGrid-n10000-d10.Rdata")

##########################################
###Extraction of the parameters
##########################################

kList = k1val

lList = l1val

rho1List = rho1val
###################################################################
#################Packages nécessaires##############################
##################################################################
library(Rcpp)
library(RMM)
library(Gmedian)
############################
#############FIchiers nécessaires
#################################
source("~/algosto/algorithmes.R")

sourceCpp("~/algosto/src/CodesRCpp.cpp")

#simNb = 5

l = 1
rho1 = 0.3

# Fit (k for mu1) Other parameters are fixed
for(r in rList){
for(k in kList){
  for(sim in 1:simNb){
      dataFile <- paste0('SimData-d', d, '-n', n, '-k', k, '-l', 1, '-rho', rho0,'-r',r , '-sim', sim,".RData")
    setwd(simDir)
    load(file = dataFile)
    fitFile <- paste0('FitParms-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
    if(!file.exists(fitFile)){
      temps_naif = system.time(
        {fitNaif <- SampleCovOnline(data$Z)})
      print(paste0("erreur CovNaif ", norm(fitNaif$Sigma - Sigma0,"F")))
      
      temps_online  = system.time({fitUsOnline <- StreamingOutlierDetection(data$Z,batch = 1)})
      temps_streaming = system.time({fitUSStreaming =StreamingOutlierDetection(data$Z,batch = ncol(data$Z))})
      print(paste0("erreur Streaming Us ", norm(fitUSStreaming$Sigma[n,,] - Sigma0,"F")))
      
      setwd(resDir)
      save(fitNaif, fitUsOnline, fitUSStreaming,temps_naif,temps_online,temps_streaming,file=fitFile)
    }else{load(fitFile)}
  }
}}
# Fit (l for Sigma1) Other parameters are fixed
simNb = 3
k = 0
rho0 = 0.3
rho1 = 0.3
for(r in rList){
  print(paste0("r =",r))
  #r = 40
  for(l in lList[5:length(lList)]){
  print(paste0("l = ",l))
    for(sim in 1:simNb){
    dataFile <- paste0('SimData-d', d, '-n', n, '-k', k, '-l', l, '-rho', rho1,'-r',r , '-sim', sim,".RData")
    setwd(simDir)
    load(file = dataFile)
    fitFile <- paste0('FitParms-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
    if(!file.exists(fitFile)){
      temps_naif = system.time(
        {fitNaif <- SampleCovOnline(data$Z)})
      print(paste0("erreur CovNaif ", norm(fitNaif$Sigma - Sigma0,"F")))
      temps_online  = system.time({fitUsOnline <- StreamingOutlierDetection(data$Z,batch = 1)})
      print(paste0("erreur Online Us ", norm(fitUsOnline$Sigma[n,,] - Sigma0,"F")))
      
      temps_streaming = system.time({fitUSStreaming =StreamingOutlierDetection(data$Z,batch = ncol(data$Z))})
      print(paste0("erreur Streaming Us ", norm(fitUSStreaming$Sigma[n,,] - Sigma0,"F")))
      
      setwd(resDir)
      save(fitNaif, fitUsOnline, fitUSStreaming,temps_naif,temps_online,temps_streaming,file=fitFile)
    }else{load(fitFile)}
  }
}}

simNb = 3
k = 0
rho0 = 0.3
rho1 = 0.3
for(r in rList){
  print(paste0("r =",r))
  #r = 40
  for(l in lList[1:4]){
    print(paste0("l = ",l))
    for(sim in 1:simNb){
      dataFile <- paste0('SimData-d', d, '-n', n, '-k', k, '-l', l, '-rho', rho1,'-r',r , '-sim', sim,".RData")
      setwd(simDir)
      load(file = dataFile)
      fitFile <- paste0('FitParms-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
      if(!file.exists(fitFile)){
        temps_naif = system.time(
          {fitNaif <- SampleCovOnline(data$Z)})
        print(paste0("erreur CovNaif ", norm(fitNaif$Sigma - Sigma0,"F")))
        temps_online  = system.time({fitUsOnline <- StreamingOutlierDetection(data$Z,batch = 1)})
        print(paste0("erreur Online Us ", norm(fitUsOnline$Sigma[n,,] - Sigma0,"F")))
        
        temps_streaming = system.time({fitUSStreaming =StreamingOutlierDetection(data$Z,batch = ncol(data$Z))})
        print(paste0("erreur Streaming Us ", norm(fitUSStreaming$Sigma[n,,] - Sigma0,"F")))
        
        setwd(resDir)
        save(fitNaif, fitUsOnline, fitUSStreaming,temps_naif,temps_online,temps_streaming,file=fitFile)
      }else{load(fitFile)}
    }
  }}



k = 0
l = 1
simNb = 1
# Fit (rho for Sigma1) Other parameters are fixed
for(r in rList){for(rho in rho1List){
  for(sim in 1:simNb){
    
    dataFile <- paste0('SimData-d', d, '-n', n, '-k', k, '-l', l, '-rho', rho,'-r',r , '-sim', sim,".RData")
    setwd(simDir)
    load(file = dataFile)
    fitFile <- paste0('FitParms-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho, '-r',r,'-sim', sim,".RData")
    if(!file.exists(fitFile)){
      temps_naif = system.time(
        {fitNaif <- SampleCovOnline(data$Z)})
      print(paste0("erreur Sigma naive ", norm(fitNaif$Sigma - Sigma0,"F")))
      temps_online  = system.time({fitUsOnline <- StreamingOutlierDetection(data$Z,batch = 1)})
      print(paste0("erreur Sigma Online Us ", norm(fitUsOnline$Sigma[n,,] - Sigma0,"F")))
      
      temps_streaming = system.time({fitUSStreaming =StreamingOutlierDetection(data$Z,batch = ncol(data$Z))})
      print(paste0("erreur Sigma Streaming Us ", norm(fitUSStreaming$Sigma[n,,] - Sigma0,"F")))
      
      setwd(resDir)
      save(fitNaif, fitUsOnline, fitUSStreaming,temps_naif,temps_online,temps_streaming,file=fitFile)
    }else{load(fitFile)}
  }
}
}

#Fit k and l rho is fixed 

rho = 0.3


# Fit (rho for Sigma1) Other parameters are fixed
for(r in rList){for(k in KLval){ 
  for(l in l1val){
   print(paste0("k = ",k,"l = ",l,"rho ",rho))
    for(sim in 1:simNb){
        
    dataFile <- paste0('SimData-d', d, '-n', n, '-k', k, '-l', l, '-rho', rho,'-r',r , '-sim', sim,".RData")
    setwd(simDir)
    load(file = dataFile)
    fitFile <- paste0('FitParms-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho, '-r',r,'-sim', sim,".RData")
    if(!file.exists(fitFile)){
      temps_naif = system.time(
        {fitNaif <- SampleCovOnline(data$Z)})
      print(paste0("erreur Sigma naive ", norm(fitNaif$Sigma - Sigma0,"F")))
      temps_online  = system.time({fitUsOnline <- StreamingOutlierDetection(data$Z,batch = 1)})
      print(paste0("erreur Sigma Online Us ", norm(fitUsOnline$Sigma[n,,] - Sigma0,"F")))
      
      temps_streaming = system.time({fitUSStreaming =StreamingOutlierDetection(data$Z,batch = ncol(data$Z))})
      print(paste0("erreur Sigma Streaming Us ", norm(fitUSStreaming$Sigma[n,,] - Sigma0,"F")))
      
      setwd(resDir)
      save(fitNaif, fitUsOnline, fitUSStreaming,temps_naif,temps_online,temps_streaming,file=fitFile)
    }else{load(fitFile)}
  }
}}
}




k = 0 ;l = 1; rho = 0.3 
# All parameters are fixed
for(r in rList){{
  for(sim in 1:simNb){
    
    dataFile <- paste0('SimData-d', d, '-n', n, '-k', k, '-l', 1, '-rho', rho,'-r',r , '-sim', sim,".RData")
    setwd(simDir)
    load(file = dataFile)
    fitFile <- paste0('FitParms-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho, '-r',r,'-sim', sim,".RData")
    if(!file.exists(fitFile)){
      temps_naif = system.time(
        {fitNaif <- SampleCovOnline(data$Z)})
      temps_online  = system.time({fitUsOnline <- StreamingOutlierDetection(data$Z,batch = 1)})
      temps_streaming = system.time({fitUSStreaming =StreamingOutlierDetection(data$Z,batch = ncol(data$Z))})
      setwd(resDir)
      save(fitNaif, fitUsOnline, fitUSStreaming,temps_naif,temps_online,temps_streaming,file=fitFile)
    }else{load(fitFile)}
  }
}
}


#Fit k and l varying rho fixed 
rho1 = 0.3
# Fit (all parms)
for(r in rList){for(k in kList){for(l in lList){
  for(sim in 1:simNb){
    print(paste0("r = ",r,"k = ",k,"l = ",l,"rho = ",rho1))
    dataFile <- paste0('SimData-d', d, '-n', n, '-k', k, '-l', 1, '-rho', rho0,'-r',r , '-sim', sim,".RData")
    setwd(simDir)
    load(file = dataFile)
    fitFile <- paste0('FitParms-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
    if(!file.exists(fitFile)){
      temps_naif = system.time(
        {fitNaif <- SampleCovOnline(data$Z)})
      print(paste0("Erreur Naif = ",norm(fitNaif$Sigma - Sigma0,"F")))
      temps_online  = system.time({fitUsOnline <- StreamingOutlierDetection(data$Z,batch = 1)})
      print(norm(fitUsOnline$Sigma[n,,] - Sigma0,"F"))
      print(paste0("Erreur online = ",norm(fitUsOnline$Sigma[n,,] - Sigma0,"F")))
      
      temps_streaming = system.time({fitUSStreaming =StreamingOutlierDetection(data$Z,batch = ncol(data$Z))})
      print(norm(fitUSStreaming$Sigma[n,,] - Sigma0,"F"))
      print(paste0("Erreur streaming = ",norm(fitUSStreaming$Sigma[n,,] - Sigma0,"F")))
      
      setwd(resDir)
      save(fitNaif, fitUsOnline, fitUSStreaming,temps_naif,temps_online,temps_streaming,file=fitFile)
    }else{load(fitFile)}
  }
}
}}