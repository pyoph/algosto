simDir = "~/Simus/DataSim"

resDir = "~/Simus/FitSim"

explDir = "~/Simus/exploitResults"

#setwd(resDir)


##########################################
###Extraction of the parameters
##########################################

kList = k1val

lList = l1valup1

rho1List = rho1val

rho1ListNeg = rho1valNeg

sim = 1
simNb = 1

erreursSigmaId = array(0,dim = c(n,length(rList),3,simNb))
erreursKLId = array(0,dim = c(n,length(rList),3,simNb))
outliersLabelsId = array(0,dim = c(n,length(rList),3,simNb))
labelsVraisId = array(0,dim = c(n,length(rList)))
faux_positifsId = array(0,dim = c(length(rList),3,simNb))

faux_negatifsId = array(0,dim = c(length(rList),3,simNb))


#Id scenario d = 10
if(d == 10) {k = 0;l=1;rho1 = 0.3}
#Near scenario d = 100
if(d == 100) {k = 0.66;l=0.82;rho1 = 0.415}
#k = 0.86;l=0.56;rho1 = 0.6
#Near scenario d = 100
#k = 0.66;l=0.82;rho1 = 0.415

for (m in seq_along(rList)){
  for(sim in(1:simNb)){
    r = rList[m]
    dataFile <- paste0('SimData-d', d, '-n', n, '-k', k, '-l', l, '-rho', rho1,'-r',r , '-sim', sim,".RData")
    
    print(dataFile)
    
    # setwd(simDir)
    # if(!file.exists(dataFile)){contParam = ParmsF1(m1, k, l, rho1)
    # data = genererEchantillon(n,n,mu1 = mu0,mu2 = contParam$mu1,Sigma1 = Sigma0,Sigma2 = contParam$Sigma1,r)
    # #save(dataFile)
    # }
    # else{load(dataFile)}
    
    contParam = ParmsF1(m1, k, l, rho1)
    data = genererEchantillon(n,n,mu1 = mu0,mu2 = contParam$mu1,Sigma1 = Sigma0,Sigma2 = contParam$Sigma1,r)
    #save(dataFile)
    labelsVraisId[,m] = data$labelsVrais
    fitFile <- paste0('FitParms-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
    temps_naif = system.time(
      
      {resNaif = SampleCovOnline(data$Z)}
    )
    fitNaif = resNaif
    for(s in (1:n)){
      erreursSigmaId[s,m,1,sim] = norm(resNaif$SigmaIter[s,,] - Sigma0,"F")
    }
    outliersLabelsId[,m,1,sim] = resNaif$outliers_labels
    
    print(paste0("Erreur naive Id ",erreursSigmaId[n,m,1,sim]))
    
    t = table(data$labelsVrais,resNaif$outliers_labels)
    if (r != 0) {faux_positifsId[m,1,sim] =  t[1,2]
    faux_negatifsId[m,1,sim] = t[2,1]}
    if(r == 0){faux_positifsId[m,1,sim] =  t[1,2]}
    
    print(paste0("faux positifs near naive ",faux_positifsId[m,1,sim]))
    
    print(paste0("faux négatifs near naive ",faux_negatifsId[m,1,sim]))
    
    
    temps_online = (
      {
        if(d == 10){resUsOnline= StreamingOutlierDetection(data$Z,batch = 1)}
        if(d == 100){resUsOnline= StreamingOutlierDetection(data$Z,batch = 1,cutoff = 1.27 * qchisq(0.95, df = d))}
      }
    )
    fitUsOnline = resUsOnline
    for(s in (1:n)){
      erreursSigmaId[s,m,2,sim] = norm(resUsOnline$Sigma[s,,] - Sigma0,"F")
      #outliersLabelsNear[s,m,2,sim] = resUsOnline$outlier_labels[s]
    }
    
    if(d == 10){outliersLabelsId[,m,2,sim] = resUsOnline$outlier_labels}
    if(d == 100){outliersLabelsId[,m,2,sim] = test_outliers(distances = resUsOnline$distances,cutoff = 1.29*qchisq(.95,df = 100))}
    
    print(paste0("Erreur us online near ",erreursSigmaNear[n,m,2,sim]))
    
    #t = table(data$labelsVrais,resUsOnline$outlier_labels)
    t = table(data$labelsVrais,outliersLabelsId[,m,2,sim])
    if (r != 0) {faux_positifsId[m,2,sim] =  t[1,2]
    faux_negatifsId[m,2,sim] = t[2,1]}
    if(r == 0){faux_positifsId[m,2,sim] = t[1,2]}
    
    
    print(paste0("faux positifs near us online ",faux_positifsId[m,2,sim]))
    
    print(paste0("faux négatifs near us online ",faux_negatifsId[m,2,sim]))
    
    
    
    temps_streaming = system.time({
      if(d == 10 ){resUsStreaming= StreamingOutlierDetection(data$Z,batch = ncol(data$Z))}
      if(d == 100){
        resUsStreaming= StreamingOutlierDetection(data$Z,batch = sqrt(ncol(data$Z)),cutoff = 1.38*qchisq(0.95,df = d))}
    })
    fitUSStreaming = resUsStreaming
    for(s in (1:n)){
      erreursSigmaId[s,m,3,sim] =norm(resUsStreaming$Sigma[s,,] - Sigma0,"F")
      #outliersLabelsNear[s,m,3,sim] = resUsStreaming$outlier_labels[s]
    }
    if(d == 10){outliersLabelsId[,m,3,sim]= resUsStreaming$outlier_labels}
    if(d == 100){
      outliersLabelsId[,m,3,sim] = test_outliers(distances = resUsStreaming$distances,cutoff = 1.38*qchisq(.95,df = 100))
    }
    print(paste0("Erreur us streaming near ",erreursSigmaId[n,m,3,sim]))
    
    #t = table(data$labelsVrais,resUsStreaming$outlier_labels)
    t = table(data$labelsVrais,outliersLabelsId[,m,3,sim])
    if(r == 0){faux_positifsId[m,3,sim] =  t[1,2]}
    if (r != 0) {faux_positifsId[m,3,sim] =  t[1,2]
    faux_negatifsId[m,3,sim] = t[2,1]}
    
    
    print(paste0("faux positifs near us streaming ",faux_positifsId[m,3,sim]))
    
    print(paste0("faux négatifs near us streaming ",faux_negatifsId[m,3,sim]))
    Sigma_naive = fitNaif$SigmaIter
    Sigma_online = fitUsOnline$Sigma
    Sigma_str = fitUSStreaming$Sigma
    
    setwd(resDir)
    #save(Sigma_naive, Sigma_online, Sigma_str,temps_naif,temps_online,temps_streaming,file=fitFile)
  }
  
}




erreursSigmaNear = array(0,dim = c(n,length(rList),3,simNb))
erreursSigmaOracle = array(0,dim = c(n,length(rList),3,simNb))
erreursKLNear = array(0,dim = c(n,length(rList),3,simNb))
outliersLabelsNear = array(0,dim = c(n,length(rList),3,simNb))
outliersLabelsNearOracle = array(0,dim = c(n,length(rList),simNb))
labelsVraisNear = array(0,dim = c(n,length(rList)))
faux_positifsNear = array(0,dim = c(length(rList),3,simNb))
faux_positifsOracle = array(0,dim = c(length(rList),simNb))
faux_negatifsNear = array(0,dim = c(length(rList),3,simNb))
faux_negatifsOracle = array(0,dim = c(length(rList),simNb))
temps = array(0,dim = c(length(rList),3,simNb)) 


#Near scenario d = 10
if(d == 10) {k = k1val[2];l=l1valup1[2];rho1 = rho1val[2]}
#Near scenario d = 100
if(d == 100) {k = 0.66;l=0.82;rho1 = 0.415}
#Near scenario d = 100
#k = 0.66;l=0.82;rho1 = 0.415

for (m in seq_along(rList)){
  for(sim in(1:simNb)){
r = rList[m]
dataFile <- paste0('SimData-d', d, '-n', n, '-k', k, '-l', l, '-rho', rho1,'-r',r , '-sim', sim,".RData")

print(dataFile)

# setwd(simDir)
# if(!file.exists(dataFile)){contParam = ParmsF1(m1, k, l, rho1)
# data = genererEchantillon(n,n,mu1 = mu0,mu2 = contParam$mu1,Sigma1 = Sigma0,Sigma2 = contParam$Sigma1,r)
# #save(dataFile)
# }
# else{load(dataFile)}

contParam = ParmsF1(m1, k, l, rho1)
data = genererEchantillon(n,n,mu1 = mu0,mu2 = contParam$mu1,Sigma1 = Sigma0,Sigma2 = contParam$Sigma1,r)
#save(dataFile)
labelsVraisNear[,m] = data$labelsVrais
fitFile <- paste0('FitParms-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
temps_naif = system.time(
  
{resNaif = SampleCovOnline(data$Z)}
)
fitNaif = resNaif
for(s in (1:n)){
erreursSigmaNear[s,m,1,sim] = norm(resNaif$SigmaIter[s,,] - Sigma0,"F")
if(t(data$Z[s,])%*% solve(Sigma0)%*% data$Z[s,] > qchisq(.95,df = d)) {
  outliersLabelsNearOracle[s,m,sim] = 1
  }

}

temps[m,1,sim] = temps_naif[1]
outliersLabelsNear[,m,1,sim] = resNaif$outliers_labels

print(paste0("Erreur naive near ",erreursSigmaNear[n,m,1,sim]))

t = table(data$labelsVrais,outliersLabelsNearOracle[,m,sim])
if (r != 0) {faux_positifsOracle[m,sim] =  t[1,2]
faux_negatifsOracle[m,sim] = t[2,1]}
if(r == 0){faux_positifsOracle[m,sim] =  t[1,2]}

print(paste0("faux positifs near naive ",faux_positifsNear[m,1,sim]))

print(paste0("faux négatifs near naive ",faux_negatifsNear[m,1,sim]))


print(paste0("faux positifs Oracle ",faux_positifsOracle[m,sim]))

print(paste0("faux négatifs Oracle ",faux_negatifsOracle[m,sim]))

t = table(data$labelsVrais,resNaif$outliers_labels)
if (r != 0) {faux_positifsNear[m,1,sim] =  t[1,2]
faux_negatifsNear[m,1,sim] = t[2,1]}
if(r == 0){faux_positifsNear[m,1,sim] =  t[1,2]}



temps_online = system.time(
  {
  if(d == 10){resUsOnline= StreamingOutlierDetection(data$Z,batch = 1)}
    if(d == 100){resUsOnline= StreamingOutlierDetection(data$Z,batch = 1,cutoff = 1.27 * qchisq(0.95, df = d))}
    }
)

temps[m,2,sim] = temps_online[1]
fitUsOnline = resUsOnline
for(s in (1:n)){
erreursSigmaNear[s,m,2,sim] = norm(resUsOnline$Sigma[s,,] - Sigma0,"F")
#outliersLabelsNear[s,m,2,sim] = resUsOnline$outlier_labels[s]
}

if(d == 10){outliersLabelsNear[,m,2,sim] = resUsOnline$outlier_labels}
if(d == 100){outliersLabelsNear[,m,2,sim] = test_outliers(distances = resUsOnline$distances,cutoff = 1.29*qchisq(.95,df = 100))}

print(paste0("Erreur us online near ",erreursSigmaNear[n,m,2,sim]))

#t = table(data$labelsVrais,resUsOnline$outlier_labels)
t = table(data$labelsVrais,outliersLabelsNear[,m,2,sim])
if (r != 0) {faux_positifsNear[m,2,sim] =  t[1,2]
faux_negatifsNear[m,2,sim] = t[2,1]}
if(r == 0){faux_positifsNear[m,2,sim] = t[1,2]}


print(paste0("faux positifs near us online ",faux_positifsNear[m,2,sim]))

print(paste0("faux négatifs near us online ",faux_negatifsNear[m,2,sim]))



temps_streaming = system.time({
if(d == 10 ){resUsStreaming= StreamingOutlierDetection(data$Z,batch = ncol(data$Z))}
  if(d == 100){
resUsStreaming= StreamingOutlierDetection(data$Z,batch = sqrt(ncol(data$Z)),cutoff = 1.38*qchisq(0.95,df = d))}
})

temps[m,3,sim] = temps_streaming[1]

fitUSStreaming = resUsStreaming
for(s in (1:n)){
erreursSigmaNear[s,m,3,sim] =norm(resUsStreaming$Sigma[s,,] - Sigma0,"F")
#outliersLabelsNear[s,m,3,sim] = resUsStreaming$outlier_labels[s]
}
if(d == 10){outliersLabelsNear[,m,3,sim]= resUsStreaming$outlier_labels}
if(d == 100){
outliersLabelsNear[,m,3,sim] = test_outliers(distances = resUsStreaming$distances,cutoff = 1.38*qchisq(.95,df = 100))
}
print(paste0("Erreur us streaming near ",erreursSigmaNear[n,m,3,sim]))

#t = table(data$labelsVrais,resUsStreaming$outlier_labels)
t = table(data$labelsVrais,outliersLabelsNear[,m,3,sim])
if(r == 0){faux_positifsNear[m,3,sim] =  t[1,2]}
if (r != 0) {faux_positifsNear[m,3,sim] =  t[1,2]
faux_negatifsNear[m,3,sim] = t[2,1]}


print(paste0("faux positifs near us streaming ",faux_positifsNear[m,3,sim]))

print(paste0("faux négatifs near us streaming ",faux_negatifsNear[m,3,sim]))
Sigma_naive = fitNaif$SigmaIter
Sigma_online = fitUsOnline$Sigma
Sigma_str = fitUSStreaming$Sigma

setwd(resDir)
#save(Sigma_naive, Sigma_online, Sigma_str,temps_naif,temps_online,temps_streaming,file=fitFile)
  }

}
setwd("~/Simus/FitSim")

file = paste0("results-k-",k,"l-",l,"rho1",rho1,".Rdata")
save(erreursSigmaNear,faux_negatifsNear,faux_positifsNear,outliersLabelsNear,labelsVraisNear,temps,file = file)
#######################Erreur Med scenario###########################

k = k1val[3];l = l1val[6];rho1 = rho1val[3]

erreursSigmaMed = array(0,dim = c(n,length(rList),3,simNb))
erreursInvSigmaMed = array(0,dim = c(n,length(rList),3,simNb))
outliersLabelsMed = array(0,dim = c(n,length(rList),3,simNb))
labelsVraisMed = array(0,dim = c(n,length(rList)))
faux_positifsMed= array(0,dim = c(length(rList),3,simNb))
faux_negatifsMed = array(0,dim = c(length(rList),3,simNb))


for (m in seq_along(rList)){
  for(sim in(1:simNb)){
    r = rList[m]
    dataFile <- paste0('SimData-d', d, '-n', n, '-k', k, '-l', l, '-rho', rho1,'-r',r , '-sim', sim,".RData")
    
    print(dataFile)
    
    setwd(simDir)
    if(!file.exists(dataFile)){contParam = ParmsF1(m1, k, l, rho1)
    data = genererEchantillon(n,n,mu1 = mu0,mu2 = contParam$mu1,Sigma1 = Sigma0,Sigma2 = contParam$Sigma1,r)
    #save(dataFile)
    }
    else{load(dataFile)}
    labelsVraisMed[,m] = data$labelsVrais
    fitFile <- paste0('FitParms-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
    temps_naif = system.time(
      
      {resNaif = SampleCovOnline(data$Z)}
    )
    fitNaif = resNaif
    for(s in (1:n)){
      erreursSigmaMed[s,m,1,sim] = norm(resNaif$SigmaIter[s,,] - Sigma0,"F")
    }
    outliersLabelsMed[,m,1,sim] = resNaif$outliers_labels[s]
    
    print(paste0("Erreur naive med ",erreursSigmaMed[n,m,1,sim]))
    
    t = table(data$labelsVrais,resNaif$outliers_labels)
    if (r != 0) {faux_positifsMed[m,1,sim] =  t[1,2]
    faux_negatifsMed[m,1,sim] = t[2,1]}
    if(r == 0){faux_positifsMed[m,1,sim] =  t[1,2]}
    
    print(paste0("faux positifs med naive ",faux_positifsMed[m,1,sim]))
    
    print(paste0("faux négatifs med naive ",faux_negatifsMed[m,1,sim]))
    
    
    temps_online = (
      {
        if(d == 10){resUsOnline= StreamingOutlierDetection(data$Z,batch = 1)}
        if(d == 100){resUsOnline= StreamingOutlierDetection(data$Z,batch = 1,cutoff = 1.27 * qchisq(0.95, df = d))}
      }
    )
    fitUsOnline = resUsOnline
    for(s in (1:n)){
      erreursSigmaMed[s,m,2,sim] = norm(resUsOnline$Sigma[s,,] - Sigma0,"F")
      #outliersLabelsNear[s,m,2,sim] = resUsOnline$outlier_labels[s]
    }
    
    if(d == 10){outliersLabelsMed[,m,2,sim] = resUsOnline$outlier_labels}
    if(d == 100){outliersLabelsMed[,m,2,sim] = test_outliers(distances = resUsOnline$distances,cutoff = 1.29*qchisq(.95,df = 100))}
    
    print(paste0("Erreur us online med ",erreursSigmaMed[n,m,2,sim]))
    
    #t = table(data$labelsVrais,resUsOnline$outlier_labels)
    t = table(data$labelsVrais,outliersLabelsMed[,m,2,sim])
    if (r != 0) {faux_positifsMed[m,2,sim] =  t[1,2]
    faux_negatifsMed[m,2,sim] = t[2,1]}
    if(r == 0){faux_positifsMed[m,2,sim] = t[1,2]}
    
    
    print(paste0("faux positifs med us online ",faux_positifsMed[m,2,sim]))
    
    print(paste0("faux négatifs med us online ",faux_negatifsMed[m,2,sim]))
    
    
    
    temps_streaming = system.time({
      if(d == 10 ){resUsStreaming= StreamingOutlierDetection(data$Z,batch = ncol(data$Z))}
      if(d == 100){
        resUsStreaming= StreamingOutlierDetection(data$Z,batch = sqrt(ncol(data$Z)),cutoff = 1.38*qchisq(0.95,df = d))}
    })
    fitUSStreaming = resUsStreaming
    for(s in (1:n)){
      erreursSigmaMed[s,m,3,sim] =norm(resUsStreaming$Sigma[s,,] - Sigma0,"F")
      #outliersLabelsNear[s,m,3,sim] = resUsStreaming$outlier_labels[s]
    }
    if(d == 10){outliersLabelsMed[,m,3,sim]= resUsStreaming$outlier_labels}
    if(d == 100){
      outliersLabelsMed[,m,3,sim] = test_outliers(distances = resUsStreaming$distances,cutoff = 1.38*qchisq(.95,df = 100))
    }
    print(paste0("Erreur us streaming med ",erreursSigmaMed[n,m,3,sim]))
    
    #t = table(data$labelsVrais,resUsStreaming$outlier_labels)
    t = table(data$labelsVrais,outliersLabelsMed[,m,3,sim])
    if(r == 0){faux_positifsMed[m,3,sim] =  t[1,2]}
    if (r != 0) {faux_positifsMed[m,3,sim] =  t[1,2]
    faux_negatifsMed[m,3,sim] = t[2,1]}
    
    
    print(paste0("faux positifs med us streaming ",faux_positifsMed[m,3,sim]))
    
    print(paste0("faux négatifs med us streaming ",faux_negatifsMed[m,3,sim]))
    Sigma_naive = fitNaif$SigmaIter
    Sigma_online = fitUsOnline$Sigma
    Sigma_str = fitUSStreaming$Sigma
    
    setwd(resDir)
    #save(Sigma_naive, Sigma_online, Sigma_str,temps_naif,temps_online,temps_streaming,file=fitFile)
  }
  
}


setwd("~/work/Simus/FitSim/")
save(erreursSigmaMed,faux_negatifsMed,faux_positifsMed,outliersLabelsMed,labelsVraisMed,file = "results_med_scenario.RData")


erreursSigmaFar = array(0,dim = c(n,length(rList),3,simNb))
erreursInvSigmaFar = array(0,dim = c(n,length(rList),3,simNb))
outliersLabelsFar = array(0,dim = c(n,length(rList),3,simNb))
labelsVraisFar = array(0,dim = c(n,length(rList)))
faux_positifsFar= array(0,dim = c(length(rList),3,simNb))
faux_negatifsFar = array(0,dim = c(length(rList),3,simNb))


#Far scenario d = 10
if(d == 10){k = 8.59;l=l1val[length(l1val)];rho1 = 0.975}
#Far scenario d = 100

if(d == 100){k = 6.57;l = 19.02;rho1 = 0.845}
#k = 6.57;l = 19.02;rho1 = 0.845

for (m in seq_along(rList)) {
  for (sim in (1:simNb)){
contParam = ParmsF1(m1, k, l, rho1)
r = rList[m]
dataFile <- paste0('SimData-d', d, '-n', n, '-k', k, '-l', l, '-rho', rho1,'-r',r , '-sim', sim,".RData")

print(dataFile)


setwd(simDir)
# 
# if(!file.exists(dataFile)){contParam = ParmsF1(m1, k, l, rho1)
# data = genererEchantillon(n,n,mu1 = mu0,mu2 = contParam$mu1,Sigma1 = Sigma0,Sigma2 = contParam$Sigma1,r)
# save(dataFile,file = dataFile)
# }
# else{load(dataFile)}
#contParam = ParmsF1(m1, k, l, rho1)
data = genererEchantillon(n,n,mu1 = mu0,mu2 = contParam$mu1,Sigma1 = Sigma0,Sigma2 = contParam$Sigma1,r)

labelsVraisFar[,m] = data$labelsVrais
fitFile <- paste0('FitParms-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")


temps_naif = system.time({
resNaif = SampleCovOnline(data$Z)})
fitNaif = resNaif
# for (s in (1:n)){
# erreursSigmaFar[s,m,1,sim] = norm(resNaif$SigmaIter[s,,] - Sigma0,"F")
# 
# VectPSigma = eigen(resNaif$SigmaIter[s,,])$vectors
# valPSigma = eigen(resNaif$SigmaIter[s,,])$values
# 
# 
# for(t in(1:d)){
#   
#   SigmaInvSqrt = 1/sqrt(valPSigma[t])*VectPSigma[,t]%*%t(VectPSigma[,t])
#   
# }
# erreursInvSigmaFar[s,m,1,sim] = norm(SigmaInvSqrt - (Sigma0)^(-0.5),"F")
# 
# }
for (s in (1:n)){
erreursSigmaFar[s,m,1,sim] = norm(resNaif$SigmaIter[s,,] - Sigma0,"F")}

outliersLabelsFar[,m,1,sim] = resNaif$outliers_labels

print(paste0("Erreur naive far ",erreursSigmaFar[n,m,1,sim]))


# if(erreursSigmaInvFar[n,m,1,sim] < 1e12){
# print(paste0("Erreur naive sqrt inv far ",erreursInvSigmaFar[n,m,1,sim]))
# }
# 
# if(erreursSigmaInvFar[n,m,1,sim] >= 1e12){
#   print("NAN")
# }



t = table(data$labelsVrais,resNaif$outliers_labels)
if (r != 0) {faux_positifsFar[m,1,sim] =  t[1,2]
faux_negatifsFar[m,1,sim] = t[2,1]}
if(r == 0){faux_positifsFar[m,1,sim] =  t[1,2]}


print(paste0("faux positifs far us naive ",faux_positifsFar[m,1,sim]))

print(paste0("faux négatifs far us naive ",faux_negatifsFar[m,1,sim]))



temps_online = system.time(
  
  {
    if(d == 10){
  resUsOnline= StreamingOutlierDetection(data$Z,batch = 1)}
    
    if(d == 100){
      resUsOnline= StreamingOutlierDetection(data$Z,batch = 1,cutoff = 1.27 * qchisq(0.95, df = d))}
})
fitUsOnline = resUsOnline
# 
# for (s in (1:n)){
#   erreursSigmaFar[s,m,2,sim] = norm(resUsOnline$Sigma[s,,] - Sigma0,"F")
#   
#   VectPSigma = eigen(resUsOnline$Sigma[s,,])$vectors
#   valPSigma = eigen(resUsOnline$Sigma[s,,])$values
#   
#   
#   for(t in(1:d)){
#     
#     SigmaInvSqrt = 1/sqrt(valPSigma[t])*VectPSigma[,t]%*%t(VectPSigma[,t])
#     
#   }
#   #outliersLabelsFar[s,m,2,sim] = resUsOnline$outlier_labels[s] 
#   
#   erreursInvSigmaFar[s,m,2,sim] = norm(SigmaInvSqrt - (Sigma0)^(-1/2),"F")

# 
# if(erreursSigmaInvFar[n,m,2,sim] < 1e12){
#   print(paste0("Erreur online us sqrt inv far ",erreursInvSigmaFar[n,m,2,sim]))
# }

# if(erreursSigmaInvFar[n,m,2,sim] >= 1e12){
#   print("NAN")
# }
   for (s in (1:n)){
     erreursSigmaFar[s,m,2,sim] = norm(resUsOnline$Sigma[s,,] - Sigma0,"F")}

if(d == 10){outliersLabelsFar[,m,2,sim]= resUsOnline$outlier_labels}
if(d == 100){
  outliersLabelsFar[,m,2,sim] = test_outliers(distances = resUsOnline$distances,cutoff = 1.38*qchisq(.95,df = 100))
}
#outliersLabelsFar[,m,2,sim] = test_outliers(distances = resUsOnline$distances,cutoff = 1.27 * qchisq(0.95, df = d))
#t = table(data$labelsVrais,resUsOnline$outlier_labels)
t = table(data$labelsVrais,outliersLabelsFar[,m,2,sim])
if(r == 0){faux_positifsFar[m,2,sim] =  t[1,2]}

if (r != 0) {faux_positifsFar[m,2,sim] =  t[1,2]
faux_negatifsFar[m,2,sim] = t[2,1]}

print(paste0("Erreur online us far ",erreursSigmaFar[n,m,2,sim]))

print(paste0("faux positifs far us online ",faux_positifsFar[m,2,sim]))

print(paste0("faux négatifs far us online ",faux_negatifsFar[m,2,sim]))

temps_streaming = system.time(
  {
  if(d == 10){
  resUsStreaming = StreamingOutlierDetection(data$Z,batch = ncol(data$Z))}
    if(d == 100){resUsStreaming= StreamingOutlierDetection(data$Z,batch = sqrt(ncol(data$Z)),cutoff = 1.38 * qchisq(0.95, df = d))}
    }
  
  )
fitUSStreaming = resUsStreaming
for(s in (1:n)){
erreursSigmaFar[s,m,3,sim] = norm(resUsStreaming$Sigma[s,,] - Sigma0,"F")
# 
# VectPSigma = eigen(resUsStreaming$Sigma[s,,])$vectors
# valPSigma = eigen(resUsStreaming$Sigma[s,,])$values
# 
# 
# for(t in(1:d)){
#   
#   SigmaInvSqrt = 1/sqrt(valPSigma[t])*VectPSigma[,t]%*%t(VectPSigma[,t])
#   
# }
# erreursInvSigmaFar[s,m,3,sim] = norm(SigmaInvSqrt - (Sigma0)^(-1/2),"F")
#outliersLabelsFar[s,m,3,sim] = resUsStreaming$outlier_labels[s]
}

if(d == 10){outliersLabelsFar[,m,3,sim]= resUsStreaming$outlier_labels}
if(d == 100){
  outliersLabelsFar[,m,3,sim] = test_outliers(distances = resUsStreaming$distances,cutoff = 1.38*qchisq(.95,df = 100))
}
# #t =table(data$labelsVrais,resUsStreaming$outlier_labels)
# if(erreursInvSigmaFar[n,m,3,sim] < 1e12)
#   {print("Erreur Sigma inv streaming ", erreursInvSigmaFar[n,m,3,sim])}
# if (erreursInvSigmaFar[n,m,3,sim] > 1e12)
#   {print("NAN")}
t =table(data$labelsVrais,outliersLabelsFar[,m,3,sim])

if (r != 0) {
faux_positifsFar[m,3,sim] =  t[1,2]
faux_negatifsFar[m,3,sim] = t[2,1]
}
if(r == 0){faux_positifsFar[m,3,sim] =  t[1,2]}


print(paste0("Erreur streaming us far ",erreursSigmaFar[n,m,3,sim]))

print(paste0("faux positifs far us streaming ",faux_positifsFar[m,3,sim]))

print(paste0("faux négatifs far us streaming ",faux_negatifsFar[m,3,sim]))



setwd(resDir)
#save(fitNaif, fitUsOnline, fitUSStreaming,temps_naif,temps_online,temps_streaming,file=fitFile)
  }}


#save(erreursSigmaNear,erreursSigmaFar,outliersLabelsNear,outliersLabelsFar,file = "res1runFarNeard10.Rdata")

setwd("~/work/Simus/FitSim/")

save(erreursSigmaFar,faux_negatifsFar,faux_positifsFar,outliersLabelsFar,labelsVraisFar,file = "results_far_scenario.Rdata")

# 
# setwd("~/algosto/resultsSelectedScenarios")
# save(erreursSigmaFar,faux_negatifsFar,faux_positifsFar,file = "erreursFarScenarios.RData")

#########################Moyenne erreurs, faux négatifs, faux positifs###############################################################
moyenne_erreursSigmaFar = apply(erreursSigmaFar,c(1,2,3),mean)
moyenne_erreursSigmaMed = apply(erreursSigmaMed,c(1,2,3),mean)
moyenne_erreursSigmaNear = apply(erreursSigmaNear,c(1,2,3),mean)
moyenne_erreursSigmaId = apply(erreursSigmaId,c(1,2,3),mean)



moyenne_faux_positifsFar = apply(faux_positifsFar,c(1,2),mean) 
moyenne_faux_positifsMed = apply(faux_positifsMed,c(1,2),mean) 
moyenne_faux_positifsNear = apply(faux_positifsNear,c(1,2),mean) 
moyenne_faux_positifsId = apply(faux_positifsId,c(1,2),mean) 
    

moyenne_faux_negatifsFar = apply(faux_negatifsFar,c(1,2),mean) 
moyenne_faux_negatifsMed = apply(faux_negatifsMed,c(1,2),mean) 
moyenne_faux_negatifsNear = apply(faux_negatifsNear,c(1,2),mean) 
moyenne_faux_negatifsId = apply(faux_negatifsId,c(1,2),mean) 


majority_vote_Far <- ifelse(apply(outliersLabelsFar, c(1,2,3), function(x) sum(x) >= 50), 1, 0)


