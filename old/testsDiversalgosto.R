simDir = "~/Simus/DataSim"

resDir = "~/Simus/FitSim"

explDir = "~/Simus/exploitResults"

#setwd(resDir)


##########################################
###Extraction of the parameters
##########################################

kList = k1val

lList = l1val

rho1List = rho1val

rho1ListNeg = rho1valNeg



test_outliers = function(distances,dim = 10,cutoff)
{
  outlab = rep(0,length(distances))
  
  
  for (i in (1:length(distances)))
  {
    if(distances[i] > cutoff){outlab[i] = 1}
  }
  return(outlab)
}

sim = 1
simNb = 1

erreursSigmaNear = array(0,dim = c(n,length(rList),3,simNb))
erreursKLNear = array(0,dim = c(n,length(rList),3,simNb))
outliersLabelsNear = array(0,dim = c(n,length(rList),3,simNb))
labelsVraisNear = array(0,dim = c(n,length(rList)))
faux_positifsNear = array(0,dim = c(length(rList),3,simNb))

faux_negatifsNear = array(0,dim = c(length(rList),3,simNb))


#Near scenario d = 10
if(d == 10) {k = 0.86;l=0.56;rho1 = 0.6}
#Near scenario d = 100
if(d == 100) {k = 0.66;l=0.82;rho1 = 0.415}
k = 0.86;l=0.56;rho1 = 0.6
#Near scenario d = 100
#k = 0.66;l=0.82;rho1 = 0.415

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
labelsVraisNear[,m] = data$labelsVrais
fitFile <- paste0('FitParms-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
temps_naif = system.time(
  
{resNaif = SampleCovOnline(data$Z)}
)
fitNaif = resNaif
for(s in (1:n)){
erreursSigmaNear[s,m,1,sim] = norm(resNaif$SigmaIter[s,,] - Sigma0,"F")
}
outliersLabelsNear[,m,1,sim] = resNaif$outliers_labels[s]

print(paste0("Erreur naive near ",erreursSigmaNear[n,m,1,sim]))

t = table(data$labelsVrais,resNaif$outliers_labels)
if (r != 0) {faux_positifsNear[m,1,sim] =  t[1,2]
faux_negatifsNear[m,1,sim] = t[2,1]}
if(r == 0){faux_positifsNear[m,1,sim] =  t[1,2]}

print(paste0("faux positifs near naive ",faux_positifsNear[m,1,sim]))

print(paste0("faux négatifs near naive ",faux_negatifsNear[m,1,sim]))


temps_online = (
  {
  if(d == 10){resUsOnline= StreamingOutlierDetection(data$Z,batch = 1)}
    if(d == 100){resUsOnline= StreamingOutlierDetection(data$Z,batch = 1,cutoff = 1.27 * qchisq(0.95, df = d))}
    }
)
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
save(Sigma_naive, Sigma_online, Sigma_str,temps_naif,temps_online,temps_streaming,file=fitFile)
  }

}
erreursSigmaFar = array(0,dim = c(n,length(rList),3,simNb))
erreursInvSigmaFar = array(0,dim = c(n,length(rList),3,simNb))
outliersLabelsFar = array(0,dim = c(n,length(rList),3,simNb))
labelsVraisFar = array(0,dim = c(n,length(rList)))
faux_positifsFar= array(0,dim = c(length(rList),3,simNb))
faux_negatifsFar = array(0,dim = c(length(rList),3,simNb))


#Far scenario d = 10
if(d == 10){k = 8.59;l=32;rho1 = 0.975}
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
contParam = ParmsF1(m1, k, l, rho1)
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
}
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
save(fitNaif, fitUsOnline, fitUSStreaming,temps_naif,temps_online,temps_streaming,file=fitFile)
  }


save(erreursSigmaNear,erreursSigmaFar,outliersLabelsNear,outliersLabelsFar,file = "res1runFarNeard10.Rdata")

