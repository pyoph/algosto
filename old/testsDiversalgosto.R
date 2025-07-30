simDir = "~/Simus/DataSim"

resDir = "~/Simus/FitSim"

explDir = "~/Simus/exploitResults"

setwd(resDir)


sim = 1

erreursSigmaNear = array(0,dim = c(n,length(rList),3))

outliersLabelsNear = array(0,dim = c(n,length(rList),3))
labelsVraisNear = array(0,dim = c(n,length(rList)))
faux_positifsNear = array(0,dim = c(length(rList),3))

faux_negatifsNear = array(0,dim = c(length(rList),3))


#Near scenario d = 10
k = 0.86;l=0.56;rho1 = 0.6
#Far scenario d = 100
k = 0.66;l=0.82;rho1 = 0.415

for (m in seq_along(rList)){
r = rList[m]
dataFile <- paste0('SimData-d', d, '-n', n, '-k', k, '-l', l, '-rho', rho1,'-r',r , '-sim', sim,".RData")

print(dataFile)

setwd(simDir)
if(!file.exists(dataFile)){contParam = ParmsF1(m1, k, l, rho1)
data = genererEchantillon(n,n,mu1 = mu0,mu2 = contParam$mu1,Sigma1 = Sigma0,Sigma2 = contParam$Sigma1,r)

}
else{load(dataFile)}
labelsVraisNear[,m] = data$labelsVrais
fitFile <- paste0('FitParms-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
temps_naif = system.time(
  
{resNaif = SampleCovOnline(data$Z)}
)
fitNaif = resNaif
for(s in (1:n)){
erreursSigmaNear[s,m,1] = norm(resNaif$SigmaIter[s,,] - Sigma0,"F")
outliersLabelsNear[,m,1] = resNaif$outliers_labels[s]
}
print(paste0("Erreur naive near ",erreursSigmaNear[n,m,1]))

t = table(data$labelsVrais,resNaif$outliers_labels)
if (r != 0) {faux_positifsNear[m,1] =  t[1,2]
faux_negatifsNear[m,1] = t[2,1]}
if(r == 0){faux_positifsNear[m,1] =  t[1,2]}

print(paste0("faux positifs near naive ",faux_positifsNear[m,1]))

print(paste0("faux négatifs near naive ",faux_negatifsNear[m,1]))


temps_online = (
  {
  if(d == 10){resUsOnline= StreamingOutlierDetection(data$Z,batch = 1)}
    if(d == 100){resUsOnline= StreamingOutlierDetection(data$Z,batch = 1,cutoff = 1.27 * qchisq(0.95, df = d))}
    }
)
fitUsOnline = resUsOnline
for(s in (1:n)){
erreursSigmaNear[s,m,2] = norm(resUsOnline$Sigma[s,,] - Sigma0,"F")
outliersLabelsNear[s,m,2] = resUsOnline$outlier_labels[s]
}
print(paste0("Erreur us online near ",erreursSigmaNear[n,m,2]))

t = table(data$labelsVrais,resUsOnline$outlier_labels)
if (r != 0) {faux_positifsNear[m,2] =  t[1,2]
faux_negatifsNear[m,2] = t[2,1]}
if(r == 0){faux_positifsNear[m,2] =  t[1,2]}


print(paste0("faux positifs near us online ",faux_positifsNear[m,2]))

print(paste0("faux négatifs near us online ",faux_negatifsNear[m,2]))



temps_streaming = system.time({
if(d == 10 ){resUsStreaming= StreamingOutlierDetection(data$Z,batch = ncol(data$Z))}
  if(d == 100){
resUsStreaming= StreamingOutlierDetection(data$Z,batch = sqrt(ncol(data$Z)),cutoff = 1.38*qchisq(0.95,df = d))}
})
fitUSStreaming = resUsStreaming
for(s in (1:n)){
erreursSigmaNear[s,m,3] =norm(resUsStreaming$Sigma[s,,] - Sigma0,"F")
outliersLabelsNear[s,m,2] = resUsStreaming$outlier_labels[s]
}
print(paste0("Erreur us streaming near ",erreursSigmaNear[n,m,3]))

t = table(data$labelsVrais,resUsStreaming$outlier_labels)
if(r == 0){faux_positifsNear[m,3] =  t[1,2]}
if (r != 0) {faux_positifsNear[m,3] =  t[1,2]
faux_negatifsNear[m,3] = t[2,1]}


print(paste0("faux positifs near us streaming ",faux_positifsNear[m,3]))

print(paste0("faux négatifs near us streaming ",faux_negatifsNear[m,3]))


setwd(resDir)
save(fitNaif, fitUsOnline, fitUSStreaming,temps_naif,temps_online,temps_streaming,file=fitFile)

}
erreursSigmaFar = array(0,dim = c(n,length(rList),3))
outliersLabelsFar = array(0,dim = c(n,length(rList),3))
labelsVraisFar = array(0,dim = c(n,length(rList)))
faux_positifsFar= array(0,dim = c(length(rList),3))
faux_negatifsFar = array(0,dim = c(length(rList),3))
labelsVraisFar = array(0,dim = c(n,length(rList)))


#Far scenario d = 10
k = 8.59;l=32;rho1 = 0.975
#Far scenario d = 100
k = 6.57;l = 19.02;rho1 = 0.845

for (m in seq_along(rList)) {
contParam = ParmsF1(m1, k, l, rho1)
r = rList[m]
dataFile <- paste0('SimData-d', d, '-n', n, '-k', k, '-l', l, '-rho', rho1,'-r',r , '-sim', sim,".RData")

print(dataFile)


setwd(simDir)

if(!file.exists(dataFile)){contParam = ParmsF1(m1, k, l, rho1)
data = genererEchantillon(n,n,mu1 = mu0,mu2 = contParam$mu1,Sigma1 = Sigma0,Sigma2 = contParam$Sigma1,r)
save(dataFile,file = dataFile)
}
else{load(dataFile)}
labelsVraisFar[,m] = data$labelsVrais


temps_naif = system.time({
resNaif = SampleCovOnline(data$Z)})
fitNaif = resNaif
for (s in (1:n)){
erreursSigmaFar[s,m,1] = norm(resNaif$SigmaIter[s,,] - Sigma0,"F")
outliersLabelsFar[s,m,1] = resNaif$outliers_labels[s]
}

print(paste0("Erreur naive far ",erreursSigmaFar[n,m,1]))

t = table(data$labelsVrais,resNaif$outliers_labels)
if (r != 0) {faux_positifsFar[m,1] =  t[1,2]
faux_negatifsFar[m,1] = t[2,1]}
if(r == 0){faux_positifsFar[m,1] =  t[1,2]}


print(paste0("faux positifs far us naive ",faux_positifsFar[m,1]))

print(paste0("faux négatifs far us naive ",faux_negatifsFar[m,1]))



temps_online = system.time(
  
  {
    if(d == 10){
  resUsOnline= StreamingOutlierDetection(data$Z,batch = 1)}
    
    if(d == 100){
      resUsOnline= StreamingOutlierDetection(data$Z,batch = 1,cutoff = 1.27 * qchisq(0.95, df = d))}
})
fitUsOnline = resUsOnline
t = table(data$labelsVrais,resUsOnline$outlier_labels)
if(r == 0){faux_positifsFar[m,2] =  t[1,2]}

if (r != 0) {faux_positifsFar[m,2] =  t[1,2]
faux_negatifsFar[m,2] = t[2,1]}

for (s in (1:n)){
  erreursSigmaFar[s,m,2] = norm(resUsOnline$Sigma[s,,] - Sigma0,"F")
  outliersLabelsFar[s,m,2] = resUsOnline$outlier_labels[s] 
}

print(paste0("Erreur online us far ",erreursSigmaFar[n,m,2]))

print(paste0("faux positifs far us online ",faux_positifsFar[m,2]))

print(paste0("faux négatifs far us online ",faux_negatifsFar[m,2]))

temps_streaming = system.time(
  {
  if(d == 10){
  resUsStreaming= StreamingOutlierDetection(data$Z,batch = ncol(data$Z))}
    if(d == 100){resUsStreaming= StreamingOutlierDetection(data$Z,batch = sqrt(ncol(data$Z)),cutoff = 1.38 * qchisq(0.95, df = d))}
    }
  
  )
fitUSStreaming = resUsStreaming
for(s in (1:n)){
erreursSigmaFar[s,m,3] = norm(resUsStreaming$Sigma[s,,] - Sigma0,"F")

outliersLabelsFar[s,m,3] = resUsStreaming$outlier_labels[s]}
t =table(data$labelsVrais,resUsStreaming$outlier_labels)

if (r != 0) {faux_positifsFar[m,3] =  t[1,2]
if(r == 0){faux_positifsFar[m,3] =  t[1,2]}

faux_negatifsFar[m,3] = t[2,1]}


print(paste0("Erreur streaming us far ",erreursSigmaFar[n,m,3]))

print(paste0("faux positifs far us streaming ",faux_positifsFar[m,3]))

print(paste0("faux négatifs far us streaming ",faux_negatifsFar[m,3]))



setwd(resDir)
save(fitNaif, fitUsOnline, fitUSStreaming,temps_naif,temps_online,temps_streaming,file=fitFile)

}

save(erreursSigmaNear,erreursSigmaFar,outliersLabelsNear,outliersLabelsFar,file = "res1runFarNeard100.Rdata")

