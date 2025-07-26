erreursSigmaNear = array(0,dim = c(n,length(rList),3))

outliersLabelsNear = array(0,dim = c(n,length(rList),3))

faux_positifsNear = array(0,dim = c(length(rList),3))

faux_negatifsNear = array(0,dim = c(length(rList),3))


k = 0.86;l=0.56;rho1 = 0.6

for (m in seq_along(rList)){
r = rList[m]
contParam = ParmsF1(m1, k, l, rho1)

data = genererEchantillon(n,n,mu1 = mu0,mu2 = contParam$mu1,Sigma1 = Sigma0,Sigma2 = contParam$Sigma1,r)

resNaif = SampleCovOnline(data$Z)
for(s in (1:n)){
erreursSigmaNear[s,m,1] = norm(resNaif$SigmaIter[s,,] - Sigma0,"F")
outliersLabelsNear[,m,1] = resNaif$outliers_labels[s]
}

t = table(data$labelsVrais,resNaif$outliers_labels)
if (r != 0) {faux_positifsNear[m,1] =  t[1,2]
faux_negatifsNear[m,1] = t[2,1]}
if(r == 0){faux_positifsNear[m,1] =  t[1,2]}
resUsOnline= StreamingOutlierDetection(data$Z,batch = 1)
for(s in (1:n)){
erreursSigmaNear[s,m,2] = norm(resUsOnline$Sigma[s,,] - Sigma0,"F")
outliersLabelsNear[s,m,2] = resUsOnline$outlier_labels[s]
}
t = table(data$labelsVrais,resUsOnline$outlier_labels)
if (r != 0) {faux_positifsNear[m,2] =  t[1,2]
faux_negatifsNear[m,2] = t[2,1]}
if(r == 0){faux_positifsNear[m,2] =  t[1,2]}
resUsStreaming= StreamingOutlierDetection(data$Z,batch = ncol(data$Z))
for(s in (1:n)){
erreursSigmaNear[s,m,3] =norm(resUsStreaming$Sigma[s,,] - Sigma0,"F")
outliersLabelsNear[,m,2] = resUsStreaming$outlier_labels[s]
}
t = table(data$labelsVrais,resUsStreaming$outlier_labels)
if(r == 0){faux_positifsNear[m,3] =  t[1,2]}
if (r != 0) {faux_positifsNear[m,3] =  t[1,2]
faux_negatifsNear[m,3] = t[2,1]}

}
erreursSigmaFar = array(0,dim = c(n,length(rList),3))
outliersLabelsFar = array(0,dim = c(n,length(rList),3))
faux_positifsFar= array(0,dim = c(length(rList),3))
faux_negatifsFar = array(0,dim = c(length(rList),3))



k = 8.59;l=32;rho1 = 0.975
for (m in seq_along(rList)) {
contParam = ParmsF1(m1, k, l, rho1)
r = rList[m]

data = genererEchantillon(n,n,mu1 = mu0,mu2 = contParam$mu1,Sigma1 = Sigma0,Sigma2 = contParam$Sigma1,r )

resNaif = SampleCovOnline(data$Z)
for (s in (1:n)){
erreursSigmaFar[s,m,1] = norm(resNaif$SigmaIter[s,,] - Sigma0,"F")
outliersLabelsFar[,m,1] = resNaif$outliers_labels[s]
}
t = table(data$labelsVrais,resNaif$outliers_labels)
if (r != 0) {faux_positifsFar[m,1] =  t[1,2]
faux_negatifsFar[m,1] = t[2,1]}
if(r == 0){faux_positifsFar[m,1] =  t[1,2]}
resUsOnline= StreamingOutlierDetection(data$Z,batch = 1)


t = table(data$labelsVrais,resUsOnline$outlier_labels)
if(r == 0){faux_positifsFar[m,2] =  t[1,2]}

if (r != 0) {faux_positifsFar[m,2] =  t[1,2]
faux_negatifsFar[m,2] = t[2,1]}

for (s in (1:n)){
  erreursSigmaFar[s,m,2] = norm(resUsOnline$Sigma[s,,] - Sigma0,"F")
  outliersLabelsFar[s,m,2] = resUsOnline$outlier_labels[s] 
}

resUsStreaming= StreamingOutlierDetection(data$Z,batch = ncol(data$Z))

for(s in (1:n)){
erreursSigmaFar[s,m,3] = norm(resUsStreaming$Sigma[s,,] - Sigma0,"F")

outliersLabelsFar[s,m,3] = resUsStreaming$outlier_labels[s]}
t =table(data$labelsVrais,resUsStreaming$outlier_labels)

if (r != 0) {faux_positifsFar[m,3] =  t[1,2]
faux_negatifsFar[m,3] = t[2,1]}

}
if(r == 0){faux_positifsFar[m,3] =  t[1,2]}

save(erreursSigmaNear,erreursSigmaFar,outliersLabelsNear,outliersLabelsFar,file = "res1runFarNear.Rdata")