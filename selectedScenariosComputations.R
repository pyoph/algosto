simDir = "~/work/Simus/DataSim"

resDir = "~/work/Simus/FitSim"

explDir = "~/work/Simus/exploitResults"

#setwd(resDir)


##########################################
###Extraction of the parameters
##########################################

kList = k1val

lList = l1valup1

rho1List = rho1val

rho1ListNeg = rho1valNeg

sim = 1
simNb = 1e2
cutoff = qchisq(p = .95, df = d)
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
    
    
    
    print(paste0("faux positifs Oracle ",faux_positifsOracle[m,sim]))
    
    print(paste0("faux négatifs Oracle ",faux_negatifsOracle[m,sim]))
    
    t = table(data$labelsVrais,resNaif$outliers_labels)
    
    if (r != 0) {
      
      faux_positifsNear[m,1,sim] =  t[1,2]
      faux_negatifsNear[m,1,sim] = t[2,1]
    }
    
    if(r == 0){faux_positifsNear[m,1,sim] =  t[1,2]}
    
    print(paste0("faux positifs near naive ",faux_positifsNear[m,1,sim]))
    
    print(paste0("faux négatifs near naive ",faux_negatifsNear[m,1,sim]))
    
    
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
      if(resUsOnline$distances[s] > cutoff){outliersLabelsNear[s,m,2,sim] = 1}
      if(resUsStreaming$distances[s] > cutoff){outliersLabelsNear[s,m,3,sim] = 1}
    }
    
    # if(d == 10){outliersLabelsNear[,m,3,sim]= resUsStreaming$outlier_labels}
    # if(d == 100){
    # outliersLabelsNear[,m,3,sim] = test_outliers(distances = resUsStreaming$distances,cutoff = 1.38*qchisq(.95,df = 100))
    # }
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
setwd("~")

file = paste0("results-k-",k,"l-",l,"rho1",rho1,"oracle",".Rdata")
save(erreursSigmaNear,faux_negatifsNear,faux_positifsNear,outliersLabelsNear,labelsVraisNear,temps,file = file)
#######################Erreur Med scenario###########################

k = k1val[3];l = l1valup1[3];rho1 = rho1val[3]

erreursSigmaMed = array(0,dim = c(n,length(rList),3,simNb))
erreursInvSigmaMed = array(0,dim = c(n,length(rList),3,simNb))
outliersLabelsMed = array(0,dim = c(n,length(rList),3,simNb))
outliersLabelsOracleMed1 = array(0,dim = c(n,length(rList),simNb))
labelsVraisMed = array(0,dim = c(n,length(rList)))
faux_positifsMed= array(0,dim = c(length(rList),3,simNb))
faux_negatifsMed = array(0,dim = c(length(rList),3,simNb))
faux_positifsOracleMed = array(0,dim = c(length(rList),simNb))
faux_negatifsOracleMed = array(0,dim = c(length(rList),simNb))


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
      if(t(data$Z[s,])%*% solve(Sigma0)%*% data$Z[s,] > qchisq(.95,df = d)) {
        outliersLabelsOracleMed1[s,m,sim] = 1
      }
      
    }
    outliersLabelsMed[,m,1,sim] = resNaif$outliers_labels
    
    print(paste0("Erreur naive med ",erreursSigmaMed[n,m,1,sim]))
    
    t = table(data$labelsVrais,resNaif$outliers_labels)
    if (r != 0) {faux_positifsMed[m,1,sim] =  t[1,2]
    faux_negatifsMed[m,1,sim] = t[2,1]}
    if(r == 0){faux_positifsMed[m,1,sim] =  t[1,2]}
    
    print(paste0("faux positifs med naive ",faux_positifsMed[m,1,sim]))
    
    print(paste0("faux négatifs med naive ",faux_negatifsMed[m,1,sim]))
    
    t = table(data$labelsVrais,outliersLabelsOracleMed1[,m,sim])
    if (r != 0) {faux_positifsOracleMed[m,sim] =  t[1,2]
    faux_negatifsOracleMed[m,sim] = t[2,1]}
    if(r == 0){faux_positifsOracleMed[m,sim] =  t[1,2]}
    
    
    print(paste0("faux positifs Oracle ",faux_positifsOracleMed[m,sim]))
    
    print(paste0("faux négatifs Oracle ",faux_negatifsOracleMed[m,sim]))
    
    
    temps_online = (
      {
        if(d == 10){resUsOnline= StreamingOutlierDetection(data$Z,batch = 1)}
        if(d == 100){resUsOnline= StreamingOutlierDetection(data$Z,batch = 1,cutoff = 1.27 * qchisq(0.95, df = d))}
      }
    )
    fitUsOnline = resUsOnline
    
    for(s in (1:n)){
      erreursSigmaMed[s,m,3,sim] =norm(resUsStreaming$Sigma[s,,] - Sigma0,"F")
      if(resUsOnline$distances[s] > cutoff){outliersLabelsMed[s,m,2,sim] = 1}
      if(resUsStreaming$distances[s] > cutoff){outliersLabelsMed[s,m,3,sim] = 1}
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


setwd("~")

file = paste0("results-k-",k,"l-",l,"rho1",rho1,".Rdata")
save(erreursSigmaMed,faux_negatifsOracleMed,faux_positifsOracleMed,faux_negatifsMed,faux_positifsMed,outliersLabelsMed,labelsVraisMed,file = file)

######################################Scen Med 2######################################
#######################Erreur Med scenario###########################

k = k1val[4];l = l1valup1[4];rho1 = rho1val[4]

erreursSigmaMed2 = array(0,dim = c(n,length(rList),3,simNb))
erreursInvSigmaMed2 = array(0,dim = c(n,length(rList),3,simNb))
outliersLabelsMed2 = array(0,dim = c(n,length(rList),3,simNb))
outliersLabelsOracleMed2 = array(0,dim = c(n,length(rList),simNb))
labelsVraisMed2 = array(0,dim = c(n,length(rList)))
faux_positifsMed2= array(0,dim = c(length(rList),3,simNb))
faux_negatifsMed2 = array(0,dim = c(length(rList),3,simNb))
faux_positifsOracleMed2 = array(0,dim = c(length(rList),simNb))
faux_negatifsOracleMed2 = array(0,dim = c(length(rList),simNb))


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
    labelsVraisMed2[,m] = data$labelsVrais
    fitFile <- paste0('FitParms-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
    temps_naif = system.time(
      
      {resNaif = SampleCovOnline(data$Z)}
    )
    fitNaif = resNaif
    for(s in (1:n)){
      erreursSigmaMed2[s,m,1,sim] = norm(resNaif$SigmaIter[s,,] - Sigma0,"F")
      if(t(data$Z[s,])%*% solve(Sigma0)%*% data$Z[s,] > qchisq(.95,df = d)) {
        outliersLabelsOracleMed2[s,m,sim] = 1
      }
      
    }
    outliersLabelsMed2[,m,1,sim] = resNaif$outliers_labels
    
    print(paste0("Erreur naive med ",erreursSigmaMed2[n,m,1,sim]))
    
    t = table(data$labelsVrais,resNaif$outliers_labels)
    if (r != 0) {faux_positifsMed2[m,1,sim] =  t[1,2]
    faux_negatifsMed2[m,1,sim] = t[2,1]}
    if(r == 0){faux_positifsMed2[m,1,sim] =  t[1,2]}
    
    print(paste0("faux positifs med naive ",faux_positifsMed2[m,1,sim]))
    
    print(paste0("faux négatifs med naive ",faux_negatifsMed2[m,1,sim]))
    
    t = table(data$labelsVrais,outliersLabelsOracleMed2[,m,sim])
    if (r != 0) {faux_positifsOracleMed2[m,sim] =  t[1,2]
    faux_negatifsOracleMed2[m,sim] = t[2,1]}
    if(r == 0){faux_positifsOracleMed2[m,sim] =  t[1,2]}
    
    
    print(paste0("faux positifs Oracle ",faux_positifsOracleMed2[m,sim]))
    
    print(paste0("faux négatifs Oracle ",faux_negatifsOracleMed2[m,sim]))
    
    
    temps_online = (
      {
        if(d == 10){resUsOnline= StreamingOutlierDetection(data$Z,batch = 1)}
        if(d == 100){resUsOnline= StreamingOutlierDetection(data$Z,batch = 1,cutoff = 1.27 * qchisq(0.95, df = d))}
      }
    )
    fitUsOnline = resUsOnline
    
    
    for(s in (1:n)){
      erreursSigmaMed2[s,m,3,sim] =norm(resUsStreaming$Sigma[s,,] - Sigma0,"F")
      if(resUsOnline$distances[s] > cutoff){outliersLabelsMed2[s,m,2,sim] = 1}
      if(resUsStreaming$distances[s] > cutoff){outliersLabelsMed2[s,m,3,sim] = 1}
    }
    
    
    if(d == 10){outliersLabelsMed2[,m,2,sim] = resUsOnline$outlier_labels}
    if(d == 100){outliersLabelsMed2[,m,2,sim] = test_outliers(distances = resUsOnline$distances,cutoff = 1.29*qchisq(.95,df = 100))}
    
    print(paste0("Erreur us online med ",erreursSigmaMed[n,m,2,sim]))
    
    #t = table(data$labelsVrais,resUsOnline$outlier_labels)
    t = table(data$labelsVrais,outliersLabelsMed2[,m,2,sim])
    if (r != 0) {faux_positifsMed2[m,2,sim] =  t[1,2]
    faux_negatifsMed2[m,2,sim] = t[2,1]}
    if(r == 0){faux_positifsMed2[m,2,sim] = t[1,2]}
    
    
    print(paste0("faux positifs med us online ",faux_positifsMed2[m,2,sim]))
    
    print(paste0("faux négatifs med us online ",faux_negatifsMed2[m,2,sim]))
    
    
    
    temps_streaming = system.time({
      if(d == 10 ){resUsStreaming= StreamingOutlierDetection(data$Z,batch = ncol(data$Z))}
      if(d == 100){
        resUsStreaming= StreamingOutlierDetection(data$Z,batch = sqrt(ncol(data$Z)),cutoff = 1.38*qchisq(0.95,df = d))}
    })
    fitUSStreaming = resUsStreaming
    for(s in (1:n)){
      erreursSigmaMed2[s,m,3,sim] =norm(resUsStreaming$Sigma[s,,] - Sigma0,"F")
      #outliersLabelsNear[s,m,3,sim] = resUsStreaming$outlier_labels[s]
    }
    if(d == 10){outliersLabelsMed2[,m,3,sim]= resUsStreaming$outlier_labels}
    if(d == 100){
      outliersLabelsMed2[,m,3,sim] = test_outliers(distances = resUsStreaming$distances,cutoff = 1.38*qchisq(.95,df = 100))
    }
    print(paste0("Erreur us streaming med ",erreursSigmaMed2[n,m,3,sim]))
    
    #t = table(data$labelsVrais,resUsStreaming$outlier_labels)
    t = table(data$labelsVrais,outliersLabelsMed2[,m,3,sim])
    if(r == 0){faux_positifsMed2[m,3,sim] =  t[1,2]}
    if (r != 0) {faux_positifsMed2[m,3,sim] =  t[1,2]
    faux_negatifsMed2[m,3,sim] = t[2,1]}
    
    
    print(paste0("faux positifs med us streaming ",faux_positifsMed2[m,3,sim]))
    
    print(paste0("faux négatifs med us streaming ",faux_negatifsMed2[m,3,sim]))
    Sigma_naive = fitNaif$SigmaIter
    Sigma_online = fitUsOnline$Sigma
    Sigma_str = fitUSStreaming$Sigma
    
    setwd(resDir)
    #save(Sigma_naive, Sigma_online, Sigma_str,temps_naif,temps_online,temps_streaming,file=fitFile)
  }
  
}


setwd("~")

file = paste0("results-k-",k,"l-",l,"rho1",rho1,".Rdata")
save(erreursSigmaMed2,faux_negatifsOracleMed2,faux_positifsOracleMed2,faux_negatifsMed2,faux_positifsMed2,outliersLabelsMed2,labelsVraisMed2,file = file)

######################################Scen Med 2######################################
#######################Erreur Med scenario###########################

k = k1val[5];l = l1valup1[5];rho1 = rho1val[5]

erreursSigmaMed3 = array(0,dim = c(n,length(rList),3,simNb))
erreursInvSigmaMed3 = array(0,dim = c(n,length(rList),3,simNb))
outliersLabelsMed3 = array(0,dim = c(n,length(rList),3,simNb))
outliersLabelsOracleMed3 = array(0,dim = c(n,length(rList),simNb))
labelsVraisMed3 = array(0,dim = c(n,length(rList)))
faux_positifsMed3= array(0,dim = c(length(rList),3,simNb))
faux_negatifsMed3 = array(0,dim = c(length(rList),3,simNb))
faux_positifsOracleMed3 = array(0,dim = c(length(rList),simNb))
faux_negatifsOracleMed3 = array(0,dim = c(length(rList),simNb))


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
    labelsVraisMed3[,m] = data$labelsVrais
    fitFile <- paste0('FitParms-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
    temps_naif = system.time(
      
      {resNaif = SampleCovOnline(data$Z)}
    )
    fitNaif = resNaif
    for(s in (1:n)){
      erreursSigmaMed3[s,m,1,sim] = norm(resNaif$SigmaIter[s,,] - Sigma0,"F")
      if(t(data$Z[s,])%*% solve(Sigma0)%*% data$Z[s,] > qchisq(.95,df = d)) {
        outliersLabelsOracleMed3[s,m,sim] = 1
      }
      
    }
    outliersLabelsMed3[,m,1,sim] = resNaif$outliers_labels
    
    print(paste0("Erreur naive med ",erreursSigmaMed2[n,m,1,sim]))
    
    t = table(data$labelsVrais,resNaif$outliers_labels)
    if (r != 0) {faux_positifsMed3[m,1,sim] =  t[1,2]
    faux_negatifsMed3[m,1,sim] = t[2,1]}
    if(r == 0){faux_positifsMed3[m,1,sim] =  t[1,2]}
    
    print(paste0("faux positifs med naive ",faux_positifsMed3[m,1,sim]))
    
    print(paste0("faux négatifs med naive ",faux_negatifsMed3[m,1,sim]))
    
    t = table(data$labelsVrais,outliersLabelsOracleMed3[,m,sim])
    if (r != 0) {faux_positifsOracleMed3[m,sim] =  t[1,2]
    faux_negatifsOracleMed3[m,sim] = t[2,1]}
    if(r == 0){faux_positifsOracleMed3[m,sim] =  t[1,2]}
    
    
    print(paste0("faux positifs Oracle ",faux_positifsOracleMed3[m,sim]))
    
    print(paste0("faux négatifs Oracle ",faux_negatifsOracleMed3[m,sim]))
    
    
    temps_online = (
      {
        if(d == 10){resUsOnline= StreamingOutlierDetection(data$Z,batch = 1)}
        if(d == 100){resUsOnline= StreamingOutlierDetection(data$Z,batch = 1,cutoff = 1.27 * qchisq(0.95, df = d))}
      }
    )
    fitUsOnline = resUsOnline
    for(s in (1:n)){
      erreursSigmaMed3[s,m,2,sim] = norm(resUsOnline$Sigma[s,,] - Sigma0,"F")
      #outliersLabelsNear[s,m,2,sim] = resUsOnline$outlier_labels[s]
    }
    
    if(d == 10){outliersLabelsMed3[,m,2,sim] = resUsOnline$outlier_labels}
    if(d == 100){outliersLabelsMed3[,m,2,sim] = test_outliers(distances = resUsOnline$distances,cutoff = 1.29*qchisq(.95,df = 100))}
    
    print(paste0("Erreur us online med ",erreursSigmaMed[n,m,2,sim]))
    
    #t = table(data$labelsVrais,resUsOnline$outlier_labels)
    t = table(data$labelsVrais,outliersLabelsMed3[,m,2,sim])
    if (r != 0) {faux_positifsMed3[m,2,sim] =  t[1,2]
    faux_negatifsMed3[m,2,sim] = t[2,1]}
    if(r == 0){faux_positifsMed3[m,2,sim] = t[1,2]}
    
    
    print(paste0("faux positifs med us online ",faux_positifsMed3[m,2,sim]))
    
    print(paste0("faux négatifs med us online ",faux_negatifsMed3[m,2,sim]))
    
    
    
    temps_streaming = system.time({
      if(d == 10 ){resUsStreaming= StreamingOutlierDetection(data$Z,batch = ncol(data$Z))}
      if(d == 100){
        resUsStreaming= StreamingOutlierDetection(data$Z,batch = sqrt(ncol(data$Z)),cutoff = 1.38*qchisq(0.95,df = d))}
    })
    fitUSStreaming = resUsStreaming
    
    # if(d == 10){outliersLabelsMed3[,m,3,sim]= resUsStreaming$outlier_labels}
    # if(d == 100){
    #   outliersLabelsMed3[,m,3,sim] = test_outliers(distances = resUsStreaming$distances,cutoff = 1.38*qchisq(.95,df = 100))
    # }
    # 
    for(s in (1:n)){
      erreursSigmaMed3[s,m,3,sim] =norm(resUsStreaming$Sigma[s,,] - Sigma0,"F")
      if(resUsOnline$distances[s] > cutoff){outliersLabelsMed3[s,m,2,sim] = 1}
      if(resUsStreaming$distances[s] > cutoff){outliersLabelsMed3[s,m,3,sim] = 1}
    }
    
    print(paste0("Erreur us streaming med ",erreursSigmaMed2[n,m,3,sim]))
    
    #t = table(data$labelsVrais,resUsStreaming$outlier_labels)
    t = table(data$labelsVrais,outliersLabelsMed3[,m,3,sim])
    if(r == 0){faux_positifsMed3[m,3,sim] =  t[1,2]}
    if (r != 0) {faux_positifsMed3[m,3,sim] =  t[1,2]
    faux_negatifsMed3[m,3,sim] = t[2,1]}
    
    
    print(paste0("faux positifs med us streaming ",faux_positifsMed3[m,3,sim]))
    
    print(paste0("faux négatifs med us streaming ",faux_negatifsMed3[m,3,sim]))
    Sigma_naive = fitNaif$SigmaIter
    Sigma_online = fitUsOnline$Sigma
    Sigma_str = fitUSStreaming$Sigma
    
    setwd(resDir)
    #save(Sigma_naive, Sigma_online, Sigma_str,temps_naif,temps_online,temps_streaming,file=fitFile)
  }
  
}


setwd("~")

file = paste0("results-k-",k,"l-",l,"rho1",rho1,".Rdata")
save(erreursSigmaMed3,faux_negatifsOracleMed3,faux_positifsOracleMed3,faux_negatifsMed3,faux_positifsMed3,outliersLabelsMed3,labelsVraisMed3,file = file)

#################Scen Med 4############################################

k = k1val[6];l = l1valup1[6];rho1 = rho1val[6]

erreursSigmaMed4 = array(0,dim = c(n,length(rList),3,simNb))
erreursInvSigmaMed4 = array(0,dim = c(n,length(rList),3,simNb))
outliersLabelsMed4 = array(0,dim = c(n,length(rList),3,simNb))
outliersLabelsOracleMed4 = array(0,dim = c(n,length(rList),simNb))
labelsVraisMed4 = array(0,dim = c(n,length(rList)))
faux_positifsMed4= array(0,dim = c(length(rList),3,simNb))
faux_negatifsMed4 = array(0,dim = c(length(rList),3,simNb))
faux_positifsOracleMed4 = array(0,dim = c(length(rList),simNb))
faux_negatifsOracleMed4 = array(0,dim = c(length(rList),simNb))


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
    labelsVraisMed4[,m] = data$labelsVrais
    fitFile <- paste0('FitParms-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
    temps_naif = system.time(
      
      {resNaif = SampleCovOnline(data$Z)}
    )
    fitNaif = resNaif
    for(s in (1:n)){
      erreursSigmaMed4[s,m,1,sim] = norm(resNaif$SigmaIter[s,,] - Sigma0,"F")
      if(t(data$Z[s,])%*% solve(Sigma0)%*% data$Z[s,] > qchisq(.95,df = d)) {
        outliersLabelsOracleMed4[s,m,sim] = 1
      }
      
    }
    outliersLabelsMed4[,m,1,sim] = resNaif$outliers_labels
    
    print(paste0("Erreur naive med ",erreursSigmaMed4[n,m,1,sim]))
    
    t = table(data$labelsVrais,resNaif$outliers_labels)
    if (r != 0) {faux_positifsMed4[m,1,sim] =  t[1,2]
    faux_negatifsMed4[m,1,sim] = t[2,1]}
    if(r == 0){faux_positifsMed4[m,1,sim] =  t[1,2]}
    
    print(paste0("faux positifs med naive ",faux_positifsMed4[m,1,sim]))
    
    print(paste0("faux négatifs med naive ",faux_negatifsMed4[m,1,sim]))
    
    t = table(data$labelsVrais,outliersLabelsOracleMed4[,m,sim])
    if (r != 0) {faux_positifsOracleMed4[m,sim] =  t[1,2]
    faux_negatifsOracleMed4[m,sim] = t[2,1]}
    if(r == 0){faux_positifsOracleMed4[m,sim] =  t[1,2]}
    
    
    print(paste0("faux positifs Oracle ",faux_positifsOracleMed4[m,sim]))
    
    print(paste0("faux négatifs Oracle ",faux_negatifsOracleMed4[m,sim]))
    
    
    temps_online = (
      {
        if(d == 10){resUsOnline= StreamingOutlierDetection(data$Z,batch = 1)}
        if(d == 100){resUsOnline= StreamingOutlierDetection(data$Z,batch = 1,cutoff = 1.27 * qchisq(0.95, df = d))}
      }
    )
    fitUsOnline = resUsOnline
    for(s in (1:n)){
      erreursSigmaMed4[s,m,2,sim] = norm(resUsOnline$Sigma[s,,] - Sigma0,"F")
      #outliersLabelsNear[s,m,2,sim] = resUsOnline$outlier_labels[s]
    }
    
    if(d == 10){outliersLabelsMed4[,m,2,sim] = resUsOnline$outlier_labels}
    if(d == 100){outliersLabelsMed4[,m,2,sim] = test_outliers(distances = resUsOnline$distances,cutoff = 1.29*qchisq(.95,df = 100))}
    
    print(paste0("Erreur us online med ",erreursSigmaMed[n,m,2,sim]))
    
    #t = table(data$labelsVrais,resUsOnline$outlier_labels)
    t = table(data$labelsVrais,outliersLabelsMed4[,m,2,sim])
    if (r != 0) {faux_positifsMed4[m,2,sim] =  t[1,2]
    faux_negatifsMed4[m,2,sim] = t[2,1]}
    if(r == 0){faux_positifsMed4[m,2,sim] = t[1,2]}
    
    
    print(paste0("faux positifs med us online ",faux_positifsMed4[m,2,sim]))
    
    print(paste0("faux négatifs med us online ",faux_negatifsMed4[m,2,sim]))
    
    
    
    temps_streaming = system.time({
      if(d == 10 ){resUsStreaming= StreamingOutlierDetection(data$Z,batch = ncol(data$Z))}
      if(d == 100){
        resUsStreaming= StreamingOutlierDetection(data$Z,batch = sqrt(ncol(data$Z)),cutoff = 1.38*qchisq(0.95,df = d))}
    })
    fitUSStreaming = resUsStreaming
    for(s in (1:n)){
      erreursSigmaMed4[s,m,3,sim] =norm(resUsStreaming$Sigma[s,,] - Sigma0,"F")
      #outliersLabelsNear[s,m,3,sim] = resUsStreaming$outlier_labels[s]
    }
    if(d == 10){outliersLabelsMed4[,m,3,sim]= resUsStreaming$outlier_labels}
    if(d == 100){
      outliersLabelsMed4[,m,3,sim] = test_outliers(distances = resUsStreaming$distances,cutoff = 1.38*qchisq(.95,df = 100))
    }
    print(paste0("Erreur us streaming med ",erreursSigmaMed4[n,m,3,sim]))
    
    #t = table(data$labelsVrais,resUsStreaming$outlier_labels)
    t = table(data$labelsVrais,outliersLabelsMed4[,m,3,sim])
    if(r == 0){faux_positifsMed4[m,3,sim] =  t[1,2]}
    if (r != 0) {faux_positifsMed4[m,3,sim] =  t[1,2]
    faux_negatifsMed4[m,3,sim] = t[2,1]}
    
    
    print(paste0("faux positifs med us streaming ",faux_positifsMed4[m,3,sim]))
    
    print(paste0("faux négatifs med us streaming ",faux_negatifsMed4[m,3,sim]))
    Sigma_naive = fitNaif$SigmaIter
    Sigma_online = fitUsOnline$Sigma
    Sigma_str = fitUSStreaming$Sigma
    
    setwd(resDir)
    #save(Sigma_naive, Sigma_online, Sigma_str,temps_naif,temps_online,temps_streaming,file=fitFile)
  }
  
}


setwd("~")

file = paste0("results-k-",k,"l-",l,"rho1",rho1,".Rdata")
save(erreursSigmaMed4,faux_negatifsOracleMed4,faux_positifsOracleMed4,faux_negatifsMed4,faux_positifsMed4,outliersLabelsMed4,labelsVraisMed4,file = file)


#################Scen Med 5############################################

k = k1val[7];l = l1valup1[7];rho1 = rho1val[7]

erreursSigmaMed5 = array(0,dim = c(n,length(rList),3,simNb))
erreursInvSigmaMed5 = array(0,dim = c(n,length(rList),3,simNb))
outliersLabelsMed5 = array(0,dim = c(n,length(rList),3,simNb))
outliersLabelsOracleMed5 = array(0,dim = c(n,length(rList),simNb))
labelsVraisMed5 = array(0,dim = c(n,length(rList)))
faux_positifsMed5= array(0,dim = c(length(rList),3,simNb))
faux_negatifsMed5 = array(0,dim = c(length(rList),3,simNb))
faux_positifsOracleMed5 = array(0,dim = c(length(rList),simNb))
faux_negatifsOracleMed5 = array(0,dim = c(length(rList),simNb))


for (m in seq_along(rList)){
  for(sim in(1:simNb)){
    r = rList[m]
    dataFile <- paste0('SimData-d', d, '-n', n, '-k', k, '-l', l, '-rho', rho1,'-r',r , '-sim', sim,".RData")
    
    print(dataFile)
    
    setwd(simDir)
    if(!file.exists(dataFile)){contParam = ParmsF1(m1, k, l, rho1)
    #save(dataFile)
    }
    data = genererEchantillon(n,n,mu1 = mu0,mu2 = contParam$mu1,Sigma1 = Sigma0,Sigma2 = contParam$Sigma1,r)
    
    #else{load(dataFile)}
    labelsVraisMed5[,m] = data$labelsVrais
    fitFile <- paste0('FitParms-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
    temps_naif = system.time(
      
      {resNaif = SampleCovOnline(data$Z)}
    )
    fitNaif = resNaif
    for(s in (1:n)){
      erreursSigmaMed5[s,m,1,sim] = norm(resNaif$SigmaIter[s,,] - Sigma0,"F")
      if(t(data$Z[s,])%*% solve(Sigma0)%*% data$Z[s,] > qchisq(.95,df = d)) {
        outliersLabelsOracleMed5[s,m,sim] = 1
      }
      
    }
    outliersLabelsMed5[,m,1,sim] = resNaif$outliers_labels
    
    print(paste0("Erreur naive med ",erreursSigmaMed5[n,m,1,sim]))
    
    t = table(data$labelsVrais,resNaif$outliers_labels)
    if (r != 0) {faux_positifsMed5[m,1,sim] =  t[1,2]
    faux_negatifsMed5[m,1,sim] = t[2,1]}
    if(r == 0){faux_positifsMed5[m,1,sim] =  t[1,2]}
    
    print(paste0("faux positifs med naive ",faux_positifsMed5[m,1,sim]))
    
    print(paste0("faux négatifs med naive ",faux_negatifsMed5[m,1,sim]))
    
    t = table(data$labelsVrais,outliersLabelsOracleMed5[,m,sim])
    if (r != 0) {faux_positifsOracleMed5[m,sim] =  t[1,2]
    faux_negatifsOracleMed5[m,sim] = t[2,1]}
    if(r == 0){faux_positifsOracleMed5[m,sim] =  t[1,2]}
    
    
    print(paste0("faux positifs Oracle ",faux_positifsOracleMed5[m,sim]))
    
    print(paste0("faux négatifs Oracle ",faux_negatifsOracleMed5[m,sim]))
    
    
    temps_online = (
      {
        if(d == 10){resUsOnline= StreamingOutlierDetection(data$Z,batch = 1)}
        if(d == 100){resUsOnline= StreamingOutlierDetection(data$Z,batch = 1,cutoff = 1.27 * qchisq(0.95, df = d))}
      }
    )
    fitUsOnline = resUsOnline
    for(s in (1:n)){
      erreursSigmaMed5[s,m,2,sim] = norm(resUsOnline$Sigma[s,,] - Sigma0,"F")
      #outliersLabelsNear[s,m,2,sim] = resUsOnline$outlier_labels[s]
    }
    
    if(d == 10){outliersLabelsMed5[,m,2,sim] = resUsOnline$outlier_labels}
    if(d == 100){outliersLabelsMed5[,m,2,sim] = test_outliers(distances = resUsOnline$distances,cutoff = 1.29*qchisq(.95,df = 100))}
    
    print(paste0("Erreur us online med ",erreursSigmaMed5[n,m,2,sim]))
    
    #t = table(data$labelsVrais,resUsOnline$outlier_labels)
    t = table(data$labelsVrais,outliersLabelsMed5[,m,2,sim])
    if (r != 0) {faux_positifsMed5[m,2,sim] =  t[1,2]
    faux_negatifsMed5[m,2,sim] = t[2,1]}
    if(r == 0){faux_positifsMed5[m,2,sim] = t[1,2]}
    
    
    print(paste0("faux positifs med us online ",faux_positifsMed5[m,2,sim]))
    
    print(paste0("faux négatifs med us online ",faux_negatifsMed5[m,2,sim]))
    
    
    
    temps_streaming = system.time({
      if(d == 10 ){resUsStreaming= StreamingOutlierDetection(data$Z,batch = ncol(data$Z))}
      if(d == 100){
        resUsStreaming= StreamingOutlierDetection(data$Z,batch = sqrt(ncol(data$Z)),cutoff = 1.38*qchisq(0.95,df = d))}
    })
    fitUSStreaming = resUsStreaming
    for(s in (1:n)){
      erreursSigmaMed5[s,m,3,sim] =norm(resUsStreaming$Sigma[s,,] - Sigma0,"F")
      #outliersLabelsNear[s,m,3,sim] = resUsStreaming$outlier_labels[s]
    }
    if(d == 10){outliersLabelsMed5[,m,3,sim]= resUsStreaming$outlier_labels}
    if(d == 100){
      outliersLabelsMed5[,m,3,sim] = test_outliers(distances = resUsStreaming$distances,cutoff = 1.38*qchisq(.95,df = 100))
    }
    print(paste0("Erreur us streaming med ",erreursSigmaMed5[n,m,3,sim]))
    
    #t = table(data$labelsVrais,resUsStreaming$outlier_labels)
    t = table(data$labelsVrais,outliersLabelsMed5[,m,3,sim])
    if(r == 0){faux_positifsMed5[m,3,sim] =  t[1,2]}
    if (r != 0) {faux_positifsMed5[m,3,sim] =  t[1,2]
    faux_negatifsMed5[m,3,sim] = t[2,1]}
    
    
    print(paste0("faux positifs med us streaming ",faux_positifsMed5[m,3,sim]))
    
    print(paste0("faux négatifs med us streaming ",faux_negatifsMed5[m,3,sim]))
    Sigma_naive = fitNaif$SigmaIter
    Sigma_online = fitUsOnline$Sigma
    Sigma_str = fitUSStreaming$Sigma
    
    setwd(resDir)
    #save(Sigma_naive, Sigma_online, Sigma_str,temps_naif,temps_online,temps_streaming,file=fitFile)
  }
  
}


setwd("~")

file = paste0("results-k-",k,"l-",l,"rho1",rho1,".Rdata")
save(erreursSigmaMed5,faux_negatifsOracleMed5,faux_positifsOracleMed5,faux_negatifsMed5,faux_positifsMed5,outliersLabelsMed5,labelsVraisMed5,file = file)

#########################Trajectoires faux négatifs, faux positifs far#############################################

faux_negatifsFarTrajMed5 = array(0,dim = c(n,length(rList),4,simNb))

faux_positifsFarTrajMed5 = array(0,dim = c(n,length(rList),4,simNb))


for(m in seq_along(rList)){
  
  for (i in (1:n))
  {
    for(l in (1:3)){
      if ( labelsVraisMed5[i,m] == 1 & outliersLabelsMed5[i,m,l,1] == 0){faux_negatifsFarTrajMed5[i,m,l,1] = 1}
      
      if (labelsVraisMed5[i,m] == 0 & outliersLabelsMed5[i,m,l,1] == 1){faux_positifsFarTrajMed5[i,m,l,1] = 1
      print(paste0("FP l =",l))}
      
    }
  }
  
  
  
}
for(m in seq_along(rList)){
  for (i in (1:n)){
    if (labelsVraisMed5[i,m] == 1 & outliersLabelsOracleMed5[i,m,1] == 0){faux_negatifsFarTrajMed5[i,m,4,1] = 1}
    
    if (labelsVraisMed5[i,m] == 0 & outliersLabelsOracleMed5[i,m,1] == 1){faux_positifsFarTrajMed5[i,m,4,1] = 1}
  }
}
faux_negatifs_cum = array(0, dim = c(n,length(rList),4))

faux_positifs_cum = array(0, dim = c(n,length(rList),4))
#Sommes cumulées
for(l in (1:4)){
  for (m in seq_along(rList)){
    faux_negatifs_cum[,m,l] <- cumsum(faux_negatifsFarTrajMed5[,m,l,1])
    faux_positifs_cum[,m,l] <- cumsum(faux_positifsFarTrajMed5[,m,l,1])}
}
taux_fn_cum = array(0,dim= c(n,length(rList),l))
taux_fp_cum = array(0,dim= c(n,length(rList),l))

for (i in (1:n))
{
  for(m in seq_along(rList)){
    n_outliers <- sum(labelsVraisMed5[1:i, m] == 1)
    n_inliers  <- sum(labelsVraisMed5[1:i, m] == 0)
    for(l in (1:4)){
      if(n_outliers > 0){
        taux_fn_cum[i,m,l] = faux_negatifs_cum[i,m,l]/n_outliers
      }
      
      
      if(n_inliers > 0){
        taux_fp_cum[i,m,l] = faux_positifs_cum[i,m,l]/n_inliers}
    }}
}


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


