simDir = "~/Simus/DataSim"

resDir = "~/Simus/FitSim"

explDir = "~/Simus/exploitResults"

##########################################
###Extraction of the parameters
##########################################

kList = k1val

lList = l1valup1

rho1List = rho1val


sim = 1
simNb = 1



k = KLval[3];l = l1val[6] ;rho1 = rho1val[1]
erreursSigmaMed5 = array(0,dim = c(n,length(rList),6,simNb))
erreursInvSigmaMed5 = array(0,dim = c(n,length(rList),6,simNb))
outliersLabelsMed5 = array(0,dim = c(n,length(rList),6,simNb))
outliersLabelsOracleMed5 = array(0,dim = c(n,length(rList),simNb))
labelsVraisMed5 = array(0,dim = c(n,length(rList)))
faux_positifsMed5= array(0,dim = c(length(rList),6,simNb))
faux_negatifsMed5 = array(0,dim = c(length(rList),6,simNb))
faux_positifsOracleMed5 = array(0,dim = c(length(rList),simNb))
faux_negatifsOracleMed5 = array(0,dim = c(length(rList),simNb))

#Pool d'outliers

id_pool <- sample(1:n)

outlier_sets <- vector("list", length(rList))

id_pool <- sample(1:n)

r_max <- max(rList)

outlier_sets = list()

for (m in seq_along(rList)) {
  n_active = floor(rList[m]/100*n)
  outlier_sets[[m]] = id_pool[1:n_active]
}

for (m in seq_along(rList)){
  for(sim in(1:simNb)){
    r = rList[m]
    dataFile <- paste0('SimData-d', d, '-n', n, '-k', k, '-l', l, '-rho', rho1,'-r',r , '-sim', sim,".RData")
    
    print(dataFile)
    
    setwd(simDir)
    if(!file.exists(dataFile)){contParam = ParmsF1(m1, k, l, rho1)
    #data = genererEchantillon(n,n,mu1 = mu0,mu2 = contParam$mu1,Sigma1 = Sigma0,Sigma2 = contParam$Sigma1,r)
    
    if(r == 0){
      
      data = genererEchantillon_new(n,d,mu1 = mu0,mu2 = contParam$mu1,Sigma1 = Sigma0,Sigma2 = contParam$Sigma1,r)
    }#save(dataFile)
    
    if(r != 0){
      
      data = genererEchantillon_new(n,d,mu1 = mu0,mu2 = contParam$mu1,Sigma1 = Sigma0,Sigma2 = contParam$Sigma1,r, id_outliers =  outlier_sets[[m]])
    }#save(dataFile)
    }
    else{load(dataFile)}
    labelsVraisMed5[,m] = data$labelsVrais
    fitFile <- paste0('FitParms-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
    temps_naif = system.time(
      
      {resNaif = SampleCovOnline(data$Z)}
    )
    fitNaif = resNaif
    erreursSigmaMed5[n,m,1,sim] = norm(resNaif$Sigma -Sigma0,"F")
    outliersLabelsMed5[,m,1,sim] = resNaif$outliers_labels
    # 
    # for(s in (1:n)){
    #   #   erreursSigmaMed3[s,m,1,sim] = norm(resNaif$SigmaIter[s,,] - Sigma0,"F")
    #   if(t(data$Z[s,])%*% solve(Sigma0)%*% data$Z[s,] > qchisq(.95,df = d)) {
    #     outliersLabelsOracleMed3[s,m,sim] = 1
    #   }
    
  }
  print(paste0("Erreur naive med ",erreursSigmaMed5[n,m,1,sim]))
  
  temps[1,sim] = temps_naif[3]
  
  
  t = table(data$labelsVrais,resNaif$outliers_labels)
  if (r != 0) {faux_positifsMed5[m,1,sim] =  t[1,2]
  faux_negatifsMed5[m,1,sim] = t[2,1]}
  if(r == 0){faux_positifsMed5[m,1,sim] =  t[1,2]}
  # 
  print(paste0("faux positifs med naive ",faux_positifsMed5[m,1,sim]))
  
  print(paste0("faux négatifs med naive ",faux_negatifsMed5[m,1,sim]))
  
  for(s in (1:n)){
    if(t(data$Z[s,])%*% solve(Sigma0)%*% data$Z[s,] > qchisq(.95,df = d)) {
      outliersLabelsOracleMed5[s,m,sim] = 1
    }
    
  }
  
  t = table(data$labelsVrais,outliersLabelsOracleMed5[,m,sim])
  if (r != 0) {faux_positifsOracleMed5[m,sim] =  t[1,2]
  faux_negatifsOracleMed5[m,sim] = t[2,1]}
  if(r == 0){faux_positifsOracleMed5[m,sim] =  t[1,2]}
  
  
  print(paste0("faux positifs Oracle ",faux_positifsOracleMed5[m,sim]))
  
  print(paste0("faux négatifs Oracle ",faux_negatifsOracleMed5[m,sim]))
  
  
  temps_online = system.time(
    {
      #if(d == 10){resUsOnline= StreamingOutlierDetection(data$Z,batch = 1)}
      #if(d == 10){
        #}
      #if(d == 100){resUsOnline= StreamingOutlierDetection(data$Z,batch = 1,cutoff = 1.27 * qchisq(0.95, df = d))}
      resUsOnline= onlineRobustVariance(data$Z,batch = 1,computeOutliers = TRUE)
      
    })
  
  temps[2,sim] = temps_online[3]
  fitUsOnline = resUsOnline
  # for(s in (1:n)){
  #   erreursSigmaMed3[s,m,2,sim] = norm(resUsOnline$Sigma[s,,] - Sigma0,"F")
  #   #outliersLabelsNear[s,m,2,sim] = resUsOnline$outlier_labels[s]
  # }
  erreursSigmaMed5[n,m,2,sim] = norm(resUsOnline$variance - Sigma0,"F")
  if(d ==10){outliersLabelsMed5[,m,2,sim] = resUsOnline$outlier_labels}
  if(d ==100){outliersLabelsMed5[,m,2,sim] = test_outliers(distances = resUsOnline$distances,cutoff = 1.29*qchisq(.95,df = 100))}
  
  print(paste0("Erreur us online med ",erreursSigmaMed5[n,m,2,sim]))
  
  #t = table(data$labelsVrais,resUsOnline$outlier_labels)
  t = table(data$labelsVrais,outliersLabelsMed5[,m,2,sim])
  
  if (nrow(t) >= 1 && ncol(t) >= 2){
  if (r != 0) {faux_positifsMed5[m,2,sim] =  t[1,2]
  faux_negatifsMed5[m,2,sim] = t[2,1]}
  if(r == 0){faux_positifsMed5[m,2,sim] = t[1,2]}
  }
  
  print(paste0("faux positifs med us online ",faux_positifsMed5[m,2,sim]))
  
  print(paste0("faux négatifs med us online ",faux_negatifsMed5[m,2,sim]))
  
  
  
  temps_streaming = system.time({
    # if(d == 10 ){
    #   #resUsStreaming= StreamingOutlierDetection(data$Z,batch = ncol(data$Z))
    # }
  #  if(d == 10 ){
      resUsStreaming= onlineRobustVariance(data$Z,computeOutliers = TRUE)
   # }
    #if(d == 100){
    #  resUsStreaming= StreamingOutlierDetection(data$Z,batch = sqrt(ncol(data$Z)))}
  })
  temps[3,sim] = temps_streaming[3]
  
  fitUSStreaming = resUsStreaming
  # for(s in (1:n)){
  #   erreursSigmaMed5[s,m,3,sim] =norm(resUsStreaming$Sigma[s,,] - Sigma0,"F")
  #   #outliersLabelsNear[s,m,3,sim] = resUsStreaming$outlier_labels[s]
  # }
  erreursSigmaMed5[n,m,3,sim] =norm(resUsStreaming$variance - Sigma0,"F")
  if(d == 10){outliersLabelsMed5[,m,3,sim]= resUsStreaming$outlier_labels}
  if(d == 100){
    outliersLabelsMed5[,m,3,sim] = test_outliers(distances = resUsStreaming$distances,cutoff = 1.38*qchisq(.95,df = 100))
  }
  print(paste0("Erreur us streaming med ",erreursSigmaMed5[n,m,3,sim]))
  
  #t = table(data$labelsVrais,resUsStreaming$outlier_labels)
  t = table(data$labelsVrais,outliersLabelsMed5[,m,3,sim])
  
  if (nrow(t) >= 1 && ncol(t) >= 2){
  if(r == 0){faux_positifsMed5[m,3,sim] =  t[1,2]}
  if (r != 0) {faux_positifsMed5[m,3,sim] =  t[1,2]
  faux_negatifsMed5[m,3,sim] = t[2,1]}
  }
  
  print(paste0("faux positifs med us streaming ",faux_positifsMed5[m,3,sim]))
  
  print(paste0("faux négatifs med us streaming ",faux_negatifsMed5[m,3,sim]))
  Sigma_naive = fitNaif$SigmaIter
  Sigma_online = fitUsOnline$Sigma
  Sigma_str = fitUSStreaming$Sigma
  
  print(paste0("Erreur mcd med  ",erreursSigmaMed5[n,m,5,sim]))
  
  
  temps_offline = system.time({
    # resOffline = OfflineOutlierD etection(data$Z)
    resOffline = offlineRobustVariance(data$Z,computeOutliers = TRUE)
  }
  )
  
  
  erreursSigmaMed5[n,m,4,sim] =norm(resOffline$variance - Sigma0,"F")
  
  print(paste0("Erreur us offline med ",erreursSigmaMed5[n,m,4,sim]))
  
  
  #outliersLabelsMed5[,m,4,sim]= resOffline$outlier_labels 
  
  # fitUsOnline = resUsOnline
  temps[4,sim] = temps_offline[3]
  
  
  print(paste0("temps offline ",temps[4,sim]))
  
  
  t = table(data$labelsVrais,resOffline$outlier_labels)
  if (r != 0) {faux_positifsMed5[m,4,sim] =  t[1,2]
  faux_negatifsMed5[m,4,sim] = t[2,1]}
  if(r == 0){faux_positifsMed5[m,4,sim] =  t[1,2]}
  
  print(paste0("faux positifs med offline ",faux_positifsMed5[m,4,sim]))
  
  print(paste0("faux négatifs med offline ",faux_negatifsMed5[m,4,sim]))
  
  
  
  temps_mcd = system.time({
    resMcd = covMcd(data$Z)
    invSigmaMCD = solve(resMcd$cov) 
    for(s in (1:n))
    {
      if (t(data$Z[s,] - resMcd$center)%*%invSigmaMCD%*%(data$Z[s,] - resMcd$center) > qchisq(.95,df = d)){outliersLabelsMed5[s,m,5,sim] = 1}
    }
    }
  )
  
  erreursSigmaMed5[n,m,5,sim] =norm(resMcd$cov - Sigma0,"F")
  
  print(paste0("Erreur mcd med ",erreursSigmaMed5[n,m,5,sim]))
  
  
  
  temps[5,sim] = temps_mcd[3]
  
  t = table(data$labelsVrais,outliersLabelsMed5[,m,5,sim] )
  if (r != 0) {faux_positifsMed5[m,5,sim] =  t[1,2]
  faux_negatifsMed5[m,5,sim] = t[2,1]}
  if(r == 0){faux_positifsMed5[m,5,sim] =  t[1,2]}
  
  print(paste0("faux positifs med MCD ",faux_positifsMed5[m,5,sim]))
  
  print(paste0("faux négatifs med MCD ",faux_negatifsMed5[m,5,sim]))
  
  
  
  print(paste0("temps MCD ",temps[5,sim]))
  
  
  temps_ogk = system.time({
    resOGK = covOGK(data$Z,sigmamu = scaleTau2)
    invSigmaOGK = solve(resOGK$cov)
    for(s in (1:n))
    {
      if (t(data$Z[s,] - resOGK$center)%*%invSigmaOGK%*%(data$Z[s,] - resOGK$center) > qchisq(.95,df = d)){outliersLabelsMed5[s,m,6,sim] = 1}
    }
    
    }
  )
  
  temps[6,sim] = temps_ogk[3]
  
  print(paste0("temps OGK ",temps[6,sim]))
  
  
  
  
  erreursSigmaMed5[n,m,6,sim] =norm(resOGK$cov - Sigma0,"F")
  
  print(paste0("Erreur OGK med ",erreursSigmaMed5[n,m,6,sim]))
  
  
  
  t = table(data$labelsVrais,outliersLabelsMed5[,m,6,sim] )
  if (r != 0) {faux_positifsMed5[m,6,sim] =  t[1,2]
  faux_negatifsMed5[m,6,sim] = t[2,1]}
  if(r == 0){faux_positifsMed5[m,6,sim] =  t[1,2]}
  
  print(paste0("faux positifs med OGK ",faux_positifsMed5[m,6,sim]))
  
  print(paste0("faux négatifs med OGK ",faux_negatifsMed5[m,6,sim]))
  
  
  
  setwd(resDir)
  #save(Sigma_naive, Sigma_online, Sigma_str,temps_naif,temps_online,temps_streaming,file=fitFile)
}



setwd("~")

file = paste0("results-k-",k,"l-",l,"rho1",rho1,"-d",d,".Rdata")
save(erreursSigmaMed5,faux_negatifsOracleMed5,faux_positifsOracleMed5,faux_negatifsMed5,faux_positifsMed5,outliersLabelsOracleMed5,outliersLabelsMed5,labelsVraisMed5,file = file)


#########################Moyenne erreurs, faux négatifs, faux positifs###############################################################
# moyenne_erreursSigmaNear = apply(erreursSigmaNear,c(1,2,3),mean)
# moyenne_erreursSigmaMed = apply(erreursSigmaMed,c(1,2,3),mean)
# moyenne_erreursSigmaMed2 = apply(erreursSigmaMed2,c(1,2,3),mean)
# moyenne_erreursSigmaMed3 = apply(erreursSigmaMed3,c(1,2,3),mean)
moyenne_erreursSigmaMed5 = apply(erreursSigmaMed5,c(1,2,3),mean)


# moyenne_faux_positifsNear = apply(faux_positifsNear,c(1,2),mean) 
# moyenne_faux_positifsMed = apply(faux_positifsMed,c(1,2),mean) 
# moyenne_faux_positifsMed2 = apply(faux_positifsMed2,c(1,2),mean) 
# moyenne_faux_positifsMed3 = apply(faux_positifsMed3,c(1,2),mean) 
moyenne_faux_positifsMed5 = apply(faux_positifsMed5,c(1,2),mean) 


# moyenne_faux_negatifsNear = apply(faux_negatifsNear,c(1,2),mean) 
# moyenne_faux_negatifsMed = apply(faux_negatifsMed,c(1,2),mean) 
# moyenne_faux_negatifsMed2 = apply(faux_negatifsMed2,c(1,2),mean) 
# moyenne_faux_negatifsMed3 = apply(faux_negatifsMed3,c(1,2),mean) 
moyenne_faux_negatifsMed5 = apply(faux_negatifsMed5,c(1,2),mean) 

# moyenne_faux_negatifsNearOracle = apply(faux_negatifsOracle,c(1,2),mean) 
# moyenne_faux_negatifsMedOracle = apply(faux_negatifsOracleMed,c(1,2),mean) 
# moyenne_faux_negatifsMed2Oracle = apply(faux_negatifsOracleMed2,c(1,2),mean) 
# moyenne_faux_negatifsMed3Oracle = apply(faux_negatifsOracleMed3,c(1,2),mean) 
moyenne_faux_negatifsMed5Oracle = apply(faux_negatifsOracleMed5,c(1,2),mean) 


# moyenne_faux_positifsNearOracle = apply(faux_positifsOracle,c(1,2),mean) 
# moyenne_faux_positifsMedOracle = apply(faux_positifsOracleMed,c(1,2),mean) 
# moyenne_faux_positifsMed2Oracle = apply(faux_positifsOracleMed2,c(1,2),mean) 
# moyenne_faux_positifsMed3Oracle = apply(faux_positifsOracleMed3,c(1,2),mean) 
moyenne_faux_negatifsMed5Oracle = apply(faux_negatifsOracleMed5,c(1,2),mean) 

############################Trajectoires##########################################################"

majority_vote_Med3 = array(0,dim = c(n,length(rList),4))
majority_vote_Med3NotOracle <- array(0, dim = c(n, length(rList), 3))
majority_vote_Med3NotOracle <- apply(outliersLabelsMed3, c(1, 2, 3), function(x) {
  ifelse(mean(x) >= 0.5, 1, 0)
})

majority_vote_Med3[,,1:3] <- majority_vote_Med3NotOracle
majority_vote_Med3Oracle <- ifelse(apply(outliersLabelsOracleMed3, c(1,2), function(x) sum(x) >= 50), 1, 0)
majority_vote_Med3[,,4] <- majority_vote_Med3Oracle
######################Trajectories#################################

faux_negatifsFarTrajMed3 = array(0,dim = c(n,length(rList),4))

faux_positifsFarTrajMed3 = array(0,dim = c(n,length(rList),4))

# 
for(m in seq_along(rList)){
  
  for (i in (1:n))
  {
    for(l in (1:3)){
      if ( labelsVraisMed3[i,m] == 1 & majority_vote_Med3[i,m,l] == 0){faux_negatifsFarTrajMed3[i,m,l] = 1}
      
      if (labelsVraisMed3[i,m] == 0 & majority_vote_Med3[i,m,l] == 1){faux_positifsFarTrajMed3[i,m,l] = 1
      print(paste0("FP l =",l))}
      
    }
  }
  
  
  
}

for(m in seq_along(rList)){
  
  for (i in (1:n))
  {
    if ( labelsVraisMed3[i,m] == 1 & outliersLabelsMed3[i,m,1,1] == 0){faux_negatifsFarTrajMed3[i,m,1] = 1}
    
    if (labelsVraisMed3[i,m] == 0 & outliersLabelsMed3[i,m,1,1] == 1){faux_positifsFarTrajMed3[i,m,1] = 1
    print(paste0("FP l =",l))}
    
  }
  
  
  
}

for(m in seq_along(rList)){
  for (i in (1:n)){
    if (labelsVraisMed3[i,m] == 1 &  majority_vote_Med3[i,m,4] == 0){faux_negatifsFarTrajMed3[i,m,4] = 1}
    
    if (labelsVraisMed3[i,m] == 0 & majority_vote_Med3[i,m,4] == 1){faux_positifsFarTrajMed3[i,m,4] = 1}
  }
}
faux_negatifs_cum = array(0, dim = c(n,length(rList),4))

faux_positifs_cum = array(0, dim = c(n,length(rList),4))
#Sommes cumulées
for(l in (1:4)){
  for (m in seq_along(rList)){
    faux_negatifs_cum[,m,l] <- cumsum(faux_negatifsFarTrajMed3[,m,l])
    faux_positifs_cum[,m,l] <- cumsum(faux_positifsFarTrajMed3[,m,l])}
}
taux_fn_cum = array(0,dim= c(n,length(rList),4))
taux_fp_cum = array(0,dim= c(n,length(rList),4))

for (i in (1:n))
{
  for(m in seq_along(rList)){
    n_outliers <- sum(labelsVraisMed3[1:i, m] == 1)
    n_inliers  <- sum(labelsVraisMed3[1:i, m] == 0)
    for(l in (1:4)){
      if(n_outliers > 0){
        taux_fn_cum[i,m,l] = faux_negatifs_cum[i,m,l]/n_outliers
      }
      
      
      if(n_inliers > 0){
        taux_fp_cum[i,m,l] = faux_positifs_cum[i,m,l]/n_inliers}
    }}
}
