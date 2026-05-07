k_values <- seq(0.1, 100, by = 0.5)



r = 30



erreursSigmaMed5 = array(0,dim = c(n,length(k_values),6,simNb))
erreursInvSigmaMed5 = array(0,dim = c(n,length(k_values),6,simNb))
outliersLabelsMed5 = array(0,dim = c(n,length(k_values),6,simNb))
outliersLabelsOracleMed5 = array(0,dim = c(n,length(k_values),simNb))
labelsVraisMed5 = array(0,dim = c(n,length(k_values)))
faux_positifsMed5= array(0,dim = c(length(k_values),6,simNb))
faux_negatifsMed5 = array(0,dim = c(length(k_values),6,simNb))
faux_positifsOracleMed5 = array(0,dim = c(length(k_values),simNb))
faux_negatifsOracleMed5 = array(0,dim = c(length(k_values),simNb))
temps = array(0,dim=c(6,simNb))


outlier_sets <- vector("list", length(rList))

id_pool <- sample(1:n)

r_max <- max(rList)

outlier_sets = list()

for (m in seq_along(rList[1:9])) {
  n_active = floor(rList[m]/100*n)
  outlier_sets[[m]] = id_pool[1:n_active]
}


m = 7

for(j in seq_along(k_values)){
  
  
  k = k_values[j]
  l = 1
  rho1 = 0.3
  
  
  dataFile <- paste0('SimData-d', d, '-n', n, '-k', k, '-l', l, '-rho', rho1,'-r',r , '-sim', sim,".RData")
  
  print(dataFile)
  
  setwd(simDir)
  contParam = ParmsF1(m1 = m1,k1 = k,l1 = l,rho1 = rho1)
  
  if(r == 0){
    
    data = genererEchantillon_new(n,d,mu1 = mu0,mu2 = contParam$mu1,Sigma1 = Sigma0,Sigma2 = contParam$Sigma1,r)
    #save(dataFile)
    
    #save(dataFile)
  }
  if(r != 0){
    
    data = genererEchantillon_new(n,d,mu1 = mu0,mu2 = contParam$mu1,Sigma1 = Sigma0,Sigma2 = contParam$Sigma1,r, id_outliers =  outlier_sets[[m]])
  }
  #else{load(dataFile)}
  labelsVraisMed5[,j] = data$labelsVrais
  
  
  
  fitFile <- paste0('FitParms-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
  temps_naif = system.time(
    
    {resNaif = SampleCovOnline(data$Z)}
  )
  fitNaif = resNaif
  erreursSigmaMed5[n,j,1,sim] = norm(resNaif$Sigma -Sigma0,"F")
  outliersLabelsMed5[,j,1,sim] = resNaif$outliers_labels
  # 
  # for(s in (1:n)){
  #   #   erreursSigmaMed3[s,m,1,sim] = norm(resNaif$SigmaIter[s,,] - Sigma0,"F")
  #   if(t(data$Z[s,])%*% solve(Sigma0)%*% data$Z[s,] > qchisq(.95,df = d)) {
  #     outliersLabelsOracleMed3[s,m,sim] = 1
  #   }
  
  
  print(paste0("Erreur naive med ",erreursSigmaMed5[n,m,1,sim]))
  
  temps[1,sim] = temps_naif[3]
  
  
  t = table(data$labelsVrais,resNaif$outliers_labels)
  if (r != 0) {faux_positifsMed5[j,1,sim] =  t[1,2]
  faux_negatifsMed5[j,1,sim] = t[2,1]}
  if(r == 0){faux_positifsMed5[j,1,sim] =  t[1,2]}
  # 
  print(paste0("faux positifs med naive ",faux_positifsMed5[j,1,sim]))
  
  print(paste0("faux négatifs med naive ",faux_negatifsMed5[j,1,sim]))
  
  for(s in (1:n)){
    if(t(data$Z[s,])%*% solve(Sigma0)%*% data$Z[s,] > qchisq(.95,df = d)) {
      outliersLabelsOracleMed5[s,j,sim] = 1
    }
    
  }
  
  t = table(data$labelsVrais,outliersLabelsOracleMed5[,j,sim])
  if (r != 0) {faux_positifsOracleMed5[j,sim] =  t[1,2]
  faux_negatifsOracleMed5[j,sim] = t[2,1]}
  if(r == 0){faux_positifsOracleMed5[j,sim] =  t[1,2]}
  
  
  print(paste0("faux positifs Oracle ",faux_positifsOracleMed5[j,sim]))
  
  print(paste0("faux négatifs Oracle ",faux_negatifsOracleMed5[j,sim]))
  
  
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
  erreursSigmaMed5[n,j,2,sim] = norm(resUsOnline$variance - Sigma0,"F")
  if(d ==10){outliersLabelsMed5[,j,2,sim] = resUsOnline$outlier_labels}
  if(d ==100){outliersLabelsMed5[,j,2,sim] = test_outliers(distances = resUsOnline$distances,cutoff = 1.29*qchisq(.95,df = 100))}
  
  print(paste0("Erreur us online med ",erreursSigmaMed5[n,m,2,sim]))
  
  #t = table(data$labelsVrais,resUsOnline$outlier_labels)
  t = table(data$labelsVrais,outliersLabelsMed5[,j,2,sim])
  
  if (nrow(t) >= 1 && ncol(t) >= 2){
    if (r != 0) {faux_positifsMed5[j,2,sim] =  t[1,2]
    faux_negatifsMed5[j,2,sim] = t[2,1]}
    if(r == 0){faux_positifsMed5[j,2,sim] = t[1,2]}
  }
  
  print(paste0("faux positifs med us online ",faux_positifsMed5[j,2,sim]))
  
  print(paste0("faux négatifs med us online ",faux_negatifsMed5[j,2,sim]))
  
  
  
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
  erreursSigmaMed5[n,j,3,sim] =norm(resUsStreaming$variance - Sigma0,"F")
  if(d == 10){outliersLabelsMed5[,j,3,sim]= resUsStreaming$outlier_labels}
  if(d == 100){
    outliersLabelsMed5[,j,3,sim] = test_outliers(distances = resUsStreaming$distances,cutoff = 1.38*qchisq(.95,df = 100))
  }
  print(paste0("Erreur us streaming med ",erreursSigmaMed5[n,j,3,sim]))
  
  #t = table(data$labelsVrais,resUsStreaming$outlier_labels)
  t = table(data$labelsVrais,outliersLabelsMed5[,j,3,sim])
  
  if (nrow(t) >= 1 && ncol(t) >= 2){
    if(r == 0){faux_positifsMed5[j,3,sim] =  t[1,2]}
    if (r != 0) {faux_positifsMed5[j,3,sim] =  t[1,2]
    faux_negatifsMed5[j,3,sim] = t[2,1]}
  }
  
  print(paste0("faux positifs med us streaming ",faux_positifsMed5[j,3,sim]))
  
  print(paste0("faux négatifs med us streaming ",faux_negatifsMed5[j,3,sim]))
  Sigma_naive = fitNaif$SigmaIter
  Sigma_online = fitUsOnline$Sigma
  Sigma_str = fitUSStreaming$Sigma
  
  print(paste0("Erreur mcd med  ",erreursSigmaMed5[n,j,5,sim]))
  
  
  temps_offline = system.time({
    # resOffline = OfflineOutlierD etection(data$Z)
    resOffline = offlineRobustVariance(data$Z,computeOutliers = TRUE)
  }
  )
  
  
  erreursSigmaMed5[n,j,4,sim] =norm(resOffline$variance - Sigma0,"F")
  
  print(paste0("Erreur us offline med ",erreursSigmaMed5[n,j,4,sim]))
  
  
  #outliersLabelsMed5[,m,4,sim]= resOffline$outlier_labels 
  
  # fitUsOnline = resUsOnline
  temps[4,sim] = temps_offline[3]
  
  
  print(paste0("temps offline ",temps[4,sim]))
  
  
  t = table(data$labelsVrais,resOffline$outlier_labels)
  if (r != 0) {faux_positifsMed5[j,4,sim] =  t[1,2]
  faux_negatifsMed5[m,4,sim] = t[2,1]}
  if(r == 0){faux_positifsMed5[j,4,sim] =  t[1,2]}
  
  print(paste0("faux positifs med offline ",faux_positifsMed5[j,4,sim]))
  
  print(paste0("faux négatifs med offline ",faux_negatifsMed5[j,4,sim]))
  
  
  
  temps_mcd = system.time({
    resMcd = covMcd(data$Z)
    invSigmaMCD = solve(resMcd$cov) 
    for(s in (1:n))
    {
      if (t(data$Z[s,] - resMcd$center)%*%invSigmaMCD%*%(data$Z[s,] - resMcd$center) > qchisq(.95,df = d)){outliersLabelsMed5[s,j,5,sim] = 1}
    }
  }
  )
  
  erreursSigmaMed5[n,j,5,sim] =norm(resMcd$cov - Sigma0,"F")
  
  print(paste0("Erreur mcd med ",erreursSigmaMed5[n,j,5,sim]))
  
  
  
  temps[5,sim] = temps_mcd[3]
  
  t = table(data$labelsVrais,outliersLabelsMed5[,j,5,sim] )
  if (r != 0) {faux_positifsMed5[j,5,sim] =  t[1,2]
  faux_negatifsMed5[j,5,sim] = t[2,1]}
  if(r == 0){faux_positifsMed5[j,5,sim] =  t[1,2]}
  
  print(paste0("faux positifs med MCD ",faux_positifsMed5[j,5,sim]))
  
  print(paste0("faux négatifs med MCD ",faux_negatifsMed5[j,5,sim]))
  
  
  
  print(paste0("temps MCD ",temps[5,sim]))
  
  
  temps_ogk = system.time({
    resOGK = covOGK(data$Z,sigmamu = scaleTau2)
    invSigmaOGK = solve(resOGK$cov)
    for(s in (1:n))
    {
      if (t(data$Z[s,] - resOGK$center)%*%invSigmaOGK%*%(data$Z[s,] - resOGK$center) > qchisq(.95,df = d)){outliersLabelsMed5[s,j,6,sim] = 1}
    }
    
  }
  )
  
  temps[6,sim] = temps_ogk[3]
  
  print(paste0("temps OGK ",temps[6,sim]))
  
  
  
  
  erreursSigmaMed5[n,j,6,sim] =norm(resOGK$cov - Sigma0,"F")
  
  print(paste0("Erreur OGK med ",erreursSigmaMed5[n,j,6,sim]))
  
  
  
  t = table(data$labelsVrais,outliersLabelsMed5[,j,6,sim] )
  if (r != 0) {faux_positifsMed5[j,6,sim] =  t[1,2]
  faux_negatifsMed5[j,6,sim] = t[2,1]}
  if(r == 0){faux_positifsMed5[j,6,sim] =  t[1,2]}
  
  print(paste0("faux positifs med OGK ",faux_positifsMed5[j,6,sim]))
  
  print(paste0("faux négatifs med OGK ",faux_negatifsMed5[j,6,sim]))
  
  
  
  setwd(resDir)
  #save(Sigma_naive, Sigma_online, Sigma_str,temps_naif,temps_online,temps_streaming,file=fitFile)
  
  
}

setwd("~")


# --- Transformation pseudo-log  ---
pseudo_log <- function(y) {
  return(log10(1+y))  
}


file <- paste0("scen-k", ".pdf")

pdf(file, width = 18, height = 6)  # largeur augmentée pour 3 plots
par(mfrow = c(1, 3), mar = c(5, 5, 2, 1))  # 1 ligne, 3 colonnes

### --- 1. Sigma erreurs (log Y) ---
plot(k_values, erreursSigmaMed5[n,, 3,1],
     type = "l", lwd = 4, col = "red",
     xlab = "", ylab = "",
     yaxt = "n", xaxt = "n",
     log = "y",
     ylim = c(1e-1, 1e2)
)

lines(k_values, erreursSigmaMed5[n,,2,1], lwd = 4, col = "blue", lty = "dashed")
lines(k_values, erreursSigmaMed5[n,,1,1], lwd = 4, col = "darkgreen", lty = "dotted")
lines(k_values, erreursSigmaMed5[n,,4,1], lwd = 4, col = "orange", lty = "dotted")
lines(k_values, erreursSigmaMed5[n,,5,1], lwd = 5, col = "black", lty = "dotted")
lines(k_values, erreursSigmaMed5[n,,6,1], lwd = 5, col = "brown", lty = "dotted")

log_ticks <- 10^seq(-1, 2, by = 1)
axis(2, at = log_ticks,
     labels = parse(text = paste0("10^", -1:2)),
     las = 1, cex.axis = 1.5)

axis(1, at = k_values, las = 1, cex.axis = 1.5)
box()

### --- 2. False negatives ---
y_red   <- faux_negatifsMed5[,3,1]/((30/100)*n)*100
y_blue  <- faux_negatifsMed5[,2,1]/((30/100)*n)*100
y_green <- faux_negatifsMed5[,1,1]/((30/100)*n)*100
#y_purp  <- faux_negatifsMed5Oracle[,1,1]/((30/100)*n)*100
y_orange <- faux_negatifsMed5[,4,1]/((30/100)*n)*100
y_grey <- faux_negatifsMed5[,5,1]/((30/100)*n)*100
y_brown <- faux_negatifsMed5[,6,1]/((30/100)*n)*100

yticks <- c(0, 1, 10, 100)

plot(k_values, pseudo_log(y_red),
     type = "l", lwd = 4, col = "red",
     xlab = "", ylab = "",
     yaxt = "n", xaxt = "n",
     ylim = pseudo_log(c(0, 100))
)

lines(k_values, pseudo_log(y_blue),  lwd = 4, col = "blue",      lty = "dashed")
lines(k_values, pseudo_log(y_green), lwd = 4, col = "darkgreen", lty = "dotted")
#lines(k_values, pseudo_log(y_purp),  lwd = 4, col = "purple4",   lty = "longdash")
lines(k_values, pseudo_log(y_orange),lwd = 4, col = "orange",    lty = "longdash")
lines(k_values, pseudo_log(y_grey),  lwd = 4, col = "black",     lty = "longdash")
lines(k_values, pseudo_log(y_brown), lwd = 4, col = "brown",     lty = "longdash")

axis(1, at = k_values, las = 1, cex.axis = 1.5)
axis(2, at = pseudo_log(yticks),
     labels = c("0", expression(10^0), expression(10^1), expression(10^2)),
     las = 1, cex.axis = 1.5)
box()


### --- 3. False positives ---
plot(k_values, faux_positifsMed5[,3,1]/((1 - 30/100)*n)*100,
     type = "l", lwd = 4, col = "red",
     xlab = "", ylab = "",
     yaxt = "n", xaxt = "n",
     ylim = c(0, 20)
)

lines(k_values, faux_positifsMed5[,2,1]/((1 - 30/100)*n)*100, lwd = 4, col = "blue", lty = "dashed")
lines(k_values, faux_positifsMed5[,1,1]/((1 - 30/100)*n)*100, lwd = 4, col = "darkgreen", lty = "dotted")
#lines(k_values, faux_positifsMed5Oracle[1:9,1]/((1 - rList[1:9]/100)*n)*100, lwd = 4, col = "purple", lty = "longdash")
lines(k_values, faux_positifsMed5[,4,1]/((1 - 30/100)*n)*100, lwd = 4, col = "orange", lty = "longdash")
lines(k_values, faux_positifsMed5[,5,1]/((1 - 30/100)*n)*100, lwd = 4, col = "black", lty = "longdash")
lines(k_values, faux_positifsMed5[,6,1]/((1 - 30/100)*n)*100, lwd = 4, col = "brown", lty = "longdash")

axis(1, at = k_values, las = 1, cex.axis = 1.5)
axis(2, at = seq(0, 20, by = 5), las = 1, cex.axis = 1.5)
box()


dev.off()




linf1_values <- seq(0.01, 0.96, by = 0.05)



r = 30



erreursSigmaMed5 = array(0,dim = c(n,length(linf1_values),6,simNb))
erreursInvSigmaMed5 = array(0,dim = c(n,length(linf1_values),6,simNb))
outliersLabelsMed5 = array(0,dim = c(n,length(linf1_values),6,simNb))
outliersLabelsOracleMed5 = array(0,dim = c(n,length(linf1_values),simNb))
labelsVraisMed5 = array(0,dim = c(n,length(linf1_values)))
faux_positifsMed5= array(0,dim = c(length(linf1_values),6,simNb))
faux_negatifsMed5 = array(0,dim = c(length(linf1_values),6,simNb))
faux_positifsOracleMed5 = array(0,dim = c(length(linf1_values),simNb))
faux_negatifsOracleMed5 = array(0,dim = c(length(linf1_values),simNb))
temps = array(0,dim=c(6,simNb))


outlier_sets <- vector("list", length(rList))

id_pool <- sample(1:n)

r_max <- max(rList)

outlier_sets = list()

for (m in seq_along(rList[1:9])) {
  n_active = floor(rList[m]/100*n)
  outlier_sets[[m]] = id_pool[1:n_active]
}


m = 7

for(j in seq_along(linf1_values)){
  
  
  k = 0
  l = linf1_values[j]
  rho1 = 0.3
  
  
  dataFile <- paste0('SimData-d', d, '-n', n, '-k', k, '-l', l, '-rho', rho1,'-r',r , '-sim', sim,".RData")
  
  print(dataFile)
  
  setwd(simDir)
  contParam = ParmsF1(m1 = m1,k1 = k,l1 = l,rho1 = rho1)
  
  if(r == 0){
    
    data = genererEchantillon_new(n,d,mu1 = mu0,mu2 = contParam$mu1,Sigma1 = Sigma0,Sigma2 = contParam$Sigma1,r)
    #save(dataFile)
    
    #save(dataFile)
  }
  if(r != 0){
    
    data = genererEchantillon_new(n,d,mu1 = mu0,mu2 = contParam$mu1,Sigma1 = Sigma0,Sigma2 = contParam$Sigma1,r, id_outliers =  outlier_sets[[m]])
  }
  #else{load(dataFile)}
  labelsVraisMed5[,j] = data$labelsVrais
  
  
  
  fitFile <- paste0('FitParms-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
  temps_naif = system.time(
    
    {resNaif = SampleCovOnline(data$Z)}
  )
  fitNaif = resNaif
  erreursSigmaMed5[n,j,1,sim] = norm(resNaif$Sigma -Sigma0,"F")
  outliersLabelsMed5[,j,1,sim] = resNaif$outliers_labels
  # 
  # for(s in (1:n)){
  #   #   erreursSigmaMed3[s,m,1,sim] = norm(resNaif$SigmaIter[s,,] - Sigma0,"F")
  #   if(t(data$Z[s,])%*% solve(Sigma0)%*% data$Z[s,] > qchisq(.95,df = d)) {
  #     outliersLabelsOracleMed3[s,m,sim] = 1
  #   }
  
  
  print(paste0("Erreur naive med ",erreursSigmaMed5[n,m,1,sim]))
  
  temps[1,sim] = temps_naif[3]
  
  
  t = table(data$labelsVrais,resNaif$outliers_labels)
  if (r != 0) {faux_positifsMed5[j,1,sim] =  t[1,2]
  faux_negatifsMed5[j,1,sim] = t[2,1]}
  if(r == 0){faux_positifsMed5[j,1,sim] =  t[1,2]}
  # 
  print(paste0("faux positifs med naive ",faux_positifsMed5[j,1,sim]))
  
  print(paste0("faux négatifs med naive ",faux_negatifsMed5[j,1,sim]))
  
  for(s in (1:n)){
    if(t(data$Z[s,])%*% solve(Sigma0)%*% data$Z[s,] > qchisq(.95,df = d)) {
      outliersLabelsOracleMed5[s,j,sim] = 1
    }
    
  }
  
  t = table(data$labelsVrais,outliersLabelsOracleMed5[,j,sim])
  if (r != 0) {faux_positifsOracleMed5[j,sim] =  t[1,2]
  faux_negatifsOracleMed5[j,sim] = t[2,1]}
  if(r == 0){faux_positifsOracleMed5[j,sim] =  t[1,2]}
  
  
  print(paste0("faux positifs Oracle ",faux_positifsOracleMed5[j,sim]))
  
  print(paste0("faux négatifs Oracle ",faux_negatifsOracleMed5[j,sim]))
  
  
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
  erreursSigmaMed5[n,j,2,sim] = norm(resUsOnline$variance - Sigma0,"F")
  if(d ==10){outliersLabelsMed5[,j,2,sim] = resUsOnline$outlier_labels}
  if(d ==100){outliersLabelsMed5[,j,2,sim] = test_outliers(distances = resUsOnline$distances,cutoff = 1.29*qchisq(.95,df = 100))}
  
  print(paste0("Erreur us online med ",erreursSigmaMed5[n,m,2,sim]))
  
  #t = table(data$labelsVrais,resUsOnline$outlier_labels)
  t = table(data$labelsVrais,outliersLabelsMed5[,j,2,sim])
  
  if (nrow(t) >= 1 && ncol(t) >= 2){
    if (r != 0) {faux_positifsMed5[j,2,sim] =  t[1,2]
    faux_negatifsMed5[j,2,sim] = t[2,1]}
    if(r == 0){faux_positifsMed5[j,2,sim] = t[1,2]}
  }
  
  print(paste0("faux positifs med us online ",faux_positifsMed5[j,2,sim]))
  
  print(paste0("faux négatifs med us online ",faux_negatifsMed5[j,2,sim]))
  
  
  
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
  erreursSigmaMed5[n,j,3,sim] =norm(resUsStreaming$variance - Sigma0,"F")
  if(d == 10){outliersLabelsMed5[,j,3,sim]= resUsStreaming$outlier_labels}
  if(d == 100){
    outliersLabelsMed5[,j,3,sim] = test_outliers(distances = resUsStreaming$distances,cutoff = 1.38*qchisq(.95,df = 100))
  }
  print(paste0("Erreur us streaming med ",erreursSigmaMed5[n,j,3,sim]))
  
  #t = table(data$labelsVrais,resUsStreaming$outlier_labels)
  t = table(data$labelsVrais,outliersLabelsMed5[,j,3,sim])
  
  if (nrow(t) >= 1 && ncol(t) >= 2){
    if(r == 0){faux_positifsMed5[j,3,sim] =  t[1,2]}
    if (r != 0) {faux_positifsMed5[j,3,sim] =  t[1,2]
    faux_negatifsMed5[j,3,sim] = t[2,1]}
  }
  
  print(paste0("faux positifs med us streaming ",faux_positifsMed5[j,3,sim]))
  
  print(paste0("faux négatifs med us streaming ",faux_negatifsMed5[j,3,sim]))
  Sigma_naive = fitNaif$SigmaIter
  Sigma_online = fitUsOnline$Sigma
  Sigma_str = fitUSStreaming$Sigma
  
  print(paste0("Erreur mcd med  ",erreursSigmaMed5[n,j,5,sim]))
  
  
  temps_offline = system.time({
    # resOffline = OfflineOutlierD etection(data$Z)
    resOffline = offlineRobustVariance(data$Z,computeOutliers = TRUE)
  }
  )
  
  
  erreursSigmaMed5[n,j,4,sim] =norm(resOffline$variance - Sigma0,"F")
  
  print(paste0("Erreur us offline med ",erreursSigmaMed5[n,j,4,sim]))
  
  
  #outliersLabelsMed5[,m,4,sim]= resOffline$outlier_labels 
  
  # fitUsOnline = resUsOnline
  temps[4,sim] = temps_offline[3]
  
  
  print(paste0("temps offline ",temps[4,sim]))
  
  
  t = table(data$labelsVrais,resOffline$outlier_labels)
  if (r != 0) {faux_positifsMed5[j,4,sim] =  t[1,2]
  faux_negatifsMed5[j,4,sim] = t[2,1]}
  if(r == 0){faux_positifsMed5[j,4,sim] =  t[1,2]}
  
  print(paste0("faux positifs med offline ",faux_positifsMed5[j,4,sim]))
  
  print(paste0("faux négatifs med offline ",faux_negatifsMed5[j,4,sim]))
  
  
  
  temps_mcd = system.time({
    resMcd = covMcd(data$Z)
    invSigmaMCD = solve(resMcd$cov) 
    for(s in (1:n))
    {
      if (t(data$Z[s,] - resMcd$center)%*%invSigmaMCD%*%(data$Z[s,] - resMcd$center) > qchisq(.95,df = d)){outliersLabelsMed5[s,j,5,sim] = 1}
    }
  }
  )
  
  erreursSigmaMed5[n,j,5,sim] =norm(resMcd$cov - Sigma0,"F")
  
  print(paste0("Erreur mcd med ",erreursSigmaMed5[n,j,5,sim]))
  
  
  
  temps[5,sim] = temps_mcd[3]
  
  t = table(data$labelsVrais,outliersLabelsMed5[,j,5,sim] )
  if (r != 0) {faux_positifsMed5[j,5,sim] =  t[1,2]
  faux_negatifsMed5[j,5,sim] = t[2,1]}
  if(r == 0){faux_positifsMed5[j,5,sim] =  t[1,2]}
  
  print(paste0("faux positifs med MCD ",faux_positifsMed5[j,5,sim]))
  
  print(paste0("faux négatifs med MCD ",faux_negatifsMed5[j,5,sim]))
  
  
  
  print(paste0("temps MCD ",temps[5,sim]))
  
  
  temps_ogk = system.time({
    resOGK = covOGK(data$Z,sigmamu = scaleTau2)
    invSigmaOGK = solve(resOGK$cov)
    for(s in (1:n))
    {
      if (t(data$Z[s,] - resOGK$center)%*%invSigmaOGK%*%(data$Z[s,] - resOGK$center) > qchisq(.95,df = d)){outliersLabelsMed5[s,j,6,sim] = 1}
    }
    
  }
  )
  
  temps[6,sim] = temps_ogk[3]
  
  print(paste0("temps OGK ",temps[6,sim]))
  
  
  
  
  erreursSigmaMed5[n,j,6,sim] =norm(resOGK$cov - Sigma0,"F")
  
  print(paste0("Erreur OGK med ",erreursSigmaMed5[n,j,6,sim]))
  
  
  
  t = table(data$labelsVrais,outliersLabelsMed5[,j,6,sim] )
  if (r != 0) {faux_positifsMed5[j,6,sim] =  t[1,2]
  faux_negatifsMed5[j,6,sim] = t[2,1]}
  if(r == 0){faux_positifsMed5[j,6,sim] =  t[1,2]}
  
  print(paste0("faux positifs med OGK ",faux_positifsMed5[j,6,sim]))
  
  print(paste0("faux négatifs med OGK ",faux_negatifsMed5[j,6,sim]))
  
  
  
  setwd(resDir)
  #save(Sigma_naive, Sigma_online, Sigma_str,temps_naif,temps_online,temps_streaming,file=fitFile)
  
  
}
setwd("~")
file <- paste0("scen-linf1", ".pdf")

pdf(file, width = 18, height = 6)  # largeur augmentée pour 3 plots
par(mfrow = c(1, 3), mar = c(5, 5, 2, 1))  # 1 ligne, 3 colonnes

### --- 1. Sigma erreurs (log Y) ---
plot(linf1_values[1:19], erreursSigmaMed5[n,1:19, 3,1],
     type = "l", lwd = 4, col = "red",
     xlab = "", ylab = "",
     yaxt = "n", xaxt = "n",
     log = "y",
     ylim = c(1e-1, 1e2)
)

lines(linf1_values[1:19], erreursSigmaMed5[n,1:19,2,1], lwd = 4, col = "blue", lty = "dashed")
lines(linf1_values[1:19], erreursSigmaMed5[n,1:19,1,1], lwd = 4, col = "darkgreen", lty = "dotted")
lines(linf1_values[1:19], erreursSigmaMed5[n,1:19,4,1], lwd = 4, col = "orange", lty = "dotted")
lines(linf1_values[1:19], erreursSigmaMed5[n,1:19,5,1], lwd = 5, col = "black", lty = "dotted")
lines(linf1_values[1:19], erreursSigmaMed5[n,1:19,6,1], lwd = 5, col = "brown", lty = "dotted")

log_ticks <- 10^seq(-1, 2, by = 1)
axis(2, at = log_ticks,
     labels = parse(text = paste0("10^", -1:2)),
     las = 1, cex.axis = 1.5)

axis(1, at = linf1_values, las = 1, cex.axis = 1.5)
box()

### --- 2. False negatives ---
y_red   <- faux_negatifsMed5[1:19,3,1]/((30/100)*n)*100
y_blue  <- faux_negatifsMed5[1:19,2,1]/((30/100)*n)*100
y_green <- faux_negatifsMed5[1:19,1,1]/((30/100)*n)*100
#y_purp  <- faux_negatifsMed5Oracle[,1,1]/((30/100)*n)*100
y_orange <- faux_negatifsMed5[1:19,4,1]/((30/100)*n)*100
y_grey <- faux_negatifsMed5[1:19,5,1]/((30/100)*n)*100
y_brown <- faux_negatifsMed5[1:19,6,1]/((30/100)*n)*100

yticks <- c(0, 1, 10, 100)

plot(linf1_values[1:19], pseudo_log(y_red),
     type = "l", lwd = 4, col = "red",
     xlab = "", ylab = "",
     yaxt = "n", xaxt = "n",
     ylim = pseudo_log(c(0, 100))
)

lines(linf1_values[1:19], pseudo_log(y_blue),  lwd = 4, col = "blue",      lty = "dashed")
lines(linf1_values[1:19], pseudo_log(y_green), lwd = 4, col = "darkgreen", lty = "dotted")
#lines(k_values, pseudo_log(y_purp),  lwd = 4, col = "purple4",   lty = "longdash")
lines(linf1_values[1:19], pseudo_log(y_orange),lwd = 4, col = "orange",    lty = "longdash")
lines(linf1_values[1:19], pseudo_log(y_grey),  lwd = 4, col = "black",     lty = "longdash")
lines(linf1_values[1:19], pseudo_log(y_brown), lwd = 4, col = "brown",     lty = "longdash")

axis(1, at = linf1_values[1:19], las = 1, cex.axis = 1.5)
axis(2, at = pseudo_log(yticks),
     labels = c("0", expression(10^0), expression(10^1), expression(10^2)),
     las = 1, cex.axis = 1.5)
box()


### --- 3. False positives ---
plot(linf1_values[1:19], faux_positifsMed5[1:19,3,1]/((1 - 30/100)*n)*100,
     type = "l", lwd = 4, col = "red",
     xlab = "", ylab = "",
     yaxt = "n", xaxt = "n",
     ylim = c(0, 20)
)

lines(linf1_values[1:19], faux_positifsMed5[1:19,2,1]/((1 - 30/100)*n)*100, lwd = 4, col = "blue", lty = "dashed")
lines(linf1_values[1:19], faux_positifsMed5[1:19,1,1]/((1 - 30/100)*n)*100, lwd = 4, col = "darkgreen", lty = "dotted")
#lines(k_values, faux_positifsMed5Oracle[1:9,1]/((1 - rList[1:9]/100)*n)*100, lwd = 4, col = "purple", lty = "longdash")
lines(linf1_values[1:19], faux_positifsMed5[1:19,4,1]/((1 - 30/100)*n)*100, lwd = 4, col = "orange", lty = "longdash")
lines(linf1_values[1:19], faux_positifsMed5[1:19,5,1]/((1 - 30/100)*n)*100, lwd = 4, col = "black", lty = "longdash")
lines(linf1_values[1:19], faux_positifsMed5[1:19,6,1]/((1 - 30/100)*n)*100, lwd = 4, col = "brown", lty = "longdash")

axis(1, at = linf1_values[1:19], las = 1, cex.axis = 1.5)
axis(2, at = seq(0, 20, by = 5), las = 1, cex.axis = 1.5)
box()


dev.off()


lup1_values <- seq(1, 1e3, by = 2)



r = 30



erreursSigmaMed5 = array(0,dim = c(n,length(lup1_values),6,simNb))
erreursInvSigmaMed5 = array(0,dim = c(n,length(lup1_values),6,simNb))
outliersLabelsMed5 = array(0,dim = c(n,length(lup1_values),6,simNb))
outliersLabelsOracleMed5 = array(0,dim = c(n,length(lup1_values),simNb))
labelsVraisMed5 = array(0,dim = c(n,length(lup1_values)))
faux_positifsMed5= array(0,dim = c(length(lup1_values),6,simNb))
faux_negatifsMed5 = array(0,dim = c(length(lup1_values),6,simNb))
faux_positifsOracleMed5 = array(0,dim = c(length(lup1_values),simNb))
faux_negatifsOracleMed5 = array(0,dim = c(length(lup1_values),simNb))
temps = array(0,dim=c(6,simNb))


outlier_sets <- vector("list", length(rList))

id_pool <- sample(1:n)

r_max <- max(rList)

outlier_sets = list()

for (m in seq_along(rList[1:9])) {
  n_active = floor(rList[m]/100*n)
  outlier_sets[[m]] = id_pool[1:n_active]
}


m = 7

for(j in seq_along(lup1_values)){
  
  
  k = 0
  l = linf1_values[j]
  rho1 = 0.3
  
  
  dataFile <- paste0('SimData-d', d, '-n', n, '-k', k, '-l', l, '-rho', rho1,'-r',r , '-sim', sim,".RData")
  
  print(dataFile)
  
  setwd(simDir)
  contParam = ParmsF1(m1 = m1,k1 = k,l1 = l,rho1 = rho1)
  
  if(r == 0){
    
    data = genererEchantillon_new(n,d,mu1 = mu0,mu2 = contParam$mu1,Sigma1 = Sigma0,Sigma2 = contParam$Sigma1,r)
    #save(dataFile)
    
    #save(dataFile)
  }
  if(r != 0){
    
    data = genererEchantillon_new(n,d,mu1 = mu0,mu2 = contParam$mu1,Sigma1 = Sigma0,Sigma2 = contParam$Sigma1,r, id_outliers =  outlier_sets[[m]])
  }
  #else{load(dataFile)}
  labelsVraisMed5[,j] = data$labelsVrais
  
  
  
  fitFile <- paste0('FitParms-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
  temps_naif = system.time(
    
    {resNaif = SampleCovOnline(data$Z)}
  )
  fitNaif = resNaif
  erreursSigmaMed5[n,j,1,sim] = norm(resNaif$Sigma -Sigma0,"F")
  outliersLabelsMed5[,j,1,sim] = resNaif$outliers_labels
  # 
  # for(s in (1:n)){
  #   #   erreursSigmaMed3[s,m,1,sim] = norm(resNaif$SigmaIter[s,,] - Sigma0,"F")
  #   if(t(data$Z[s,])%*% solve(Sigma0)%*% data$Z[s,] > qchisq(.95,df = d)) {
  #     outliersLabelsOracleMed3[s,m,sim] = 1
  #   }
  
  
  print(paste0("Erreur naive med ",erreursSigmaMed5[n,m,1,sim]))
  
  temps[1,sim] = temps_naif[3]
  
  
  t = table(data$labelsVrais,resNaif$outliers_labels)
  if (r != 0) {faux_positifsMed5[j,1,sim] =  t[1,2]
  faux_negatifsMed5[j,1,sim] = t[2,1]}
  if(r == 0){faux_positifsMed5[j,1,sim] =  t[1,2]}
  # 
  print(paste0("faux positifs med naive ",faux_positifsMed5[j,1,sim]))
  
  print(paste0("faux négatifs med naive ",faux_negatifsMed5[j,1,sim]))
  
  for(s in (1:n)){
    if(t(data$Z[s,])%*% solve(Sigma0)%*% data$Z[s,] > qchisq(.95,df = d)) {
      outliersLabelsOracleMed5[s,j,sim] = 1
    }
    
  }
  
  t = table(data$labelsVrais,outliersLabelsOracleMed5[,j,sim])
  if (r != 0) {faux_positifsOracleMed5[j,sim] =  t[1,2]
  faux_negatifsOracleMed5[j,sim] = t[2,1]}
  if(r == 0){faux_positifsOracleMed5[j,sim] =  t[1,2]}
  
  
  print(paste0("faux positifs Oracle ",faux_positifsOracleMed5[j,sim]))
  
  print(paste0("faux négatifs Oracle ",faux_negatifsOracleMed5[j,sim]))
  
  
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
  erreursSigmaMed5[n,j,2,sim] = norm(resUsOnline$variance - Sigma0,"F")
  if(d ==10){outliersLabelsMed5[,j,2,sim] = resUsOnline$outlier_labels}
  if(d ==100){outliersLabelsMed5[,j,2,sim] = test_outliers(distances = resUsOnline$distances,cutoff = 1.29*qchisq(.95,df = 100))}
  
  print(paste0("Erreur us online med ",erreursSigmaMed5[n,m,2,sim]))
  
  #t = table(data$labelsVrais,resUsOnline$outlier_labels)
  t = table(data$labelsVrais,outliersLabelsMed5[,j,2,sim])
  
  if (nrow(t) >= 1 && ncol(t) >= 2){
    if (r != 0) {faux_positifsMed5[j,2,sim] =  t[1,2]
    faux_negatifsMed5[j,2,sim] = t[2,1]}
    if(r == 0){faux_positifsMed5[j,2,sim] = t[1,2]}
  }
  
  print(paste0("faux positifs med us online ",faux_positifsMed5[j,2,sim]))
  
  print(paste0("faux négatifs med us online ",faux_negatifsMed5[j,2,sim]))
  
  
  
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
  erreursSigmaMed5[n,j,3,sim] =norm(resUsStreaming$variance - Sigma0,"F")
  if(d == 10){outliersLabelsMed5[,j,3,sim]= resUsStreaming$outlier_labels}
  if(d == 100){
    outliersLabelsMed5[,j,3,sim] = test_outliers(distances = resUsStreaming$distances,cutoff = 1.38*qchisq(.95,df = 100))
  }
  print(paste0("Erreur us streaming med ",erreursSigmaMed5[n,j,3,sim]))
  
  #t = table(data$labelsVrais,resUsStreaming$outlier_labels)
  t = table(data$labelsVrais,outliersLabelsMed5[,j,3,sim])
  
  if (nrow(t) >= 1 && ncol(t) >= 2){
    if(r == 0){faux_positifsMed5[j,3,sim] =  t[1,2]}
    if (r != 0) {faux_positifsMed5[j,3,sim] =  t[1,2]
    faux_negatifsMed5[j,3,sim] = t[2,1]}
  }
  
  print(paste0("faux positifs med us streaming ",faux_positifsMed5[j,3,sim]))
  
  print(paste0("faux négatifs med us streaming ",faux_negatifsMed5[j,3,sim]))
  Sigma_naive = fitNaif$SigmaIter
  Sigma_online = fitUsOnline$Sigma
  Sigma_str = fitUSStreaming$Sigma
  
  print(paste0("Erreur mcd med  ",erreursSigmaMed5[n,j,5,sim]))
  
  
  temps_offline = system.time({
    # resOffline = OfflineOutlierD etection(data$Z)
    resOffline = offlineRobustVariance(data$Z,computeOutliers = TRUE)
  }
  )
  
  
  erreursSigmaMed5[n,j,4,sim] =norm(resOffline$variance - Sigma0,"F")
  
  print(paste0("Erreur us offline med ",erreursSigmaMed5[n,j,4,sim]))
  
  
  #outliersLabelsMed5[,m,4,sim]= resOffline$outlier_labels 
  
  # fitUsOnline = resUsOnline
  temps[4,sim] = temps_offline[3]
  
  
  print(paste0("temps offline ",temps[4,sim]))
  
  
  t = table(data$labelsVrais,resOffline$outlier_labels)
  if (r != 0) {faux_positifsMed5[j,4,sim] =  t[1,2]
  faux_negatifsMed5[j,4,sim] = t[2,1]}
  if(r == 0){faux_positifsMed5[j,4,sim] =  t[1,2]}
  
  print(paste0("faux positifs med offline ",faux_positifsMed5[j,4,sim]))
  
  print(paste0("faux négatifs med offline ",faux_negatifsMed5[j,4,sim]))
  
  
  
  temps_mcd = system.time({
    resMcd = covMcd(data$Z)
    invSigmaMCD = solve(resMcd$cov) 
    for(s in (1:n))
    {
      if (t(data$Z[s,] - resMcd$center)%*%invSigmaMCD%*%(data$Z[s,] - resMcd$center) > qchisq(.95,df = d)){outliersLabelsMed5[s,j,5,sim] = 1}
    }
  }
  )
  
  erreursSigmaMed5[n,j,5,sim] =norm(resMcd$cov - Sigma0,"F")
  
  print(paste0("Erreur mcd med ",erreursSigmaMed5[n,j,5,sim]))
  
  
  
  temps[5,sim] = temps_mcd[3]
  
  t = table(data$labelsVrais,outliersLabelsMed5[,j,5,sim] )
  if (r != 0) {faux_positifsMed5[j,5,sim] =  t[1,2]
  faux_negatifsMed5[j,5,sim] = t[2,1]}
  if(r == 0){faux_positifsMed5[j,5,sim] =  t[1,2]}
  
  print(paste0("faux positifs med MCD ",faux_positifsMed5[j,5,sim]))
  
  print(paste0("faux négatifs med MCD ",faux_negatifsMed5[j,5,sim]))
  
  
  
  print(paste0("temps MCD ",temps[5,sim]))
  
  
  temps_ogk = system.time({
    resOGK = covOGK(data$Z,sigmamu = scaleTau2)
    invSigmaOGK = solve(resOGK$cov)
    for(s in (1:n))
    {
      if (t(data$Z[s,] - resOGK$center)%*%invSigmaOGK%*%(data$Z[s,] - resOGK$center) > qchisq(.95,df = d)){outliersLabelsMed5[s,j,6,sim] = 1}
    }
    
  }
  )
  
  temps[6,sim] = temps_ogk[3]
  
  print(paste0("temps OGK ",temps[6,sim]))
  
  
  
  
  erreursSigmaMed5[n,j,6,sim] =norm(resOGK$cov - Sigma0,"F")
  
  print(paste0("Erreur OGK med ",erreursSigmaMed5[n,j,6,sim]))
  
  
  
  t = table(data$labelsVrais,outliersLabelsMed5[,j,6,sim] )
  if (r != 0) {faux_positifsMed5[j,6,sim] =  t[1,2]
  faux_negatifsMed5[j,6,sim] = t[2,1]}
  if(r == 0){faux_positifsMed5[j,6,sim] =  t[1,2]}
  
  print(paste0("faux positifs med OGK ",faux_positifsMed5[j,6,sim]))
  
  print(paste0("faux négatifs med OGK ",faux_negatifsMed5[j,6,sim]))
  
  
  
  setwd(resDir)
  #save(Sigma_naive, Sigma_online, Sigma_str,temps_naif,temps_online,temps_streaming,file=fitFile)
  
  
}

rho1_values <- seq(-0.9, 0.9, by = 0.05)



r = 30



erreursSigmaMed5 = array(0,dim = c(n,length(rho1__values),6,simNb))
erreursInvSigmaMed5 = array(0,dim = c(n,length(rho1__values),6,simNb))
outliersLabelsMed5 = array(0,dim = c(n,length(rho1__values),6,simNb))
outliersLabelsOracleMed5 = array(0,dim = c(n,length(rho1__values),simNb))
labelsVraisMed5 = array(0,dim = c(n,length(rho1__values)))
faux_positifsMed5= array(0,dim = c(length(rho1__values),6,simNb))
faux_negatifsMed5 = array(0,dim = c(length(rho1__values),6,simNb))
faux_positifsOracleMed5 = array(0,dim = c(length(rho1__values),simNb))
faux_negatifsOracleMed5 = array(0,dim = c(length(rho1__values),simNb))
temps = array(0,dim=c(6,simNb))


outlier_sets <- vector("list", length(rList))

id_pool <- sample(1:n)

r_max <- max(rList)

outlier_sets = list()

for (m in seq_along(rList[1:9])) {
  n_active = floor(rList[m]/100*n)
  outlier_sets[[m]] = id_pool[1:n_active]
}


m = 7

for(j in seq_along(rho1_values)){
  
  
  k = 0
  l = 1
  rho1 = rho1_values[j]
  
  
  dataFile <- paste0('SimData-d', d, '-n', n, '-k', k, '-l', l, '-rho', rho1,'-r',r , '-sim', sim,".RData")
  
  print(dataFile)
  
  setwd(simDir)
  contParam = ParmsF1(m1 = m1,k1 = k,l1 = l,rho1 = rho1)
  
  if(r == 0){
    
    data = genererEchantillon_new(n,d,mu1 = mu0,mu2 = contParam$mu1,Sigma1 = Sigma0,Sigma2 = contParam$Sigma1,r)
    #save(dataFile)
    
    #save(dataFile)
  }
  if(r != 0){
    
    data = genererEchantillon_new(n,d,mu1 = mu0,mu2 = contParam$mu1,Sigma1 = Sigma0,Sigma2 = contParam$Sigma1,r, id_outliers =  outlier_sets[[m]])
  }
  #else{load(dataFile)}
  labelsVraisMed5[,j] = data$labelsVrais
  
  
  
  fitFile <- paste0('FitParms-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
  temps_naif = system.time(
    
    {resNaif = SampleCovOnline(data$Z)}
  )
  fitNaif = resNaif
  erreursSigmaMed5[n,j,1,sim] = norm(resNaif$Sigma -Sigma0,"F")
  outliersLabelsMed5[,j,1,sim] = resNaif$outliers_labels
  # 
  # for(s in (1:n)){
  #   #   erreursSigmaMed3[s,m,1,sim] = norm(resNaif$SigmaIter[s,,] - Sigma0,"F")
  #   if(t(data$Z[s,])%*% solve(Sigma0)%*% data$Z[s,] > qchisq(.95,df = d)) {
  #     outliersLabelsOracleMed3[s,m,sim] = 1
  #   }
  
  
  print(paste0("Erreur naive med ",erreursSigmaMed5[n,m,1,sim]))
  
  temps[1,sim] = temps_naif[3]
  
  
  t = table(data$labelsVrais,resNaif$outliers_labels)
  if (r != 0) {faux_positifsMed5[j,1,sim] =  t[1,2]
  faux_negatifsMed5[j,1,sim] = t[2,1]}
  if(r == 0){faux_positifsMed5[j,1,sim] =  t[1,2]}
  # 
  print(paste0("faux positifs med naive ",faux_positifsMed5[j,1,sim]))
  
  print(paste0("faux négatifs med naive ",faux_negatifsMed5[j,1,sim]))
  
  for(s in (1:n)){
    if(t(data$Z[s,])%*% solve(Sigma0)%*% data$Z[s,] > qchisq(.95,df = d)) {
      outliersLabelsOracleMed5[s,j,sim] = 1
    }
    
  }
  
  t = table(data$labelsVrais,outliersLabelsOracleMed5[,j,sim])
  if (r != 0) {faux_positifsOracleMed5[j,sim] =  t[1,2]
  faux_negatifsOracleMed5[m,sim] = t[2,1]}
  if(r == 0){faux_positifsOracleMed5[j,sim] =  t[1,2]}
  
  
  print(paste0("faux positifs Oracle ",faux_positifsOracleMed5[j,sim]))
  
  print(paste0("faux négatifs Oracle ",faux_negatifsOracleMed5[j,sim]))
  
  
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
  erreursSigmaMed5[n,j,2,sim] = norm(resUsOnline$variance - Sigma0,"F")
  if(d ==10){outliersLabelsMed5[,j,2,sim] = resUsOnline$outlier_labels}
  if(d ==100){outliersLabelsMed5[,j,2,sim] = test_outliers(distances = resUsOnline$distances,cutoff = 1.29*qchisq(.95,df = 100))}
  
  print(paste0("Erreur us online med ",erreursSigmaMed5[n,m,2,sim]))
  
  #t = table(data$labelsVrais,resUsOnline$outlier_labels)
  t = table(data$labelsVrais,outliersLabelsMed5[,j,2,sim])
  
  if (nrow(t) >= 1 && ncol(t) >= 2){
    if (r != 0) {faux_positifsMed5[j,2,sim] =  t[1,2]
    faux_negatifsMed5[j,2,sim] = t[2,1]}
    if(r == 0){faux_positifsMed5[j,2,sim] = t[1,2]}
  }
  
  print(paste0("faux positifs med us online ",faux_positifsMed5[j,2,sim]))
  
  print(paste0("faux négatifs med us online ",faux_negatifsMed5[j,2,sim]))
  
  
  
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
  erreursSigmaMed5[n,j,3,sim] =norm(resUsStreaming$variance - Sigma0,"F")
  if(d == 10){outliersLabelsMed5[,j,3,sim]= resUsStreaming$outlier_labels}
  if(d == 100){
    outliersLabelsMed5[,j,3,sim] = test_outliers(distances = resUsStreaming$distances,cutoff = 1.38*qchisq(.95,df = 100))
  }
  print(paste0("Erreur us streaming med ",erreursSigmaMed5[n,j,3,sim]))
  
  #t = table(data$labelsVrais,resUsStreaming$outlier_labels)
  t = table(data$labelsVrais,outliersLabelsMed5[,j,3,sim])
  
  if (nrow(t) >= 1 && ncol(t) >= 2){
    if(r == 0){faux_positifsMed5[j,3,sim] =  t[1,2]}
    if (r != 0) {faux_positifsMed5[j,3,sim] =  t[1,2]
    faux_negatifsMed5[j,3,sim] = t[2,1]}
  }
  
  print(paste0("faux positifs med us streaming ",faux_positifsMed5[j,3,sim]))
  
  print(paste0("faux négatifs med us streaming ",faux_negatifsMed5[j,3,sim]))
  Sigma_naive = fitNaif$SigmaIter
  Sigma_online = fitUsOnline$Sigma
  Sigma_str = fitUSStreaming$Sigma
  
  print(paste0("Erreur mcd med  ",erreursSigmaMed5[n,j,5,sim]))
  
  
  temps_offline = system.time({
    # resOffline = OfflineOutlierD etection(data$Z)
    resOffline = offlineRobustVariance(data$Z,computeOutliers = TRUE)
  }
  )
  
  
  erreursSigmaMed5[n,j,4,sim] =norm(resOffline$variance - Sigma0,"F")
  
  print(paste0("Erreur us offline med ",erreursSigmaMed5[n,j,4,sim]))
  
  
  #outliersLabelsMed5[,m,4,sim]= resOffline$outlier_labels 
  
  # fitUsOnline = resUsOnline
  temps[4,sim] = temps_offline[3]
  
  
  print(paste0("temps offline ",temps[4,sim]))
  
  
  t = table(data$labelsVrais,resOffline$outlier_labels)
  if (r != 0) {faux_positifsMed5[j,4,sim] =  t[1,2]
  faux_negatifsMed5[j,4,sim] = t[2,1]}
  if(r == 0){faux_positifsMed5[j,4,sim] =  t[1,2]}
  
  print(paste0("faux positifs med offline ",faux_positifsMed5[j,4,sim]))
  
  print(paste0("faux négatifs med offline ",faux_negatifsMed5[j,4,sim]))
  
  
  
  temps_mcd = system.time({
    resMcd = covMcd(data$Z)
    invSigmaMCD = solve(resMcd$cov) 
    for(s in (1:n))
    {
      if (t(data$Z[s,] - resMcd$center)%*%invSigmaMCD%*%(data$Z[s,] - resMcd$center) > qchisq(.95,df = d)){outliersLabelsMed5[s,j,5,sim] = 1}
    }
  }
  )
  
  erreursSigmaMed5[n,j,5,sim] =norm(resMcd$cov - Sigma0,"F")
  
  print(paste0("Erreur mcd med ",erreursSigmaMed5[n,j,5,sim]))
  
  
  
  temps[5,sim] = temps_mcd[3]
  
  t = table(data$labelsVrais,outliersLabelsMed5[,j,5,sim] )
  if (r != 0) {faux_positifsMed5[j,5,sim] =  t[1,2]
  faux_negatifsMed5[j,5,sim] = t[2,1]}
  if(r == 0){faux_positifsMed5[j,5,sim] =  t[1,2]}
  
  print(paste0("faux positifs med MCD ",faux_positifsMed5[j,5,sim]))
  
  print(paste0("faux négatifs med MCD ",faux_negatifsMed5[j,5,sim]))
  
  
  
  print(paste0("temps MCD ",temps[5,sim]))
  
  
  temps_ogk = system.time({
    resOGK = covOGK(data$Z,sigmamu = scaleTau2)
    invSigmaOGK = solve(resOGK$cov)
    for(s in (1:n))
    {
      if (t(data$Z[s,] - resOGK$center)%*%invSigmaOGK%*%(data$Z[s,] - resOGK$center) > qchisq(.95,df = d)){outliersLabelsMed5[s,j,6,sim] = 1}
    }
    
  }
  )
  
  temps[6,sim] = temps_ogk[3]
  
  print(paste0("temps OGK ",temps[6,sim]))
  
  
  
  
  erreursSigmaMed5[n,j,6,sim] =norm(resOGK$cov - Sigma0,"F")
  
  print(paste0("Erreur OGK med ",erreursSigmaMed5[n,j,6,sim]))
  
  
  
  t = table(data$labelsVrais,outliersLabelsMed5[,j,6,sim] )
  if (r != 0) {faux_positifsMed5[j,6,sim] =  t[1,2]
  faux_negatifsMed5[j,6,sim] = t[2,1]}
  if(r == 0){faux_positifsMed5[j,6,sim] =  t[1,2]}
  
  print(paste0("faux positifs med OGK ",faux_positifsMed5[j,6,sim]))
  
  print(paste0("faux négatifs med OGK ",faux_negatifsMed5[j,6,sim]))
  
  
  
  setwd(resDir)
  #save(Sigma_naive, Sigma_online, Sigma_str,temps_naif,temps_online,temps_streaming,file=fitFile)
  
  
}
