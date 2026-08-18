scenarios = c(scenarios_1_param,scenarios_2_param)
#############################Final Frobenius norm error, false positives, false negatives V2##############################

# Pour la norme de Frobeniusmethodes = c("SampleNaiveQuantonlinecorr","OnlineUsQuantonlinecorr","StreamingUsonlineQuantcorr","OfflinewithQuantcorr","OGK","MCD")
methodes_add  = c("SampleNaivewithoutonlinequantilecorr","OnlineUswithoutQuantonlinecorr","StreamingUswithoutQuantonlinecorr","OfflineUswithoutQuantcorr","OracleRD","OracleQC","SampleRaw","OnlRaw","StrmRaw","OfflRaw","OGKRD","OGKQC","MCDRD","MCDQC")
methode_oracle = c("Oracle")
methodes_frob <- c(
  "SampleNaiveQuantonlinecorr",
  "OnlineUsQuantonlinecorr",
  "StreamingUsonlineQuantcorr",
  "OfflinewithQuantcorr",
  "OGK",
  "MCD"
)


# Méthodes QC
methodes_qc <- c(
  "SampleNaiveQuantonlinecorr",
  "OnlineUsQuantonlinecorr",
  "StreamingUsonlineQuantcorr",
  "OfflinewithQuantcorr",
  "OGKQC",
  "MCDQC",
  "OracleQC"
)

# Méthodes RD
methodes_rd <- c(
  "SampleNaivewithoutonlinequantilecorr",
  "OnlineUswithoutQuantonlinecorr",
  "StreamingUswithoutQuantonlinecorr",
  "OfflineUswithoutQuantcorr",
  "OGKRD",
  "MCDRD",
  "OracleRD"
)

# Méthodes Raw
methodes_raw <- c(
  "SampleRaw",
  "OnlRaw",
  "StrmRaw",
  "OfflRaw",
  "OGK",
  "MCD",
  "Oracle"
)


# Ordre global des méthodes
all_methodes <- c(
  methodes_qc,
  methodes_rd,
  methodes_raw
)


## Méthodes à afficher
idxFrob <- which(all_methodes %in% methodes_frob)
idxAll  <- seq_along(all_methodes)
idxAUC  <- which(all_methodes %in% methodes_frob)



#############################
# Couleurs
#############################

cols <- c(
  # QC
  "darkgreen", # SampleNaiveQuantonlinecorr
  "blue",      # OnlineUsQuantonlinecorr
  "red",       # StreamingUsonlineQuantcorr
  "orange",    # OfflinewithQuantcorr
  "brown",     # OGKQC
  "black",     # MCDQC
  "purple4",   # OracleQC
  
  # RD
  "darkgreen", # SampleNaivewithoutonlinequantilecorr
  "blue",      # OnlineUswithoutQuantonlinecorr
  "red",       # StreamingUswithoutQuantonlinecorr
  "orange",    # OfflineUswithoutQuantcorr
  "brown",     # OGKRD
  "black",     # MCDRD
  "purple4",   # OracleRD
  
  # Raw
  "darkgreen", # SampleRaw
  "blue",      # OnlRaw
  "red",       # StrmRaw
  "orange",    # OfflRaw
  "brown",     # OGK
  "black",     # MCD
  "purple4"    # Oracle
)



#############################
# Symboles
#############################

pchs <- rep(NA, length(all_methodes))


# Étoiles : correction quantile (QC)
pchs[all_methodes %in% methodes_qc] <- 8


# Carrés : rescale distances (RD)
pchs[all_methodes %in% methodes_rd] <- 15


# Triangles : Raw
pchs[all_methodes %in% methodes_raw] <- 17



# Transparence
cols_alpha <- adjustcolor(cols, alpha.f = 0.65)
for (sc in scenarios){
  k = sc$k
  l = sc$l
  rho1 = sc$rho1
  
  erreursSigmaPlot = array(0,dim = c(length(rList[1:13]),length(all_methodes)))
  faux_positifsPlot= array(0,dim = c(length(rList[1:13]),length(all_methodes)))
  faux_negatifsPlot= array(0,dim = c(length(rList[1:13]),length(all_methodes)))
  ariPlot = array(0,dim = c(length(rList[1:13]),length(all_methodes)))
  aucPlot = array(0,dim = c(length(rList[1:13]),length(all_methodes)))
  propHorsDiagPlot = array(0,dim = c(length(rList[1:13]),length(all_methodes)))
  
    
  for (m in seq_along(rList[1:13])){
    
    r = rList[m]
    
    
    for(j in seq_along(all_methodes)){
      
      
      methode = all_methodes[j]
      
      setwd(criteres)
      
      critFile <- paste0(
        'Crit-',methode, "-d",d,
        '-n', n,
        '-k', k,
        '-l', l,
        '-rho', rho1,
        '-r', r,
        '-mean',
        ".RData"
      )
      load(critFile)
      if(methode %in% methodes){
        
      erreursSigmaPlot[m,j] = crit_mean$erreurFrob
      #faux_negatifsPlot[m,j] = crit_mean$FN
      faux_positifsPlot[m,j] = crit_mean$FP
      #ariPlot[m,j] = crit_mean$ARI
      #propHorsDiagPlot[m,j] = crit_mean$prop_hors_diag
      
      if( r!= 0){
      aucPlot[m,j] = crit_mean$AUC
      faux_negatifsPlot[m,j] = crit_mean$FN
      faux_positifsPlot[m,j] = crit_mean$FP
      ariPlot[m,j] = crit_mean$ARI
      propHorsDiagPlot[m,j] = crit_mean$prop_hors_diag
      
      }
      
      }
      
      
      if (methode %in% c("Oracle")){
                faux_positifsPlot[m,j] = crit_mean$FP
        if(r != 0){
        
        aucPlot[m,j] = crit_mean$AUC
        faux_negatifsPlot[m,j] = crit_mean$FN
        propHorsDiagPlot[m,j] = crit_mean$prop_hors_diag
        ariPlot[m,j] = crit_mean$ARI
        
        }
      }
      
      if (methode %in% methodes_add){
        if(r == 0){
          faux_positifsPlot[m,j] = crit_mean$FP
          
        }
        if(r != 0){
        faux_negatifsPlot[m,j] = crit_mean$FN
        faux_positifsPlot[m,j] = crit_mean$FP
        ariPlot[m,j] = crit_mean$ARI
        propHorsDiagPlot[m,j] = crit_mean$prop_hors_diag}
      }
      
    } }
  
  
  
  
  
  #########################Graphiques##########################
  
  
  setwd(figures)
  
  file <- paste0("scen-k", k, "-l", l, "-rho1", rho1,".pdf")
  
  pdf(file, width = 25, height = 4)
  par(mfrow = c(1,5), mar = c(4,4,2,1))
  ############################Erreur norme de Frobenius#################
  plot(rList[1:13],
       erreursSigmaPlot[,idxFrob[1]],
       type="l",
       log="y",
       lwd=2,
       col=cols_alpha[idxFrob[1]],
       ylim=c(1e-1,1e2),
       xaxt="n", yaxt="n",
       xlab="", ylab="")
  
  for(i in idxFrob[-1]){
    lines(rList[1:13],
          erreursSigmaPlot[,i],
          lwd=2,
          col=cols_alpha[i])
  }
  
  axis(1,
       at=rList[1:13],
       las=1,
       cex.axis=1.8)
  
  axis(2,
       at=10^seq(-1,2),
       labels=parse(text=paste0("10^",-1:2)),
       las=1,
       cex.axis=1.8)
  
  box()
  
    ###########################Faux négatifs#######################
  
  plot(rList[2:13],
       faux_negatifsPlot[2:13,1]/((rList[2:13]/100)*n)*100,
       type="l",
       lwd=2,
       col=cols_alpha[1],
       ylim=c(0,100),
       xaxt="n",
       yaxt="n",
       xlab="",
       ylab="")
  
  points(rList[2:13],
         faux_negatifsPlot[2:13,1]/((rList[2:13]/100)*n)*100,
         pch=pchs[1],
         col=cols_alpha[1],
         cex=1.1)
  
  for(i in idxAll[-1]){
    
    lines(rList[2:13],
          faux_negatifsPlot[2:13,i]/((rList[2:13]/100)*n)*100,
          lwd=2,
          col=cols_alpha[i])
    
    points(rList[2:13],
           faux_negatifsPlot[2:13,i]/((rList[2:13]/100)*n)*100,
           pch=pchs[i],
           col=cols_alpha[i],
           cex=1.1)
  }
  
  axis(1,at=rList[-1],las=1,cex.axis=1.8)
  axis(2,las=1,cex.axis=1.8)
  box()  

  
  #######################Faux positifs###########################
  plot(rList[1:13],
       faux_positifsPlot[1:13,1]/((1-rList[1:13]/100)*n)*100,
       type="l",
       lwd=2,
       col=cols_alpha[1],
       ylim=c(0,20),
       xaxt="n",
       yaxt="n",
       xlab="",
       ylab="")
  
  points(rList[1:13],
         faux_positifsPlot[1:13,1]/((1-rList[1:13]/100)*n)*100,
         pch=pchs[1],
         col=cols_alpha[1],
         cex=1.1)
  
  for(i in idxAll[-1]){
    
    lines(rList[1:13],
          faux_positifsPlot[1:13,i]/((1-rList[1:13]/100)*n)*100,
          lwd=2,
          col=cols_alpha[i])
    
    points(rList[1:13],
           faux_positifsPlot[1:13,i]/((1-rList[1:13]/100)*n)*100,
           pch=pchs[i],
           col=cols_alpha[i],
           cex=1.1)
  }
  
  axis(1,at=rList[1:13],las=1,cex.axis=1.8)
  axis(2,las=1,cex.axis=1.8)
  box()
  
  ##############################AUC############################
  plot(rList[2:13],
       aucPlot[2:13,idxAUC[1]],
       type="l",
       lwd=2,
       col=cols_alpha[idxAUC[1]],
       ylim=c(0,1),
       xaxt="n",
       yaxt="n",
       xlab="",
       ylab="")
  
  
  for(i in idxAUC[-1]){
    
    lines(rList[2:13],
          aucPlot[2:13,i],
          lwd=2,
          col=cols_alpha[i])
  }
  
  
  axis(1,
       at=rList[2:13],
       las=1,
       cex.axis=1.8)
  
  axis(2,
       at=seq(0,1,0.1),
       las=1,
       cex.axis=1.8)
  
  box()
  #############################ARI#####################################
  # 
  # plot(rList[1:13],
  #      ariPlot[1:13,1],
  #      type="l",
  #      lwd=2,
  #      col=cols_alpha[1],
  #      ylim=c(-1,1),
  #      xaxt="n",
  #      yaxt="n",
  #      xlab="",
  #      ylab="")
  # 
  # points(rList[1:13],
  #        ariPlot[1:13,1],
  #        pch=pchs[1],
  #        col=cols_alpha[1],
  #        cex=1.1)
  # 
  # for(i in idxAll[-1]){
  #   
  #   lines(rList[1:13],
  #         ariPlot[1:13,i],
  #         lwd=2,
  #         col=cols_alpha[i])
  #   
  #   points(rList[1:13],
  #          ariPlot[1:13,i],
  #          pch=pchs[i],
  #          col=cols_alpha[i],
  #          cex=1.1)
  # }
  # 
  # axis(1,at=rList[1:13],las=1,cex.axis=1.8)
  # axis(2,at=seq(-1,1,0.1),las=1,cex.axis=1.8)
  # box()    
  ##########################Accuracy####################
  
  
  plot(rList[2:13],
       1 - propHorsDiagPlot[2:13,1],
       type = "l",
       lwd = 2,
       col = cols_alpha[1],
       ylim = c(0,1),
       xaxt = "n",
       yaxt = "n",
       xlab = "",
       ylab = "")
  
  points(rList[2:13],
         1 - propHorsDiagPlot[2:13,1],
         pch = pchs[1],
         col = cols_alpha[1],
         cex = 1.1)
  
  for(i in idxAll[-1]){
    
    lines(rList[2:13],
          1 - propHorsDiagPlot[2:13,i],
          lwd = 2,
          col = cols_alpha[i])
    
    points(rList[2:13],
           1 - propHorsDiagPlot[2:13,i],
           pch = pchs[i],
           col = cols_alpha[i],
           cex = 1.1)
  }
  
  axis(1,
       at = rList[2:13],
       las = 1,
       cex.axis = 1.8)
  
  axis(2,
       las = 1,
       cex.axis = 1.8)
  
  box()
  
  dev.off()
  
  
  
  }
    
  
####################################Temps de calculs#########

temps_calcul = array(0,dim = c(length(methodes_frob),simNb))

####Extraction

k = 0; l = 0.01;rho1 = 0.3; r= 2


for(sim in 1:simNb){
for (m in seq_along(methodes_frob)){
  
  methode = methodes_frob[m]
  
  setwd(resAlgo)
  fitFile <- paste0(
    "Fit-", methode,
    "-d", d,
    "-n", n,
    "-k", k,
    "-l", l,
    "-rho", rho1,
    "-r", r,
    "-sim", sim,
    ".RData"
  )
  
  load(fitFile)
  
  temps_calcul[m,sim] = resultats$temps[3]
  
  
  
}}


setwd(figures)


noms <- c(
  "Sple",
  "Onl",
  "Strm",
  "Offl",
  "OGK",
  "MCD"
)

cols <- c(
  "darkgreen", # Sample naive
  "blue",      # Online
  "red",       # Streaming
  "orange",    # Offline
  "brown",     # OGK
  "black"      # MCD
)


file <- paste0("temps_calcul-n", n, "-d", d,".pdf")

pdf(file, width = 25, height = 7)


boxplot(
  t(temps_calcul),
  names = noms,
  col = cols,
  las = 2,
  log = "y",
  ylab = "",
  xlab = "",
  main = "",
  yaxt = "n",
  ylim = c(10^-2, 10^2),
  cex.axis = 2.5
)

axis(
  2,
  at = 10^seq(-2, 2),
  labels = parse(text = paste0("10^", -2:2)),
  las = 1,
  cex.axis = 1.9
)


dev.off()

#######################Trajectoires##########################

methodes_online_quantile = c(
  "SampleNaiveQuantonlinecorr",
  "OnlineUsQuantonlinecorr",
  "StreamingUsonlineQuantcorr",
  "OracleQC"
)

methodes_online_rescale = c(
  "SampleNaivewithoutonlinequantilecorr",
  "OnlineUswithoutQuantonlinecorr",
  "StreamingUswithoutQuantonlinecorr",
  "OracleRD"
)

methodes_online_raw = c(
  "SampleRaw",
  "OnlRaw",
  "StrmRaw",
  "Oracle"
)

methodes_online = c(
  methodes_online_quantile,
  methodes_online_rescale,
  methodes_online_raw
)


# Symboles associés aux familles
pch_methodes = rep(NA,length(methodes_online))

# Correction quantile : étoile
pch_methodes[methodes_online %in% methodes_online_quantile] = 8

# Rescale distance : carré
pch_methodes[methodes_online %in% methodes_online_rescale] = 15

# Raw : triangle
pch_methodes[methodes_online %in% methodes_online_raw] = 17



# Symboles associés aux familles
pch_methodes = rep(NA,length(methodes_online))

pch_methodes[methodes_online %in% methodes_online_quantile] = 15 # carré
pch_methodes[methodes_online %in% methodes_online_rescale]  = 8  # étoile
pch_methodes[methodes_online %in% methodes_online_raw]      = 17 # triangle



for(sc in scenarios){
  
  k = sc$k
  l = sc$l
  rho1 = sc$rho1
  
  
  for(j in seq_along(rList[1:13])){
    
    r = rList[j]
    
    outlabTraj = array(0, dim = c(n,length(methodes_online)))
    
    
    setwd(SimDir)
    
    dataFile = paste0(
      'SimData-d', d,
      '-n', n,
      '-k', k,
      '-l', l,
      '-rho', rho1,
      '-r', r,
      "-sim", sim,
      ".RData"
    )
    
    print(dataFile)
    
    load(dataFile)
    
    
    Z_clean = data$Z[data$labelsVrais == 0,]
    
    if(r != 0){
      Z_cont = data$Z[data$labelsVrais == 1,]
    }
    
    labels = data$labelsVrais
    
    
    nboutliers = r/100*n
    nbinliers = (1-r/100)*n
    
    
    distinliers = rep(0,nbinliers)
    
    invSigma0 = solve(Sigma0)
    
    
    for(m in 1:nbinliers){
      
      distinliers[m] =
        t(Z_clean[m,]-mu0)%*%
        invSigma0%*%
        (Z_clean[m,]-mu0)
      
    }
    
    
    if(r != 0){
      
      distoutliers = rep(0,nboutliers)
      
      for(m in 1:nboutliers){
        
        distoutliers[m] =
          t(Z_cont[m,]-mu0)%*%
          invSigma0%*%
          (Z_cont[m,]-mu0)
        
      }
    }
    
    
    # ==========================
    # Chargement des résultats
    # ==========================
    
    for(s in seq_along(methodes_online)){
      
      methode = methodes_online[s]
      
      setwd(resAlgo)
      
      fitFile = paste0(
        'Fit-',methode,
        '-d', d,
        '-n', n,
        '-k', k,
        '-l', l,
        '-rho', rho1,
        '-r', r,
        '-sim', sim,
        ".RData"
      )
      
      load(fitFile)
      
      outlabTraj[,s] = resultats$outliers_labels
      
      if(methode == "StreamingUsonlineQuantcorr"){
        diststrmqc = resultats$distances
      }
      
    }
    
    
    
    # ==========================
    # Taux
    # ==========================
    
    rates_samplecov_wcorr =
      compute_rates(outlabTraj[,1],labels)
    
    rates_online_corr =
      compute_rates(outlabTraj[,2],labels)
    
    rates_strm_corr =
      compute_rates(outlabTraj[,3],labels)
    
    rates_oracle_qc =
      compute_rates(outlabTraj[,4],labels)
    
    
    rates_samplecov_wocorr =
      compute_rates(outlabTraj[,5],labels)
    
    rates_online_wocorr =
      compute_rates(outlabTraj[,6],labels)
    
    rates_strm_wocorr =
      compute_rates(outlabTraj[,7],labels)
    
    rates_oracle_rd =
      compute_rates(outlabTraj[,8],labels)
    
    
    rates_samplecov_raw =
      compute_rates(outlabTraj[,9],labels)
    
    rates_online_raw =
      compute_rates(outlabTraj[,10],labels)
    
    rates_strm_raw =
      compute_rates(outlabTraj[,11],labels)
    
    rates_oracle =
      compute_rates(outlabTraj[,12],labels)
    
    
    
    # ==========================
    # Figure
    # ==========================
    
    setwd(figures)
    
    nom_fichier = paste0(
      "trajectories_k-",k,
      "-l",l,
      "-rho1",rho1,
      "-r",r,
      "-sim",sim,
      ".pdf"
    )
    
    
    pdf(nom_fichier,width=14,height=10)
    
    par(mfrow=c(2,2),
        mar=c(4,4,2,1))
    
    
    x_vals = 1:length(rates_strm_corr$FN_rate)
    
    # =====================================================
    # BOXplot distances
    # =====================================================
    
    if(r != 0){
      
      boxplot(distinliers,
              distoutliers,
              
              col=c("lightblue","darkred","red"),
              names=c("Inliers","Outliers"),
              ylim=c(0,30),
              main="Mahalanobis distances",
              ylab="Distance",
              cex.axis=1.5,
              cex.lab=1.5,
              cex.main=1.5)
      
    }
    
    
    # =====================================================
    # FALSE NEGATIVE RATE
    # =====================================================
    
    plot(x_vals,
         rates_strm_corr$FN_rate*100,
         type="l",
         lwd=2,
         col="red",
         ylim=c(0,100),
         xlab="",
         ylab="",
         main="False Negative Rate",
         xaxt="n",
         yaxt="n")
    
    
    # fonction pour ajouter les symboles espacés
    add_points <- function(x,y,pch,col){
      id <- seq(1,length(x),by=1000)
      points(x[id],
             y[id],
             pch=pch,
             col=col,
             cex=1.3)
    }
    
    
    
    # Quantile
    
    lines(x_vals,
          rates_samplecov_wcorr$FN_rate*100,
          col="darkgreen",
          lwd=2)
    add_points(x_vals,
               rates_samplecov_wcorr$FN_rate*100,
               pch_methodes[1],
               "darkgreen")
    
    
    lines(x_vals,
          rates_online_corr$FN_rate*100,
          col="blue",
          lwd=2)
    add_points(x_vals,
               rates_online_corr$FN_rate*100,
               pch_methodes[2],
               "blue")
    
    
    lines(x_vals,
          rates_strm_corr$FN_rate*100,
          col="red",
          lwd=2)
    add_points(x_vals,
               rates_strm_corr$FN_rate*100,
               pch_methodes[3],
               "red")
    
    
    lines(x_vals,
          rates_oracle_qc$FN_rate*100,
          col="purple",
          lwd=2)
    add_points(x_vals,
               rates_oracle_qc$FN_rate*100,
               pch_methodes[4],
               "purple")
    
    
    
    # Rescale
    
    lines(x_vals,
          rates_samplecov_wocorr$FN_rate*100,
          col="darkgreen",
          lwd=2)
    add_points(x_vals,
               rates_samplecov_wocorr$FN_rate*100,
               pch_methodes[5],
               "darkgreen")
    
    
    lines(x_vals,
          rates_online_wocorr$FN_rate*100,
          col="blue",
          lwd=2)
    add_points(x_vals,
               rates_online_wocorr$FN_rate*100,
               pch_methodes[6],
               "blue")
    
    
    lines(x_vals,
          rates_strm_wocorr$FN_rate*100,
          col="red",
          lwd=2)
    add_points(x_vals,
               rates_strm_wocorr$FN_rate*100,
               pch_methodes[7],
               "red")
    
    
    lines(x_vals,
          rates_oracle_rd$FN_rate*100,
          col="purple",
          lwd=2)
    add_points(x_vals,
               rates_oracle_rd$FN_rate*100,
               pch_methodes[8],
               "purple")
    
    
    
    # Raw
    
    lines(x_vals,
          rates_samplecov_raw$FN_rate*100,
          col="darkgreen",
          lwd=2)
    add_points(x_vals,
               rates_samplecov_raw$FN_rate*100,
               pch_methodes[9],
               "darkgreen")
    
    
    lines(x_vals,
          rates_online_raw$FN_rate*100,
          col="blue",
          lwd=2)
    add_points(x_vals,
               rates_online_raw$FN_rate*100,
               pch_methodes[10],
               "blue")
    
    
    lines(x_vals,
          rates_strm_raw$FN_rate*100,
          col="red",
          lwd=2)
    add_points(x_vals,
               rates_strm_raw$FN_rate*100,
               pch_methodes[11],
               "red")
    
    
    lines(x_vals,
          rates_oracle$FN_rate*100,
          col="purple",
          lwd=2)
    add_points(x_vals,
               rates_oracle$FN_rate*100,
               pch_methodes[12],
               "purple")
    
    
    
    axis(2,
         las=1,
         cex.axis=1.8)
    
    axis(1,
         at=seq(1000,max(x_vals),by=1000),
         las=1,
         cex.axis=1.8)
    
    box()
    # =====================================================
    # FALSE POSITIVE RATE
    # =====================================================
    
    plot(x_vals,
         rates_strm_corr$FP_rate*100,
         type="l",
         lwd=2,
         col="red",
         ylim=c(0,20),
         xlab="",
         ylab="",
         main="False Positive Rate",
         xaxt="n",
         yaxt="n")
    
    
    # fonction pour ajouter les points espacés
    add_points <- function(x,y,pch,col){
      id <- seq(1,length(x),by=1000)
      points(x[id],
             y[id],
             pch=pch,
             col=col,
             cex=1.3)
    }
    
    
    # Quantile
    lines(x_vals,rates_samplecov_wcorr$FP_rate*100,
          col="darkgreen",lwd=2)
    add_points(x_vals,rates_samplecov_wcorr$FP_rate*100,
               pch_methodes[1],"darkgreen")
    
    
    lines(x_vals,rates_online_corr$FP_rate*100,
          col="blue",lwd=2)
    add_points(x_vals,rates_online_corr$FP_rate*100,
               pch_methodes[2],"blue")
    
    
    lines(x_vals,rates_strm_corr$FP_rate*100,
          col="red",lwd=2)
    add_points(x_vals,rates_strm_corr$FP_rate*100,
               pch_methodes[3],"red")
    
    
    lines(x_vals,rates_oracle_qc$FP_rate*100,
          col="purple",lwd=2)
    add_points(x_vals,rates_oracle_qc$FP_rate*100,
               pch_methodes[4],"purple")
    
    
    
    # Rescale
    lines(x_vals,rates_samplecov_wocorr$FP_rate*100,
          col="darkgreen",lwd=2)
    add_points(x_vals,rates_samplecov_wocorr$FP_rate*100,
               pch_methodes[5],"darkgreen")
    
    
    lines(x_vals,rates_online_wocorr$FP_rate*100,
          col="blue",lwd=2)
    add_points(x_vals,rates_online_wocorr$FP_rate*100,
               pch_methodes[6],"blue")
    
    
    lines(x_vals,rates_strm_wocorr$FP_rate*100,
          col="red",lwd=2)
    add_points(x_vals,rates_strm_wocorr$FP_rate*100,
               pch_methodes[7],"red")
    
    
    lines(x_vals,rates_oracle_rd$FP_rate*100,
          col="purple",lwd=2)
    add_points(x_vals,rates_oracle_rd$FP_rate*100,
               pch_methodes[8],"purple")
    
    
    
    # Raw
    lines(x_vals,rates_samplecov_raw$FP_rate*100,
          col="darkgreen",lwd=2)
    add_points(x_vals,rates_samplecov_raw$FP_rate*100,
               pch_methodes[9],"darkgreen")
    
    
    lines(x_vals,rates_online_raw$FP_rate*100,
          col="blue",lwd=2)
    add_points(x_vals,rates_online_raw$FP_rate*100,
               pch_methodes[10],"blue")
    
    
    lines(x_vals,rates_strm_raw$FP_rate*100,
          col="red",lwd=2)
    add_points(x_vals,rates_strm_raw$FP_rate*100,
               pch_methodes[11],"red")
    
    
    lines(x_vals,rates_oracle$FP_rate*100,
          col="purple",lwd=2)
    add_points(x_vals,rates_oracle$FP_rate*100,
               pch_methodes[12],"purple")
    
    
    
    axis(2,
         las=1,
         cex.axis=1.8)
    
    axis(1,
         at=seq(1000,max(x_vals),by=1000),
         las=1,
         cex.axis=1.8)
    
    box()
    dev.off()
    
  }
}