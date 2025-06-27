
calcule_tout = function(cutoff = qchisq(.95,df = d),contamin = "moyenne",nbrows  = n,nb_runs = nbruns,cluster = FALSE,reduction_dim = FALSE){
  
  
  methodes = c("sampleCovOnline","online","streaming")
  
  taux_contamination <- c(0, 2, 5, 10, 15, 20, 25, 30, 40)
  
  
  rmseSigmaRec = array(0, dim = c(nbrows, length(taux_contamination), length(methodes),nb_runs))
  rmseMedRec = array(0, dim = c(nbrows, length(taux_contamination), length(methodes),nb_runs))
  distancesRec = array(0, dim = c(nbrows, length(taux_contamination), length(methodes),nb_runs))
  outliersLabelsRec = array(0, dim = c(nbrows, length(taux_contamination), length(methodes),nb_runs))
  labelsVraisRec = array(0, dim = c(nbrows, length(taux_contamination)))
  temps_calcul = array(0, dim = c( length(taux_contamination), length(methodes),nb_runs))
  faux_positifsRec = array(0, dim = c(nbrows, length(taux_contamination), length(methodes),nb_runs))
  faux_negatifsRec = array(0, dim = c(nbrows, length(taux_contamination), length(methodes),nb_runs))
  taux_OutliersVraisRec = array(0, dim = c(nbrows, length(taux_contamination), length(methodes),nb_runs))
  taux_OutliersDetectesVraisRec = array(0, dim = c(nbrows, length(taux_contamination), length(methodes),nb_runs))
  taux_OutliersDetectesRec = array(0, dim = c(nbrows, length(taux_contamination), length(methodes),nb_runs))
  aucRec = array(0, dim = c(length(taux_contamination),length(methodes), nb_runs))
  
  
  
  for (k in seq_along(taux_contamination))
  {
    
    r = taux_contamination[k]
    print(paste("r = ",r))
    sigmaSq0 <- (1:d); sigmaSq0 <- sigmaSq0 / mean(sigmaSq0)
    SigmaContamin <- diag(sqrt(sigmaSq0)) %*% toeplitz(0.8^(0:(d-1))) %*% diag(sqrt(sigmaSq0))
    
    data <- genererEchantillon(n,d,mu1,mu2 = 30*rep(1/sqrt(d), d),p1 = 1- r/100,r/100,Sigma1,Sigma2 = 2*SigmaContamin,contamin,cluster)
    
    Z = data$Z
    
    labelsVraisRec[,k] = data$labelsVrais
    #lvc = labelsVraisRec[,k]
    #print("nombre outliers =" ,sum(lvc == 1))
    #compt = 1  
    for (j in 1:nbruns)
    {  
      #compt_meth = 1    
      for (l in seq_along(methodes))
      {
        m = methodes[l]
        
        if(m == "sampleCovOnline")
        {
          temps = system.time(
            {resultats = SampleCovOnline(Z)
            mu_hat = resultats$mean
            Sigma = resultats$Sigma  
            mu_hatIter = resultats$meanIter
            SigmaIter = resultats$SigmaIter
            
            
            
            distances = rep(0, nrow(Z))
            outliers_labels = rep(0,nrow(Z))
            
            for (i in (1:nrow(Z))){
              distances[i] = mahalanobis_generalizedRcpp(Z[i,],mu_hatIter[i,],eigen(SigmaIter[i,,])$vectors, eigen(SigmaIter[i,,])$values)
              S = distances[i]
              
              if (S > cutoff) {outliers_labels[i] = 1}
              
            }
            resultats$distances = distances
            resultats$outlier_labels = outliers_labels
            }
          )
              }
        
        
        if(m == "samplecovOffline")
        {
          
          resultats <- list()
          mu_hat = colMeans(Z)
          Sigma = cov(Z)  
          
          distances = rep(0, nrow(Z))
          outliers_labels = rep(0,nrow(Z))
          
          for (i in (1:nrow(Z))){
            distances[i] = mahalanobis_generalizedRcpp(Z[i,],mu_hat,eigen(Sigma)$vectors, eigen(Sigma)$values)
            S = distances[i]
            
            if (S > cutoff) {outliers_labels[i] = 1}
            
          }
          resultats$distances = distances
          resultats$outlier_labels = outliers_labels
        }
        
        if(m == "samplecovTrimmed")
        {
          temps = system.time(
            {resultats = update_mean_Sigma2Trimmed(Z,cutoff = qchisq(.95,df = d))}
          )
          mu_hat = resultats$mean
          Sigma = resultats$Sigma2  
          
          mu_hatIter = resultats$mean_iter
          SigmaIter = resultats$Sigma2_iter
          
          distances = rep(0, nrow(Z))
          outliers_labels = rep(0,nrow(Z))
          
          for (i in (1:nrow(Z))){
            
            distances[i] = mahalanobis_generalizedRcpp(Z[i,],mu_hat,eigen(Sigma)$vectors, eigen(Sigma)$values)
            S = distances[i]
            
            if (S > cutoff) {outliers_labels[i] = 1}
            
          }
          resultats$distances = distances
          resultats$outlier_labels = outliers_labels
        }  
        
        if(m == "OGK")
        {
          
          temps = system.time(
            {resultats = covOGK(Z,sigmamu = s_mad)})
          mu_hat <- resultats$center
          Sigma <- resultats$cov
          
          
          distances = rep(0, nrow(Z))
          outliers_labels = rep(0,nrow(Z))
          
          for (i in (1:nrow(Z))){
            distances[i] = mahalanobis_generalizedRcpp(Z[i,],mu_hat,eigen(Sigma)$vectors, eigen(Sigma)$values)
            S = distances[i]
            
            if (S > cutoff) {outliers_labels[i] = 1}
            
          }
          resultats$distances = distances
          resultats$outlier_labels = outliers_labels
          
        }
        
        if(m == "FASTMCD"){
          temps = system.time({resultats = covMcd(Z)}
          )
          
          mu_hat <- resultats$center
          Sigma <- resultats$cov
          
        }
        
        
        if(m == "comedianeOffline") {
          temps = system.time(
            {   resultats = covComed(Z,n.iter = 0)}
          )
          mu_hat = resultats$center
          Sigma = resultats$cov
          distances = rep(0, nrow(Z))
          outliers_labels = rep(0,nrow(Z))
          
          for (i in (1:nrow(Z))){
            distances[i] = mahalanobis_generalizedRcpp(Z[i,],mu_hat,eigen(Sigma)$vectors, eigen(Sigma)$values)
            S = distances[i]
            
            if (S > cutoff) {outliers_labels[i] = 1}
            
          }
          resultats$distances = distances
          resultats$outlier_labels = outliers_labels
          
          
          
          distances = rep(0, nrow(Z))
          outliers_labels = rep(0,nrow(Z))
          
          for (i in (1:nrow(Z))){
            distances[i] = mahalanobis_generalizedRcpp(Z[i,],mu_hat,eigen(Sigma)$vectors, eigen(Sigma)$values)
            S = distances[i]
            
            if (S > cutoff) {outliers_labels[i] = 1}
            
          }
          resultats$distances = distances
          resultats$outlier_labels = outliers_labels
          
          
          outliers_labels = resultats$outlier_labels
          distances = resultats$distances
          
          
          
        }  
        
        if(m == "comedianeOfflineShrinkage") {
          resultats = list()
          temps = system.time(
            {   resultats = covComed(Z,n.iter = 0)}
          )
          mu_hat <- shrinkage_med(Z)$muShrink
          Sigma <- shrinkage_SCCM(Z,k = 1)$SCCM_shrinked
          
          
          
          distances = rep(0, nrow(Z))
          outliers_labels = rep(0,nrow(Z))
          
          for (i in (1:nrow(Z))){
            distances[i] = mahalanobis_generalizedRcpp(Z[i,],mu_hat,eigen(Sigma)$vectors, eigen(Sigma)$values)
            S = distances[i]
            
            if (S > cutoff) {outliers_labels[i] = 1}
            
          }
          resultats$distances = distances
          resultats$outlier_labels = outliers_labels
          
          
          
          
        } 
        
        
        
        
        if(m == "offline"){
          
          temps = system.time(
            {   resultats = OfflineOutlierDetection(Z)}
          )
          mu_hat = resultats$median
          Sigma = resultats$variance
          
          outliers_labels = resultats$outlier_labels
          distances = resultats$distances
          
        }
        if (m == "online") {
          if (ncol(Z) == 100) {
            temps <- system.time(
              resultats <- StreamingOutlierDetection(Z, batch = 1, cutoff = 1.27 * qchisq(0.95, df = ncol(Z)))
              
            )
          } else if(ncol(Z) == 10){
            temps <- system.time(
              resultats <- StreamingOutlierDetection(Z, batch = 1)
            )
          }
          mu_hat = resultats$moyennem
          mu_hatIter = resultats$miter
          SigmaIter = resultats$Sigma
          Sigma = resultats$Sigma[nrow(Z),,]
          
          distances = resultats$distances
          outliers_labels = resultats$outlier_labels
          
          
        }
        
        if (m == "streaming")
        {
          if(ncol(Z)== 100){
          print(paste0("d = ",ncol(Z)))
            temps <- system.time(
            resultats <- StreamingOutlierDetection(Z, batch = sqrt(ncol(Z)),cutoff = 1.38*qchisq(0.95,df = d))
          )}
        } else if(ncol(Z) == 10) {
          temps <- system.time(
            resultats <- StreamingOutlierDetection(Z, batch = ncol(Z))
          )
          
          mu_hat = resultats$moyennem
          mu_hatIter = resultats$miter
          Sigma = resultats$Sigma[nrow(Z),,]
          SigmaIter = resultats$Sigma
          outliers_labels = resultats$outlier_labels
          distances = resultats$distances
        }    
        # }
        # tc <- table(data$labelsVrais[1:nrow(Z)], as.numeric(outliers_labels))
        # tc <- safe_access_tc(tc)
        # 
        # if ((tc["0", "0"] + tc["0", "1"]) != 0) {
        #   if (m == "offline"){fp_offline[compt] <- (tc["0", "1"] / (tc["0", "1"] + tc["0", "0"])) * 100}
        #   if (m == "online"){fp_online[compt] <- (tc["0", "1"] / (tc["0", "1"] + tc["0", "0"])) * 100}
        #   if (m == "streaming"){fp_streaming[compt] <- (tc["0", "1"] / (tc["0", "1"] + tc["0", "0"])) * 100}
        # }
        # 
        #faux_positifsRec[,compt_meth,compt,j] = cumulativeOutlierDetection(resultats,data,r,"Shifted Gaussian contamination scenario")
        
        #faux_negatifsRec[,compt_meth,compt,j] = cumulativeOutlierDetection(resultats,data,r,"Shifted Gaussian contamination scenario")$taux_outliers_detectes_vraiss
        
        if( !m %in% c("sampleCovOffline","comedianeOffline","comedianeOfflineShrinkage","OGK","FASTMCD","offline")){
          print(paste("méthode ",m," "))
          print(paste("rmseSigma ",rmseSigmaRec[nrow(Z),k,l,j]))
           for(i in (1:nrow(Z)))
          { 
            rmseMedRec[i,k,l,j] = norm(mu_hatIter[i] - mu1,"2")
            rmseSigmaRec[i,k,l,j] = norm(SigmaIter[i,  ,] - Sigma1,"F")}
         
        }  
        
        
        faux_positifsRec[,k,l,j] = cumsum(outliers_labels == 1 & data$labelsVrais == 0)
        faux_negatifsRec[,k,l,j] = cumsum(outliers_labels == 0 & data$labelsVrais == 1)
        temps_calcul[k,l,j] = temps["elapsed"] 
        cat("k =", k, "l =", l, "j =", j, 
            "-- temps_calcul =", temps_calcul[k, l, j], "secondes\n")  
        print(paste0("k = ",k))
        print(paste0("l = ",l))
        tc <- table(data$labelsVrais[1:(nrow(Z))], as.numeric(outliers_labels)[1:(nrow(Z))])
        if ("0" %in% rownames(tc)){
          if((tc["0","0"] + tc["0","1"]) != 0)
          {fp  <-(tc["0", "1"]/(tc["0", "1"] + tc["0", "0"]))*100
          print(paste0("FP ",fp))}
        }
        
      #  print(paste0("FP ",faux_positifsRec[1e4,k,l,j]))    
        if (length(unique(data$labelsVrais)) == 2) {
          auc <- as.numeric(pROC::auc(pROC::roc(data$labelsVrais, distances))) * 100
        } else {
          auc <- 50  # Valeur par défaut pour un cas non exploitable
        }
        aucRec[k,l,j] = auc
        distancesRec[,k,l,j] = distances
        outliersLabelsRec[,k,l,j] = outliers_labels
        #taux_OutliersDetectesVraisRec[,k,l,j] = cumulativeOutlierDetection(resultats,data,pourcentage = r,"Shifted Gaussian Contamination scenario")$taux_outliers_detectes_vrais 
      }
      
      #print(fp[compt])
      #correction = qchisq(.5,df = d)/(median(sqrt(distances)))^2
      
      # if(m == "offline"){
      # distances_corr = correction*distances}
      # if (m == "online" || m == "streaming"){
      #   facteurs = correctionDistanceMahalanobis(distances,Z,methode = "streaming")
      #  distances_corr = hadamard.prod(facteurs,distances)
      # }
      # outliers_corr = rep(0,nrow(Z))
      # #cutoff = qchisq(.95,df = d)
      # cutoff = mahalanobis_cutoff(Sigma,mu_hat)
      # distances_corr = resultats$distances
      # for (i in (1:nrow(Z))){
      # if(distances_corr[i] > cutoff) {outliers_corr[i] = 1}
      # }
      #   #print("ok Mahalanobis")
      #   tc <- table(data$labelsVrais[1:nrow(Z)], as.numeric(outliers_corr))
      #   tc <- safe_access_tc(tc)
      #   
      #   
      #   if ((tc["0", "0"] + tc["0", "1"]) != 0) {
      #     if (m == "offline")  {fp_offline_corr[compt] = (tc["0", "1"] / (tc["0", "1"] + tc["0", "0"])) * 100}
      #     if (m == "online")  {fp_online_corr[compt] = (tc["0", "1"] / (tc["0", "1"] + tc["0", "0"])) * 100}
      #     if (m == "streaming")  {fp_streaming_corr[compt] = (tc["0", "1"] / (tc["0", "1"] + tc["0", "0"])) * 100}
      #   } 
      #   
      #      
      #   #print(fp_corr[compt])
      # 
      #compt_meth = compt_meth + 1
    }
    

      
    
  }
  
  
  return (list(rmseSigmaRec= rmseSigmaRec,rmseMedRec = rmseMedRec,temps_calcul = temps_calcul,faux_positifsRec = faux_positifsRec, faux_negatifsRec = faux_negatifsRec,taux_OutliersVraisRec = taux_OutliersVraisRec,taux_OutliersDetectesVraisRec  = taux_OutliersDetectesVraisRec,taux_OutliersDetectesRec = taux_OutliersDetectesRec,aucRec = aucRec,distancesRec = distances,outliersLabelsRec = outliersLabelsRec,labelsVraisRec = labelsVraisRec))
}
