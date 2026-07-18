
compute_criteres = function(variance,outlab,distances,labels_vrais,SigmaTrue = Sigma0,r){
  
  erreurFrob <- norm(variance - SigmaTrue, "F")
  
  
  conf_mat <- table(Prediction = outlab, Vrai = labels_vrais)
  
  FP <- sum(outlab == 1 & labels_vrais == 0)
  FN <- sum(outlab == 0 & labels_vrais == 1)
  
  
  
  
  
  
  prop_hors_diag <- (FP + FN) / sum(conf_mat)
  
  
  if(r != 0){
    
    # ARI
    ari <- adjustedRandIndex(
      labels_vrais,
      outlab
    )
    
    
    #ari = ARI_manual(labels_vrais,outlab)  
    #print(paste0("ARI ",ari))
    # AUC
    auc_val <- as.numeric(auc(roc(labels_vrais, distances, direction='<')))
    
    #auc_val <- auc_manual(as.numeric(distances),labels_vrais)$auc
    
    #print(paste0("AUC ",auc_val))
    
  }
  else{
    auc_val = .5
    ari = 0
  }
  
  return(list(
    erreurFrob = erreurFrob,
    FP = FP,
    FN = FN,
    ARI = ari,
    AUC = auc_val,
    prop_hors_diag = prop_hors_diag
  ))
}





