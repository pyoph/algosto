offlineComputationTimes <- function(Z,
                                    Ninit = 1e3,
                                    batch = 10,
                                    machine = 2,
                                    repRes = res_SMD,cm = 2) {
  
  n <- nrow(Z)
  
  # Tailles des échantillons :
  # Ninit, Ninit + batch, Ninit + 2*batch, ..., n
  sizes <- seq(
    from = Ninit,
    to = n,
    by = batch
  )
  
  # Ajouter n si la dernière étape ne tombe pas exactement sur n
  if (tail(sizes, 1) != n) {
    sizes <- c(sizes, n)
  }
  
  # Temps cumulés, enregistrés observation par observation
  cum_times_MCD <- numeric(n)
  cum_times_OGK <- numeric(n)
  
  # Temps cumulés au niveau des batches
  cumulative_MCD <- 0
  cumulative_OGK <- 0
  
  #outliers labels
  
  outl_labels_mcd <- rep(0,n)
  outl_labels_ogk <- rep(0,n)
  
  # Première observation du batch courant
  start_idx <- 1
  
  for (i in seq_along(sizes)) {
    
    ni <- sizes[i]
    
    # On prend les lignes 1,...,ni
    Zi <- Z[1:ni, , drop = FALSE]
    
    # -----------------
    # MCD
    # -----------------
    t0 <- Sys.time()
    mcd <- covMcd(Zi)
    dist_mcd = mcd$mah
    cutoff = calcule_cutoff(Zi,type = "quantcorr",n = nrow(Zi),c_m = cm)
    outl_labels_mcd[1:ni] = as.numeric(dist_mcd > cutoff)
    
    
    current_time_MCD <- as.numeric(Sys.time() - t0, units = "secs")
    
    # -----------------
    # OGK
    # -----------------
    t0 <- Sys.time()
    ogk <- covOGK(Zi, sigmamu = scaleTau2)
    dist_ogk = ogk$distances
    cutoff = calcule_cutoff(Zi,type = "quantcorr",n = nrow(Zi))
    outl_labels_ogk[1:ni] = as.numeric(dist_ogk > cutoff)
    current_time_OGK <- as.numeric(Sys.time() - t0, units = "secs")
    
    # -----------------
    # Temps cumulés
    # -----------------
    cumulative_MCD <- cumulative_MCD + current_time_MCD
    cumulative_OGK <- cumulative_OGK + current_time_OGK
    
    # Toutes les observations correspondant à ce batch
    # reçoivent le même temps cumulé
    cum_times_MCD[start_idx:ni] <- cumulative_MCD
    cum_times_OGK[start_idx:ni] <- cumulative_OGK
    
    # Batch suivant
    start_idx <- ni + 1
  }
  
  setwd(repRes)
  
  fitFile = paste0("Fit-MCDonlineQC-","machine-",machine,".RData")
  
  resultats = list(outliers_labels = outl_labels_mcd,cum_times = cum_times_MCD)
  
  save(resultats,file = fitFile)
  
  fitFile = paste0("Fit-OGKonlineQC-","machine-",machine,".RData")
  
  resultats = list(outliers_labels = outl_labels_ogk,cum_times = cum_times_OGK)
  
  save(resultats,file = fitFile)
  
  
  
  result <- data.frame( N = 1:n, time_MCD = cum_times_MCD, time_OGK = cum_times_OGK, outlMCD = outl_labels_mcd, outlOGK = outl_labels_ogk ) 
  return(result)
}


# -----------------------------
# Lancer le test
# -----------------------------
times <- offlineComputationTimes(
  Z,
  Ninit = 1000,
  batch = 10
)

# Afficher les résultats
print(times)

