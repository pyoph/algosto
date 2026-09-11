# Gradient partiel
gradient_partiel <- function(x, s, h, eps = 1e-12) {
  proj <- s * (x - h)          # P_S(X - h)
  norme <- sqrt(sum(proj^2))
  if (norme < eps) return(rep(0, length(h)))
  proj / norme
}

# Algorithme récursif + moyenné (Polyak-Ruppert)
algo_geometric_median <- function(X, S, gamma_const = sqrt(ncol(X)), alpha = 0.8,
                                  max_iter = nrow(X), m0 = NULL,
                                  verbose = FALSE) {
  n <- nrow(X)
  d <- ncol(X)
  miter <- matrix(0, nrow = n, ncol = d)  
  
  if (is.null(m0)) {
    m <- rep(0, d)
  } else {
    m <- m0
  }
  m_bar <- m
  
  for (t in 1:max_iter) {
    # Échantillonnage aléatoire d'une observation (avec remise)
    idx <- sample.int(n, 1)
    x <- X[idx, ]
    s <- S[idx, ]
    
    # Pas décroissant
    gamma_t <- gamma_const / t^alpha
    
    # Gradient
    grad <- gradient_partiel(x, s, m)
    
    # Mise à jour récursive
    m_new <- m + gamma_t * grad
    
    # Moyenne de Polyak-Ruppert
    m_bar <- m_bar + (m_new - m_bar) / (t + 2)
    
    m <- m_new
    miter[t,] = m_bar
    if (verbose && t %% 1000 == 0) {
      cat("Itération", t, "\n")
    }
  }
  
  list(m = m, m_bar = m_bar,miter = miter)
}


library(robustbase)  # pour covMcd, covOGK, scaleTau2
library(Gmedian)     # pour Gmedian (Vardi-Zhang)

setwd(rep_stations)

# Liste des fichiers de stations (par exemple .RData contenant Z et S)
stations <- list.files(pattern = "\\.RData$")

log_file <- "erreurs.log"
if (file.exists(log_file)) file.remove(log_file)

for (station in stations) {
  
  
  nom_station <- tools::file_path_sans_ext(station)
  cat("\n=== Traitement de", nom_station, "===\n")
  
  setwd(rep_stations)
  # suppose que Z et S sont chargés
  # Charger les données de la station
  load(station)
  
  
  
  # ------------------------------------------------------------
  # Construire la base complète de référence
  # ------------------------------------------------------------
  
  lignes_completes <- rowSums(S) == ncol(S)
  
  Z_complete <- Z[lignes_completes, , drop = FALSE]
  
  S_complete <- matrix(
    1,
    nrow = nrow(Z_complete),
    ncol = ncol(Z_complete)
  )
  
  cat(
    "Nombre de lignes complètes :", 
    nrow(Z_complete), 
    "\n"
  )
  
  setwd(resGeomMed)
  
  # --- 1) Médiane géométrique  ---
  fitFile <- paste0("Fit-med_geom-", nom_station, ".RData")
  tryCatch({
    res <- algo_geometric_median(Z, S, verbose = FALSE)
    
    
    res_ref <- algo_geometric_median(
      Z_complete,
      S_complete,
      verbose = FALSE
    )
    
    
    resultats <- list(m = res$m, m_bar = res$m_bar, miter = res$miter  , m_ref = res_ref$m)
    
    save(resultats, file = fitFile)
    
    
    cat("OK : med_geom\n")
  }, error = function(e) {
    msg <- paste(Sys.time(), "- Erreur med_geom -", nom_station, ":", conditionMessage(e))
    write(msg, file = log_file, append = TRUE)
    cat("Erreur med_geom :", conditionMessage(e), "\n")
  })
  
  # --- 2) OGK ---
  fitFile <- paste0("Fit-OGK-", nom_station, ".RData")
  tryCatch({

    res_ogk <- covOGK(Z, sigmamu = scaleTau2)
    
    
    # --- Référence sur données complètes ---
    res_ogk_ref <- covOGK(
      Z_complete,
      sigmamu = scaleTau2
    )
    resultats <- list(m = res_ogk$center, m_bar = res_ogk$center,    m_ref = res_ogk_ref$center)
    
    save(resultats, file = fitFile)
    cat("OK : OGK\n")
  }, error = function(e) {
    msg <- paste(Sys.time(), "- Erreur OGK -", nom_station, ":", conditionMessage(e))
    write(msg, file = log_file, append = TRUE)
    cat("Erreur OGK :", conditionMessage(e), "\n")
  })
  
  # --- 3) MCD ---
  fitFile <- paste0("Fit-MCD-", nom_station, ".RData")
  tryCatch({
    # Imputation simple par médiane
   
    res_mcd <- covMcd(Z)
    
    res_mcd_ref <- covMcd(
      Z_complete
    )
    
    resultats <- list(m = res_mcd$center, m_bar = res_mcd$center,    m_ref = res_mcd_ref$center)
    
    save(resultats, file = fitFile)
    cat("OK : MCD\n")
  }, error = function(e) {
    msg <- paste(Sys.time(), "- Erreur MCD -", nom_station, ":", conditionMessage(e))
    write(msg, file = log_file, append = TRUE)
    cat("Erreur MCD :", conditionMessage(e), "\n")
  })
  
  # --- 4) Vardi-Zhang (médiane géométrique classique) ---
  fitFile <- paste0("Fit-VardiZhang-", nom_station, ".RData")
  tryCatch({

    med_vz <- Weiszfeld(Z)
    
    med_vz_ref <- Weiszfeld(
      Z_complete
    )
    
    resultats <- list(m = med_vz$median, m_bar = med_vz$median,    m_ref =   med_vz_ref$median)
    
    save(resultats, file = fitFile)
    cat("OK : Vardi-Zhang \n")
  }, error = function(e) {
    msg <- paste(Sys.time(), "- Erreur Vardi-Zhang -", nom_station, ":", conditionMessage(e))
    write(msg, file = log_file, append = TRUE)
    cat("Erreur Vardi-Zhang :", conditionMessage(e), "\n")
  })
  
  # Nettoyer
  rm(Z, S)
}
