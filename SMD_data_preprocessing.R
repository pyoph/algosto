smd_data_dir = "~/smd_data_dir"

url <- "https://github.com/NetManAIOps/OmniAnomaly/archive/refs/heads/master.zip"
dest <- "smd.zip"
setwd(smd_data_dir)

download.file(url, destfile = dest, mode = "wb")

unzip(dest)
# chercher tous les fichiers SMD
files <- list.files(".", recursive = TRUE, full.names = TRUE)

train_files  <- grep("train/machine", files, value = TRUE)

test_files  <- grep("test/machine", files, value = TRUE)
label_files <- grep("test_label/machine", files, value = TRUE)



length(train_files)
length(label_files)


smd_data <- list()

# Train
for (i in seq_along(train_files)) {
  
  smd_data[[i]] <- load_smd_machine(
    train_files[i],
    label_files[i]
  )
  
}

# Test
offset <- length(train_files)

for (i in seq_along(test_files)) {
  
  smd_data[[offset + i]] <- load_smd_machine(
    test_files[i],
    label_files[i]
  )
  
}

nbmachines = length(smd_data)

keep_columns_ok = list()

for(j in 20:nbmachines){
  
  cat("\n====================\n")
  cat("Machine", j, "\n")
  cat("====================\n")
  
  
  Z <- smd_data[[j]]$Z
  labels <- smd_data[[j]]$labels
  Z <- as.matrix(Z)
  
  # #Filtre 1 : enlever les colonnes qui font planter les méthodes
  res <- tryCatch({
    #   
    clean_cols <- sequential_clean_columns(Z)
    #   
    list(
      Z = clean_cols$Z,
      keep = clean_cols$keep
    )
    #   
  }, error = function(e) {
    #   
    cat("Erreur détectée -> retry sans colonne 1\n")
    #   
    Z2 <- Z[, -1, drop = FALSE]
    #   
    clean_cols <- sequential_clean_columns(Z2)
    #   
    list(
      Z = clean_cols$Z,
      keep = clean_cols$keep + 1
    )
  })
  # 
  Z <- res$Z
  keep_columns_ok[[j]] <- res$keep
  # 
  # 
  # 
  # #Filtre 2 : enlever les colonnes avec MAD ~ 0
  # 
  # 
  tol = 1e-3
  # 
  mad_cols <- apply(Z, 2, mad, na.rm = TRUE)
  # 
  keep <- which(mad_cols > tol)
  # 
  removed <- setdiff(1:ncol(Z), keep)
  # 
  if (length(removed) > 0) {
    cat("Colonnes retirées (MAD ~ 0):", removed, "\n")
  }
  
  keep_columns_ok[[j]] = keep
  
  
  Z = Z[, keep_columns_ok[[j]], drop = FALSE]
  
  
  
  #Z = Z[,keep_columns_ok[[j]],drop = FALSE]
  col_before <- colnames(Z)
  
  Z = remove_high_corr(Z,threshold = 0.99)
  col_after <- colnames(Z)
  keep_idx <- match(col_after, col_before)
  
  # mise à jour des colonnes globales gardées
  keep_columns_ok[[j]] <- keep_columns_ok[[j]][keep_idx]
  #Z = Z[,keep_columns_ok[[j]],drop = FALSE]
  labels <- labels[1:nrow(Z)]
  Z_clean = Z[labels == 0,]
  
  data_smd_mach = paste0("data_machine-",j,".RData")
  
  save(Z,Z_clean,labels,file = data_smd_mach)
}