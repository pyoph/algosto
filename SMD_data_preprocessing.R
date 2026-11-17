url =  "https://www.kaggle.com/api/v1/datasets/download/mgusat/smd-onmiad"


options(timeout = 600)  # 10 minutes

dest <- "archive.zip"
setwd(smd_data_dir)

download.file(url, destfile = dest, mode = "wb")

unzip(dest)

# download all SMD files

files <- list.files(".", recursive = TRUE, full.names = TRUE)

train_files  <- grep("train/machine", files, value = TRUE)

test_files  <- grep("test/machine", files, value = TRUE)
label_files <- grep("test_label/machine", files, value = TRUE)



length(train_files)
length(label_files)


load_smd_machine <- function(test_path, label_path) {
  
  Z <- as.matrix(
    read.table(test_path, sep = ",")
  )
  
  labels <- scan(label_path)
  
  list(
    Z = Z,
    labels = labels,
    name = basename(test_path)
  )
}


remove_high_corr <- function(Z, threshold = 0.7)
{
  C <- cor(Z)
  
  C[lower.tri(C, diag = TRUE)] <- 0
  
  remove <- which(apply(abs(C) > threshold, 2, any))
  
  if(length(remove) > 0)
  {
    cat("Removing cols:", remove, "\n")
  }
  
  Z[, -remove, drop = FALSE]
  return(Z)
}



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

nbmachines <- length(smd_data)

batch_size = 10

for (j in 1:nbmachines) {
  
  cat("\n====================\n")
  cat("Machine", j, "\n")
  cat("====================\n")
  
  Z <- smd_data[[j]]$Z
  labels <- smd_data[[j]]$labels
  
  Z <- as.matrix(Z)
  
  # Keep the first 23,000 observations 
  Z <- Z[1:23000, ]
  labels <- labels[1:23000]
  
  # ============================================================
  # 1. Variable filtering
  # ============================================================
  tol <- 1e-2
  
  # Keep columns that:
  # - do not contain non-finite values (NA, NaN, Inf)
  # - have more than one unique value
  # - have an empirical standard deviation > tol
  # - have a MAD > tol
  
  valid_cols <- apply(Z, 2, function(col) {
    if (any(!is.finite(col))) return(FALSE)
    if (length(unique(col)) <= 1) return(FALSE)
    if (sd(col, na.rm = TRUE) < tol) return(FALSE)
    if (mad(col, na.rm = TRUE) < tol) return(FALSE)
    return(TRUE)
  })
  
  Z <- Z[, valid_cols, drop = FALSE]
  
  # ============================================================
  # Redundancy removal (correlation threshold = 0.7) and linear dependent columns
  # ============================================================
  
  Z <- remove_high_corr(
    Z,
    threshold = 0.7
  )
  

  ###################Detect linearity#######
  
  # QR with pivoting to detect linearly dependent columns
  qr_Z <- qr(Z, tol = 1e-2)
  
  
  
  # Keep only the columns that are linearly independent
  keep_cols <- sort(qr_Z$pivot[1:qr_Z$rank])
  Z = Z[,keep_cols]
  
  
  
  # ============================================================
  # 3. Robust standardization (median and MAD)
  # ============================================================
  # Remove too low MAD
  tol <- 1e-2
  
  # Identifier les colonnes à supprimer (MAD ou sd ~ 0, ou non-finis)
  mad_vals <- apply(Z, 2, mad, na.rm = TRUE)
  sd_vals  <- apply(Z, 2, sd,  na.rm = TRUE)
  
  valid_cols <- which(
    is.finite(mad_vals) & mad_vals > tol &
    is.finite(sd_vals)  & sd_vals  > tol
  )
  
  
  Z <- Z[, valid_cols, drop = FALSE]
  

  # SRobust standardization
  Z <- apply(Z, 2, function(x) {
    med <- median(x, na.rm = TRUE)
    mad_val <- mad(x, na.rm = TRUE)
    (x - med) / mad_val
  })
  
  Z <- as.matrix(Z)
  
  if (any(!is.finite(Z))) {
    warning(paste("Machine", j, ": non-finite values detected after standardization. Replacing by 0."))
    Z[!is.finite(Z)] <- 0
  }  
  # ============================================================
  # Extract uncontaminated data (for covariance reference)
  # ============================================================
  Z_clean <- Z[labels == 0, , drop = FALSE]
  
  # ============================================================
  # Verification of dimensions
  # ============================================================
  cat("\n--- Vérification des dimensions ---\n")
  cat("Z       :", nrow(Z), "lignes x", ncol(Z), "colonnes\n")
  cat("labels  :", length(labels), "éléments\n")
  cat("Z_clean :", nrow(Z_clean), "lignes x", ncol(Z_clean), "colonnes\n")
  
  if(nrow(Z) != 23000) stop("Erreur machine ", j, " : Z ne contient pas 23000 lignes.")
  if(length(labels) != 23000) stop("Erreur machine ", j, " : labels ne contient pas 23000 éléments.")
  
  # ============================================================
  # Save
  # ============================================================
  data_smd_mach <- paste0("data_machine-", j, ".RData")
  
  save(
    Z,
    Z_clean,
    labels,
    batch_size,
    file = data_smd_mach
  )
  
  
}
