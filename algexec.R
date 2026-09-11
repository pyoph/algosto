# ============================================================
# APPEL DE L'ALGORITHME POUR UNE STATION
# ============================================================

# ------------------------------------------------------------
# 1. Choisir la station
# ------------------------------------------------------------

st <- "Aotizhongxin"

# ------------------------------------------------------------
# 2. Récupérer la matrice X
#    (NA remplacés par 0)
# ------------------------------------------------------------

X <- X_list[[st]]

# Vérification
cat("Station :", st, "\n")
cat("Dimensions de X :", nrow(X), "x", ncol(X), "\n")
cat("Nombre de NA dans X :", sum(is.na(X)), "\n")


# ------------------------------------------------------------
# 3. Construire S à partir des données RAW
#
# S[i,j] = 1 si la valeur était observée
# S[i,j] = 0 si la valeur était NA dans raw
# ------------------------------------------------------------

raw_st <- raw[
  station == st,
  ..vars
]

S <- !is.na(
  as.matrix(raw_st)
)

S <- matrix(
  as.integer(S),
  nrow = nrow(S),
  ncol = ncol(S)
)

# Vérification
cat("Nombre de valeurs observées :", sum(S), "\n")
cat("Nombre de valeurs manquantes :", sum(S == 0), "\n")


# ------------------------------------------------------------
# 4. Initialisation m0
# ------------------------------------------------------------

m0 <- rep(
  0,
  ncol(X)
)


# ------------------------------------------------------------
# 5. Lancer l'algorithme de l'article
#
# gamma_n = gamma0 * n^(-alpha)
# ------------------------------------------------------------

set.seed(123)

result <- geometric_median_SA(
  X = X,
  S = S,
  gamma0 = 1,
  alpha = 0.8,
  m0 = m0
)


# ------------------------------------------------------------
# 6. Récupérer les estimateurs finaux
# ------------------------------------------------------------

m_final <- result$final_recursive

mbar_final <- result$final_averaged


# ------------------------------------------------------------
# 7. Moyenne empirique de X
# ------------------------------------------------------------

mean_X <- colMeans(X)


# ------------------------------------------------------------
# 8. Affichage
# ------------------------------------------------------------

cat("\n========================================\n")
cat("STATION :", st, "\n")
cat("========================================\n")

cat("\nMoyenne empirique :\n")
print(mean_X)

cat("\nMédiane géométrique récursive m_n :\n")
print(m_final)

cat("\nMédiane géométrique moyennée mbar_n :\n")
print(mbar_final)


# ------------------------------------------------------------
# 9. Tableau de comparaison
# ------------------------------------------------------------

comparison <- data.frame(
  variable = vars,
  mean = mean_X,
  recursive = m_final,
  averaged = mbar_final
)

print(comparison)


# ============================================================
# PLOTS : MOYENNE vs GEOMETRIC MEDIAN
# ============================================================

# ------------------------------------------------------------
# 1. Comparaison des 3 vecteurs
# ------------------------------------------------------------

plot(
  mean_X,
  type = "b",
  pch = 16,
  lty = 1,
  xaxt = "n",
  xlab = "Variable",
  ylab = "Valeur",
  main = paste(
    "Empirical mean vs geometric median -",
    st
  )
)

lines(
  m_final,
  type = "b",
  pch = 1,
  lty = 2
)

lines(
  mbar_final,
  type = "b",
  pch = 2,
  lty = 3
)

axis(
  1,
  at = 1:length(vars),
  labels = vars,
  las = 2,
  cex.axis = 0.7
)

legend(
  "topright",
  legend = c(
    "Empirical mean",
    "Recursive m_n",
    "Averaged mbar_n"
  ),
  lty = c(1, 2, 3),
  pch = c(16, 1, 2)
)


# ============================================================
# 2. Distance entre m_n et la moyenne empirique
# ============================================================

dist_recursive_mean <- sqrt(
  sum(
    (m_final - mean_X)^2
  )
)

dist_averaged_mean <- sqrt(
  sum(
    (mbar_final - mean_X)^2
  )
)

cat(
  "Distance m_n - moyenne =",
  dist_recursive_mean,
  "\n"
)

cat(
  "Distance mbar_n - moyenne =",
  dist_averaged_mean,
  "\n"
)


# ============================================================
# 3. Comparaison variable par variable
# ============================================================

comparison <- data.frame(
  variable = vars,
  mean = mean_X,
  recursive = m_final,
  averaged = mbar_final
)

matplot(
  t(
    comparison[, c(
      "mean",
      "recursive",
      "averaged"
    )]
  ),
  type = "b",
  pch = c(16, 1, 2),
  lty = c(1, 2, 3),
  xaxt = "n",
  xlab = "Variable",
  ylab = "Valeur",
  main = paste(
    "Comparison -",
    st
  )
)

axis(
  1,
  at = 1:length(vars),
  labels = vars,
  las = 2,
  cex.axis = 0.7
)

legend(
  "topright",
  legend = c(
    "Empirical mean",
    "Recursive m_n",
    "Averaged mbar_n"
  ),
  lty = c(1, 2, 3),
  pch = c(16, 1, 2)
)