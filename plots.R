################################################################
# 0) CONFIGURATION GLOBALE
################################################################

scenarios  <- c(scenarios_1_param)     # ou c(scenarios_1_param, scenarios_2_param)
lwd_value  <- 4
alpha_val  <- 0.8             # transparence globale

################################################################
# 1) MÉTADONNÉES DES MÉTHODES (courbes)
################################################################


methods_df <- data.frame(
  method = c(
    # --- QC ---
    "SampleNaiveQuantonlinecorr",
    "OnlineUsQuantonlinecorr",
    "StreamingUsonlineQuantcorr",
    "OfflinewithQuantcorr",
    "OGKQC",
    "MCDQC",
    "OracleQC",
    # --- RD ---
    "SampleNaivewithoutonlinequantilecorr",
    "OnlineUswithoutQuantonlinecorr",
    "StreamingUswithoutQuantonlinecorr",
    "OfflineUswithoutQuantcorr",
    "OGKRD",
    "MCDRD",
    "OracleRD",
    # --- Raw ---
    "SampleRaw",
    "OnlRaw",
    "StrmRaw",
    "OfflRaw",
    "OGK",
    "MCD",
    "Oracle"
  ),
  group = c(rep("QC",  7),
            rep("RD",  7),
            rep("Raw", 7)),
  color = c(
    # QC
    "darkgreen", "#FF6B6B", "red", "purple4", "brown", "black", "blue",
    # RD
    "darkgreen", "#FF6B6B", "red", "purple4", "brown", "black", "blue",
    # Raw
    "darkgreen", "#FF6B6B", "red", "purple4", "brown", "black", "blue"
  ),
  pch = c(
    rep(8,  7),   # QC  : étoiles
    rep(15, 7),   # RD  : carrés
    rep(17, 7)    # Raw : triangles
  ),
  lty = c(
    # QC
    1, 1, 1, 5, 5, 5, 3,
    
    # RD
    1, 1, 1, 5, 5, 5, 3,
    # Raw
    1, 1, 1, 5, 5, 5, 3
  ),
  stringsAsFactors = FALSE
)

methods_df$color_alpha <- adjustcolor(methods_df$color, alpha.f = alpha_val)
rownames(methods_df)  <- methods_df$method

# Sous-ensembles (construits depuis le df, plus de vecteurs parallèles)
subset_frob <- methods_df[methods_df$method %in%
                            methods_df$method[methods_df$group %in% c("QC", "RD", "Raw")][1:6], ]
# plus lisible : liste explicite
methods_frob <- c("SampleNaiveQuantonlinecorr",
                  "OnlineUsQuantonlinecorr",
                  "StreamingUsonlineQuantcorr",
                  "OfflinewithQuantcorr",
                  "OGK", "MCD")

methodList <- c("SampleNaivewithoutonlinequantilecorr",
                "OfflineUswithoutQuantcorr",
                "OnlineUswithoutQuantonlinecorr",
                "StreamingUswithoutQuantonlinecorr",
                "MCD", "OGK", "Oracle")

methods_auc <- c(methods_frob, "Oracle")

subset_frob <- methods_df[methods_df$method %in% methods_frob, ]
subset_list <- methods_df[methods_df$method %in% methodList,  ]
subset_auc  <- methods_df[methods_df$method %in% methods_auc, ]

# Vérifications
stopifnot(nrow(methods_df) == 21,
          !any(is.na(methods_df$color)),
          !any(is.na(methods_df$pch)))

# Index (dans all_methodes = methods_df$method) des sous-ensembles :
idxFrob <- which(methods_df$method %in% methods_frob)
idxList <- which(methods_df$method %in% methodList)
idxAUC  <- which(methods_df$method %in% methods_auc)
all_methodes <- methods_df$method

################################################################
# 2) FONCTION DE TRAÇAGE GÉNÉRIQUE
################################################################
# trace une courbe pour chaque ligne de `subset_df`
# mat       : matrice [nR x length(all_methodes)]
# x         : abscisses
# subset_df : sous-ensemble de methods_df
# ylim, log, ... : passés à plot()

plot_lines <- function(mat, x, subset_df, rows = seq_along(x),
                       ylim, log = "", xaxt = "n", yaxt = "n",
                       xlab = "", ylab = "", lwd = lwd_value) {
  
  idx <- which(methods_df$method %in% subset_df$method)
  y   <- mat[rows, idx[1]]
  
  plot(x, y,
       type = "l", log = log, lwd = lwd,
       col  = methods_df$color_alpha[idx[1]],
       lty  = methods_df$lty[idx[1]],
       ylim = ylim, xaxt = xaxt, yaxt = yaxt,
       xlab = xlab, ylab = ylab)
  
  if (length(idx) > 1) {
    for (k in 2:length(idx)) {
      i <- idx[k]
      lines(x, mat[rows, i],
            lwd = lwd,
            col = methods_df$color_alpha[i],
            lty = methods_df$lty[i])
    }
  }
}
################################################################
# 3) BOUCLE SUR LES SCÉNARIOS
################################################################

for (sc in scenarios) {
  
  k    <- sc$k
  l    <- sc$l
  rho1 <- sc$rho1
  
  nR <- length(rList[1:13])
  
  erreursSigmaPlot  <- array(0, dim = c(nR, length(all_methodes)))
  faux_positifsPlot <- array(0, dim = c(nR, length(all_methodes)))
  faux_negatifsPlot <- array(0, dim = c(nR, length(all_methodes)))
  ariPlot           <- array(0, dim = c(nR, length(all_methodes)))
  aucPlot           <- array(0, dim = c(nR, length(all_methodes)))
  propHorsDiagPlot  <- array(0, dim = c(nR, length(all_methodes)))
  
  for (m in seq_along(rList[1:13])) {
    r <- rList[m]
    
    for (j in seq_along(all_methodes)) {
      methode <- all_methodes[j]
      
      setwd(criteres)
      critFile <- paste0('Crit-', methode, "-d", d, '-n', n,
                         '-k', k, '-l', l, '-rho', rho1,
                         '-r', r, '-mean', ".RData")
      load(critFile)
      
      if (methode %in% methods_frob) {
        erreursSigmaPlot[m, j]  <- crit_mean$erreurFrob
        faux_positifsPlot[m, j] <- crit_mean$FP
        if (r != 0) {
          aucPlot[m, j]           <- crit_mean$AUC
          faux_negatifsPlot[m, j] <- crit_mean$FN
          ariPlot[m, j]           <- crit_mean$ARI
          propHorsDiagPlot[m, j]  <- crit_mean$prop_hors_diag
        }
      }
      
      if (methode == "Oracle") {
        faux_positifsPlot[m, j] <- crit_mean$FP
        if (r != 0) {
          aucPlot[m, j]           <- crit_mean$AUC
          faux_negatifsPlot[m, j] <- crit_mean$FN
          propHorsDiagPlot[m, j]  <- crit_mean$prop_hors_diag
          ariPlot[m, j]           <- crit_mean$ARI
        }
      }
      
      if (methode %in% c(methods_frob,
                         "SampleNaivewithoutonlinequantilecorr",
                         "OnlineUswithoutQuantonlinecorr",
                         "StreamingUswithoutQuantonlinecorr",
                         "OfflineUswithoutQuantcorr",
                         "SampleRaw", "OnlRaw", "StrmRaw", "OfflRaw")) {
        if (r == 0) {
          faux_positifsPlot[m, j] <- crit_mean$FP
        } else {
          faux_negatifsPlot[m, j] <- crit_mean$FN
          faux_positifsPlot[m, j] <- crit_mean$FP
          ariPlot[m, j]           <- crit_mean$ARI
          propHorsDiagPlot[m, j]  <- crit_mean$prop_hors_diag
        }
      }
    }
  }
  
  ##############################################################
  # GRAPHIQUES
  ##############################################################
  setwd(figures)
  file <- paste0("scen-k", k, "-l", l, "-rho1", rho1, ".pdf")
  
  pdf(file, width = 25, height = 4)
  par(mfrow = c(1,5), mar = c(4,4,2,1))
  
  # 1) Frobenius -------------------------------------------------
  plot_lines(erreursSigmaPlot, rList[1:13], subset_frob,
             ylim = c(1e-1, 1e2), log = "y")
  axis(1, at = rList[1:13], las = 1, cex.axis = 1.8)
  axis(2, at = 10^seq(-1,2),
       labels = parse(text = paste0("10^", -1:2)),
       las = 1, cex.axis = 1.8)
  box()
  
  # 2) FN --------------------------------------------------------
  fn_rate <- function(i) faux_negatifsPlot[2:13, i] /
    ((rList[2:13]/100) * n) * 100
  plot(rList[2:13], fn_rate(idxList[1]), type = "l", lwd = lwd_value,
       col = methods_df$color_alpha[idxList[1]],
       lty = methods_df$lty[idxList[1]],
       ylim = c(0,100), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  for (k2 in 2:length(idxList)) {
    i <- idxList[k2]
    lines(rList[2:13], fn_rate(i), lwd = lwd_value,
          col = methods_df$color_alpha[i],
          lty = methods_df$lty[i])
  }
  axis(1, at = rList[-1], las = 1, cex.axis = 1.8)
  axis(2, las = 1, cex.axis = 1.8)
  box()
  
  # 3) FP --------------------------------------------------------
  fp_rate <- function(i) faux_positifsPlot[1:13, i] /
    ((1 - rList[1:13]/100) * n) * 100
  plot(rList[1:13], fp_rate(idxList[1]), type = "l", lwd = lwd_value,
       col = methods_df$color_alpha[idxList[1]],
       lty = methods_df$lty[idxList[1]],
       ylim = c(0,20), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  for (k2 in 2:length(idxList)) {
    i <- idxList[k2]
    lines(rList[1:13], fp_rate(i), lwd = lwd_value,
          col = methods_df$color_alpha[i],
          lty = methods_df$lty[i])
  }
  axis(1, at = rList[1:13], las = 1, cex.axis = 1.8)
  axis(2, las = 1, cex.axis = 1.8)
  box()
  
  # 4) AUC -------------------------------------------------------
  plot_lines(aucPlot[2:13,], rList[2:13], subset_auc, ylim = c(0,1))
  axis(1, at = rList[-1], las = 1, cex.axis = 1.8)
  axis(2, las = 1, cex.axis = 1.8)
  box()
  
  # 5) Prop hors diag --------------------------------------------
  prop_inv <- 1 - propHorsDiagPlot
  plot(rList[2:13], prop_inv[2:13, idxList[1]],
       type = "l", lwd = lwd_value,
       col = methods_df$color_alpha[idxList[1]],
       lty = methods_df$lty[idxList[1]],
       ylim = c(0,1), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  for (k2 in 2:length(idxList)) {
    i <- idxList[k2]
    lines(rList[2:13], prop_inv[2:13, i],
          lwd = lwd_value,
          col = methods_df$color_alpha[i],
          lty = methods_df$lty[i])
  }
  axis(1, at = rList[1:13], las = 1, cex.axis = 1.8)
  axis(2, las = 1, cex.axis = 1.8)
  box()
  
  dev.off()
}


################################################################
# 4) BOXPLOTS RELIÉS
################################################################

# --------------------------------------------------------------
# Métadonnées boxplot (6 méthodes)
# --------------------------------------------------------------


# On prend les 6 premières méthodes de la famille QC comme référence
# pour les noms courts : Sample naive, Online, Streaming, Offline, OGK, MCD
ref_methods <- c(
  "SampleNaiveQuantonlinecorr",
  "OnlineUsQuantonlinecorr",
  "StreamingUsonlineQuantcorr",
  "OfflinewithQuantcorr",
  "OGKQC",
  "MCDQC"
)

box_df <- data.frame(
  method      = c("Sample naive", "Online", "Streaming",
                  "Offline", "OGK", "MCD"),
  ref         = ref_methods,
  stringsAsFactors = FALSE
)

# Récupération des attributs visuels depuis methods_df
box_df$color       <- methods_df$color[match(box_df$ref, methods_df$method)]
box_df$color_alpha <- methods_df$color_alpha[match(box_df$ref, methods_df$method)]
box_df$lty         <- methods_df$lty[match(box_df$ref, methods_df$method)]
box_df$pch         <- 20   # seul élément propre au boxplot : le symbole des médianes

# --------------------------------------------------------------
# Préparation des données (matrices 100 x 6)
# --------------------------------------------------------------
list_configs <- list(
  t(temps_calcul_n1e4_d100),
  t(temps_calcul_n1e4_d10),
    t(temps_calcul_n1e5_d10)
)
for (m in seq_along(list_configs)) {
  colnames(list_configs[[m]]) <- box_df$method
}

labels_config <- c(
  expression(n == 10^4 * "," ~ d == 100),
  expression(n == 10^4 * "," ~ d == 10),
    expression(n == 10^5 * "," ~ d == 10)
)

n_meth     <- nrow(box_df)
method_gap <- 0.65
group_gap  <- 2.5

positions <- c(
  1:6  * method_gap,
  1:6  * method_gap + 6 * method_gap + group_gap,
  1:6  * method_gap + 2 * (6 * method_gap + group_gap)
)

group_centers <- c(mean(positions[1:6]),
                   mean(positions[7:12]),
                   mean(positions[13:18]))

pow_min <- -1; pow_max <- 2
y_ticks <- 10^(pow_min:pow_max)

# --------------------------------------------------------------
# PDF
# --------------------------------------------------------------
setwd(figures)
cairo_pdf("boxplot_temps_calcul.pdf", width = 12, height = 6)

par(mar = c(5,5,2,1), mgp = c(2.5, 0.8, 0))

# Groupe 1 (initialise le repère)
boxplot(list_configs[[1]],
        at = positions[1:6], col = box_df$color,
        names = rep("", 6), xaxt = "n", yaxt = "n",
        ylab = "", xlim = c(min(positions)-1, max(positions)+1),
        ylim = c(10^pow_min, 10^pow_max), log = "y",
        outline = FALSE, boxwex = 1.5)

# Groupes 2 et 3
for (g in 2:3) {
  boxplot(list_configs[[g]],
          at = positions[((g-1)*6+1):(g*6)],
          col = box_df$color, names = rep("", 6),
          xaxt = "n", yaxt = "n", add = TRUE,
          outline = FALSE, boxwex = 1.5)
}

# Lignes reliant les médianes
for (m in 1:n_meth) {
  medians <- sapply(list_configs, function(cfg) median(cfg[, m], na.rm = TRUE))
  x_m     <- positions[c(m, m + 6, m + 12)]
  
  lines(x_m, medians, col = box_df$color[m], lwd = 2, lty = 1)
  points(x_m, medians, pch = box_df$pch[m],
         col = box_df$color[m], cex = 1.4)
}

axis(1, at = group_centers, labels = labels_config,
     tick = FALSE, cex.axis = 1.1)
axis(2, at = y_ticks,
     labels = parse(text = paste0("10^", pow_min:pow_max)),
     las = 1, cex.axis = 1.1)
box()

dev.off()
