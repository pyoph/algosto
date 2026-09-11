library(dplyr)
library(ggplot2)

resultats_erreurs <- list()

for (station in stations) {
  
  nom_station <- tools::file_path_sans_ext(station)
  
  for (m in methodes) {
    
    setwd(critMissing)
    
    critFile <- paste0(
      "Crit-",
      m,"-",
      nom_station,
      ".RData"
    )
    
    if (!file.exists(critFile)) {
      next
    }
    
    load(critFile)
    
    # --------------------------------------------------------
    # Geometric median : UNIQUEMENT l'estimateur moyenné
    # --------------------------------------------------------
    
    if (m == "med_geom") {
      
      resultats_erreurs[[length(resultats_erreurs) + 1]] <-
        data.frame(
          station = nom_station,
          methode = "Geometric median",
          erreur = crit$erreur_avg
        )
      
    } else {
      
      # ------------------------------------------------------
      # Autres méthodes : estimateur standard
      # ------------------------------------------------------
      
      resultats_erreurs[[length(resultats_erreurs) + 1]] <-
        data.frame(
          station = nom_station,
          methode = m,
          erreur = crit$erreur_m
        )
    }
  }
}

erreurs <- bind_rows(resultats_erreurs)

print(erreurs)


# ============================================================
# GRAPHIQUE
# ============================================================
ggplot(
  erreurs,
  aes(
    x = station,
    y = erreur,
    color = methode,
    group = methode
  )
) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c(
      "Geometric median" = "red",
      "OGK" = "brown",
      "MCD" = "black",
      "VardiZhang" = "darkgreen"
    )
  ) +
  scale_y_log10() +
  labs(
    x = "Station",
    y = "Quadratic error",
    color = "Method"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )


for (station in stations) {
  
  nom_station <- tools::file_path_sans_ext(station)
  
  setwd(critMissing)
  
  critFile <- paste0(
    "Crit-",
    "med_geom","-",
    nom_station,
    ".RData"
  )
  
  load(critFile)
  
  figFile <- 
    paste0(
      "Estimation-error-",
      nom_station,
      ".pdf"
    )
  
  
  setwd(fig_missing)
  
  pdf(
    file = figFile,
    width = 8,
    height = 6
  )
  
  plot(
    crit$erreur_iter,
    type = "l",
    col = "red",
    lwd = 2,
    xlab = "Iteration",
    main = paste(
      "Estimation error -",
      nom_station
    )
  )
  
  dev.off()
  
  cat(
    "Figure enregistrée :",
    figFile,
    "\n"
  )
}