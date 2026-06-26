
####################################Scenarios##############################################"
scenarios <- list(
  #scenarios <- list(
  

  
  # =====================================================
  # DECENTRAGE + CONCENTRATION
  # =====================================================
  
  # Décentrage + concentration faible
  list(k = 5, l = 0.5, rho1 = 0.3),
  
  # Décentrage + concentration moyenne
  list(k = 10, l = 0.1, rho1 = 0.3),
  
  # Décentrage + concentration forte
  list(k = 50, l = 0.01, rho1 = 0.3),
  
  
  # =====================================================
  # DILATATION
  # =====================================================
  
  # Dilatation faible
  #list(k = 0, l = 2, rho1 = 0.3),
  
  # Dilatation moyenne
  #list(k = 0, l = 10, rho1 = 0.3),
  
  # Dilatation forte
  #list(k = 0, l = 100, rho1 = 0.3),
  
  
  # =====================================================
  # DECENTRAGE + DILATATION
  # =====================================================
  
  # Décentrage + dilatation faible
  list(k = 5, l = 2, rho1 = 0.3),
  
  # Décentrage + dilatation moyenne
  list(k = 10, l = 10, rho1 = 0.3),
  
  # Décentrage + dilatation forte
  list(k = 50, l = 100, rho1 = 0.3),
  
  
  
  
  
  # =====================================================
  # DEFORMATION + CONCENTRATION
  # =====================================================
  
  # Déformation + concentration faible
  list(k = 0, l = 0.5, rho1 = 0.5),
  
  # Déformation + concentration moyenne
  list(k = 0, l = 0.1, rho1 = 0.7),
  
  # Déformation + concentration forte
  list(k = 0, l = 0.01, rho1 = 0.95),
  
  
  # =====================================================
  # DEFORMATION + DILATATION
  # =====================================================
  
  # Déformation + dilatation faible
  list(k = 0, l = 2, rho1 = 0.5),
  
  # Déformation + dilatation moyenne
  list(k = 0, l = 10, rho1 = 0.7),
  
  # Déformation + dilatation forte
  list(k = 0, l = 100, rho1 = 0.95)
  #
  
  
  # =====================================================
  # CAS EXTREMES
  # =====================================================
  # 
  # # Anomalie extrême concentrée
  #list(k = 100, l = 0.01, rho1 = 0.99),
  # 
  # # Anomalie extrême dilatée
  #list(k = 100, l = 100, rho1 = 0.99)
  # 
)


#####################k, l and rho1 varying alone########################################

scenarios_1_param = list(
  list(k = 5, l = 1, rho1 = 0.3),
  

  list(k = 10, l = 1, rho1 = 0.3),
  
  list(k = 50, l = 1, rho1 = 0.3),
 
  list(k = 0, l = 2, rho1 = 0.3),
  
  list(k = 0, l = 10, rho1 = 0.3),
  
  list(k = 0, l = 1e2, rho1 = 0.3),
  list(k = 0, l = 1, rho1 = 0.5),
  
  list(k = 0, l = 1, rho1 = 0.7),
  
  list(k = 0, l = 1, rho1 = 0.95)
  
)

setwd("~/algosto")

save(
  scenarios,
  scenarios_k,
  scenarios_l,
  scenarios_rho1,
  file = "scenarios.RData"
)