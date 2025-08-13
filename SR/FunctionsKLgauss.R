KL <- function(parms1, parms2){
  invSigma2 <- solve(parms2$Sigma)
  0.5*(log(det(parms2$Sigma)/det(parms1$Sigma)) - d + sum(diag(invSigma2%*%parms1$Sigma)) +
         t(parms2$mu-parms1$mu)%*%invSigma2%*%(parms2$mu-parms1$mu))[1, 1]
}

FrobeniusNormError = function(Sigmahat,Sigma){
  return (norm(Sigmahat - Sigma,"F"))
}

taux_detection <- function(vrais_labels, labels_predits) {
  # S'assurer que les entrées sont binaires
  if (!all(vrais_labels %in% c(0, 1)) || !all(labels_predits %in% c(0, 1))) {
    stop("Les vecteurs doivent contenir uniquement 0 (normal) ou 1 (outlier).")
  }
  
  # Calcul des éléments de la matrice de confusion
  VP <- sum(vrais_labels == 1 & labels_predits == 1)  # Vrai Positifs
  FP <- sum(vrais_labels == 0 & labels_predits == 1)  # Faux Positifs
  FN <- sum(vrais_labels == 1 & labels_predits == 0)  # Faux Négatifs
  VN <- sum(vrais_labels == 0 & labels_predits == 0)  # Vrai Négatifs
  
  # Taux
  TPR <- if ((VP + FN) > 0) VP / (VP + FN) else NA  # Sensibilité
  FPR <- if ((FP + VN) > 0) FP / (FP + VN) else NA  # 1 - Spécificité
  
  return(list(
    TPR = TPR,
    FPR = FPR,
    VP = VP,
    FP = FP,
    FN = FN,
    VN = VN
  ))
}

# Contamination parms: F1
ParmsF1 <- function(m1, k1, l1, rho1){
  d <- length(m1)
  mu1 <- k1*m1
  sigmaSq1 <- l1*sigmaSq0
  Sigma1 <- diag(sqrt(sigmaSq1)) %*% toeplitz(rho1^(0:(d-1))) %*% diag(sqrt(sigmaSq1))
  return(list(mu1=mu1, Sigma1=Sigma1))
}

