# KL divergence btw F0 and F1

rm(list=ls()); par(pch=20, lwd=2)
library(mvtnorm)
source('~/algosto/FunctionsKLgauss.R')

# Dims
d <- 10; KLval <- c(0, 1, 5, 10, 25, 50, 100)

# Parms
#load(paste0('~/algosto/SimParmsGrid-d', d, '.Rdata'))
# l1val <- sort(l1val[which(l1val >= 1)])
# rho1val <- sort(rho1val[which(rho1val >= rho0)])
rbind(KLval, k1val, l1valup1, rho1val)
klNb <- length(KLval)

# Sigma0
Sigma0inv <- solve(Sigma0)
Sigma0eig <- eigen(Sigma0)
Sigma0invHalf <- Sigma0eig$vectors%*%diag(1/sqrt(Sigma0eig$values))%*%t(Sigma0eig$vectors)
Sigma0half <- Sigma0eig$vectors%*%diag(sqrt(Sigma0eig$values))%*%t(Sigma0eig$vectors)

#par(mfrow=c(2, 1))
B <- 1e6

klcol = c("black","darkgreen","blue","purple","orange","red","brown")


  
# Distribution (log scale)
U <- matrix(rnorm(d*B), B, d)
X0 <- rmvnorm(B, mean=mu0, sigma=Sigma0)
d0 <- rowSums((X0%*%Sigma0invHalf)^2)
f0 <- density(log10(d0))
xMin <- 1e-1; xMax <- 1e11
file = "khi2density.pdf"
pdf(file,width = 8,height = 6)
plot(f0, xlim=c(log10(xMin), log10(xMax)), ylim=c(0, 2.2),col = klcol[1], main='')  
# curve(dchisq(x, df=d), from=1e-2, to=max(d0), add=TRUE, lty=2)
text(0.98*f0$x[which.max(f0$y)], 1.02*max(f0$y), labels='0', col=klcol[1])
q0_95 <- quantile(d0, probs = .95)
abline(v=log10(quantile(d0, probs=.95)), lty=2)
text(x = log10(quantile(d0, probs=.95))-0.05, y = 2.0, labels = paste0("", round(q0_95, 2)), pos = 4, cex = 0.8)
compt = 1
for(kl in 2:(klNb-2)){
  parms1 <- ParmsF1(m1=m1, k1=k1val[kl], l1=l1valup1[kl], rho1=rho1val[kl])
  X1 <- rmvnorm(B, mean=parms1$mu, sigma=parms1$Sigma)
  d1 <- rowSums((X1%*%Sigma0invHalf)^2)
  f1 <- density(log10(d1))
  lines(f1, col= klcol[compt+1]) 
  #q1_95 <- quantile(d1, probs = .95)
  #abline(v=log10(quantile(d1, probs=.95)), lty=2,col = klcol[(compt+1)])
  # text(
  #   x = log10(quantile(d1, probs = .95)) - 0.05,
  #   y = 2.0 - 0.2 * compt,
  #   labels = paste0(round(q1_95, 2)),  # texte numérique seulement
  #   col = klcol[compt + 1],           # couleur appliquée ici
  #   pos = 4,
  #   cex = 0.8
  # )
  text(0.98*f1$x[which.max(f1$y)], 1.02*max(f1$y), labels= kl-1, col=klcol[compt+1])
compt = compt + 1
  }
dev.off()

# 
# 
# par(mar = c(4, 4, 2, 1))  # marges plus petites : bas, gauche, haut, droite
# 
# plot.new()
# q0 <- quantile(d0, probs = .95)
# 
# abline(v = q0, lty = 2, col = 1)
# 
# # Stockage des quantiles
# q95_vals <- numeric(klNb)
# q95_vals[1] <- q0
# 
# # --- Boucle sur les distributions alternatives
# for (kl in 2:klNb) {
#   parms1 <- ParmsF1(m1 = m1, k1 = k1val[kl], l1 = l1val[kl], rho1 = rho1val[kl])
#   X1 <- rmvnorm(B, mean = parms1$mu, sigma = parms1$Sigma)
#   d1 <- rowSums((X1 %*% Sigma0invHalf)^2)
#   f1 <- density(d1)
#   
#   lines(f1, col = kl)
#   text(0.98 * f1$x[which.max(f1$y)], 1.02 * max(f1$y),
#        labels = KLval[kl], col = kl)
#   
#   q95_vals[kl] <- quantile(d1, probs = .95)
# }
# 
# legend("topright",
#        legend = c("H0", paste0("KL=", KLval[2:klNb])),
#        col = 1:klNb, lty = 1, cex = 0.8)
# 
# # --- Nouveau graphique : quantiles à 95 %
# plot(1:klNb, q95_vals, type = "b", pch = 19, col = "blue",
#      xlab = "Index du modèle (kl)",
#      ylab = "Quantile à 95 %",
#      main = "Quantiles 95% des distributions")
# grid()
# text(1:klNb, q95_vals,
#      labels = round(q95_vals, 2),
#      pos = 3, cex = 0.8, col = "darkblue")