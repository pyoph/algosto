###################################################
#Affiche contamination scenarios
##################################################

rate = 5

######KL = 1##########


k =k1val[2];l=l1valup1[2];rho1 = rho1val[2];



contParam = ParmsF1(m1, k, l, rho1)

data = genererEchantillon(n,n,mu1 = mu0,mu2 = contParam$mu1,Sigma1 = Sigma0,Sigma2 = contParam$Sigma1,rate )

Z = data$Z
# Création d'un dataframe pour ggplot
df <- data.frame(X1 = data$Z[,1], X2 = data$Z[,2], 
                 Label = factor(data$labelsVrais, levels = c(0, 1)))
# Définir les couleurs selon le label
cols <- ifelse(df$Label == 1, "red", "blue")

# Définir les transparences (alpha = 0.5)
# converti avec adjustcolor()
cols <- adjustcolor(cols, alpha.f = 0.5)

pch_vals <- ifelse(df$Label == 1, 4, 19)


file = paste0("contaminScen","-k",k,"-l",l,"-rho1",rho1,"-r",rate,".pdf")

pdf(file, width = 8, height = 6)  # Ouvre un device PDF

# Créer le plot principal
plot(df$X1, df$X2,
     col = cols,
     pch = pch_vals ,
     cex = 1.2,       # taille des points
     xlim = c(a, b),
     ylim = c(a, b),
     xlab = "",
     ylab = "",
     #main = paste("Rate:", rate, "%, k =", k, ", l =", l, ", rho1 =", rho1),
     main = "",
     xaxt = "n", yaxt = "n"  
)

axis(1,cex.axis = 1.8)

axis(2,cex.axis = 1.8)
dev.off()


######KL = 5##########

k =k1val[3];l=l1valup1[3];rho1 = rho1val[3];

contParam = ParmsF1(m1, k, l, rho1)

data = genererEchantillon(n,n,mu1 = mu0,mu2 = contParam$mu1,Sigma1 = Sigma0,Sigma2 = contParam$Sigma1,rate )

Z = data$Z
# Création d'un dataframe pour ggplot
df <- data.frame(X1 = data$Z[,1], X2 = data$Z[,2], 
                 Label = factor(data$labelsVrais, levels = c(0, 1)))
# Définir les couleurs selon le label
cols <- ifelse(df$Label == 1, "red", "blue")

# Définir les transparences (alpha = 0.5)
# converti avec adjustcolor()
cols <- adjustcolor(cols, alpha.f = 0.5)

pch_vals <- ifelse(df$Label == 1, 4, 19)


file = paste0("contaminScen","-k",k,"-l",l,"-rho1",rho1,"-r",rate,".pdf")

pdf(file, width = 8, height = 6)  # Ouvre un device PDF

# Créer le plot principal
plot(df$X1, df$X2,
     col = cols,
     pch = pch_vals,        # points pleins
     cex = 1.2,       # taille des points
     xlim = c(a, b),
     ylim = c(a, b),
     xlab = "",
     ylab = "",
     #main = paste("Rate:", rate, "%, k =", k, ", l =", l, ", rho1 =", rho1),
     main = "",
     xaxt = "n", yaxt = "n"  
)

axis(1,cex.axis = 1.8)

axis(2,cex.axis = 1.8)
dev.off()

#################KL 10####################

k =k1val[4];l=l1valup1[4];rho1 = rho1val[4];

contParam = ParmsF1(m1, k, l, rho1)

data = genererEchantillon(n,n,mu1 = mu0,mu2 = contParam$mu1,Sigma1 = Sigma0,Sigma2 = contParam$Sigma1,rate )

Z = data$Z
# Création d'un dataframe pour ggplot
df <- data.frame(X1 = data$Z[,1], X2 = data$Z[,2], 
                 Label = factor(data$labelsVrais, levels = c(0, 1)))
# Définir les couleurs selon le label
cols <- ifelse(df$Label == 1, "red", "blue")

# Définir les transparences (alpha = 0.5)
# converti avec adjustcolor()
cols <- adjustcolor(cols, alpha.f = 0.5)

pch_vals <- ifelse(df$Label == 1, 4, 19)

file = paste0("contaminScen","-k",k,"-l",l,"-rho1",rho1,"-r",rate,".pdf")

pdf(file, width = 8, height = 6)  # Ouvre un device PDF

# Créer le plot principal
plot(df$X1, df$X2,
     col = cols,
     pch = pch_vals,        # points pleins
     cex = 1.2,       # taille des points
     xlim = c(a, b),
     ylim = c(a, b),
     xlab = "",
     ylab = "",
     #main = paste("Rate:", rate, "%, k =", k, ", l =", l, ", rho1 =", rho1),
     main = "",
     xaxt = "n", yaxt = "n"  
)

axis(1,cex.axis = 1.8)

axis(2,cex.axis = 1.8)
dev.off()


#################KL 25####################

k =k1val[5];l=l1valup1[5];rho1 = rho1val[5];

contParam = ParmsF1(m1, k, l, rho1)

data = genererEchantillon(n,n,mu1 = mu0,mu2 = contParam$mu1,Sigma1 = Sigma0,Sigma2 = contParam$Sigma1,rate )

Z = data$Z
# Création d'un dataframe pour ggplot
df <- data.frame(X1 = data$Z[,1], X2 = data$Z[,2], 
                 Label = factor(data$labelsVrais, levels = c(0, 1)))
# Définir les couleurs selon le label
cols <- ifelse(df$Label == 1, "red", "blue")

# Définir les transparences (alpha = 0.5)
# converti avec adjustcolor()
cols <- adjustcolor(cols, alpha.f = 0.5)

pch_vals <- ifelse(df$Label == 1, 4, 19)

file = paste0("contaminScen","-k",k,"-l",l,"-rho1",rho1,"-r",rate,".pdf")

setwd("~")

pdf(file, width = 8, height = 6)  # Ouvre un device PDF

# Créer le plot principal
plot(df$X1, df$X2,
     col = cols,
     pch = pch_vals,
     cex = 1.2, 
     xlim = c(a, b),
     ylim = c(a, b),
     xlab = "",
     ylab = "",
     #main = paste("Rate:", rate, "%, k =", k, ", l =", l, ", rho1 =", rho1),
     main = "",
     xaxt = "n", yaxt = "n"  
)

axis(1,cex.axis = 1.8)

axis(2,cex.axis = 1.8)
dev.off()


#########################Khi2 densities########################################
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

scen_names = c("d","c","b","a")

# Distribution (log scale)
U <- matrix(rnorm(d*B), B, d)
X0 <- rmvnorm(B, mean=mu0, sigma=Sigma0)
d0 <- rowSums((X0%*%Sigma0invHalf)^2)
f0 <- density(log10(d0))
xMin <- 1e-1; xMax <- 1e6
setwd("~")
file = "khi2density.pdf"
pdf(file,width = 8,height = 6)
plot(f0, xlim=c(log10(xMin), log10(xMax)), ylim=c(0, 2.2),col = klcol[1], main='')  
# curve(dchisq(x, df=d), from=1e-2, to=max(d0), add=TRUE, lty=2)
#text(0.98*f0$x[which.max(f0$y)], 1.02*max(f0$y), labels='0', col=klcol[1])
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
  if(kl != 1){
    text(0.98*f1$x[which.max(f1$y)], 1.02*max(f1$y), labels= scen_names[kl-1], col=klcol[compt+1])}
  compt = compt + 1
}
dev.off()


#######################################################False positives and false negatives and Sigma error final########################################

# --- Transformation pseudo-log  ---
pseudo_log <- function(y) {
  return(log10(1+y))  
}


#######Near scenario###############

k = k1val[2]; l = l1valup1[2]; rho1 = rho1val[2]

##Sigma error

setwd("~")

file <- paste0("erreurs_Sigma_final-k", k, "-l", l, "-rho1", rho1, ".pdf")

pdf(file, width = 8, height = 6)

# Plot principal — axe Y en log (automatique)
plot(rList[1:9], moyenne_erreursSigmaNear[n,1:9 , 3],
     type = "l", lwd = 4, col = "red",
     xlab = "", ylab = "",
     yaxt = "n", xaxt = "n",
     log = "y",                       # ✅ active l’échelle log pour Y
     ylim = c(1e-1, 1e2)              # borne Y cohérente
)

# Autres courbes
lines(rList[1:9], moyenne_erreursSigmaNear[n,1:9 , 3], lwd = 4, col = "blue", lty = "dashed")
lines(rList[1:9], moyenne_erreursSigmaNear[n,1:9 , 3], lwd = 4, col = "darkgreen", lty = "dotted")

# Axe Y logarithmique lisible
log_ticks <- 10^seq(-1, 2, by = 1)
axis(2, at = log_ticks,
     labels = parse(text = paste0("10^", -1:2)),
     las = 1, cex.axis = 2.1)

# Axe X (je suppose que tu veux une échelle régulière)
x_ticks <- rList
axis(1, at = x_ticks, las = 1, cex.axis = 2.3)

box()  # bordure
dev.off()



#False negatives

file = paste0("false_negatives_final-k", k, "-l", l, "-rho1", rho1, ".pdf")
setwd("~")
pdf(file, width = 8, height = 6)
# --- Données transformées ---
epsilon <- 0
y_red   <- moyenne_faux_negatifsNear[2:9,3]/((rList[2:9]/100)*n)*100
y_blue  <- moyenne_faux_negatifsNear[2:9,2]/((rList[2:9]/100)*n)*100
y_green <- moyenne_faux_negatifsNear[2:9,1]/((rList[2:9]/100)*n)*100
y_purp  <- moyenne_faux_negatifsOracle[2:9,1]/((rList[2:9]/100)*n)*100


# --- Échelle des ticks ---
yticks <- c(0, 1, 10, 100)

# --- Plot principal ---
plot(rList[2:9], pseudo_log(y_red),
     type = "l", lwd = 4, col = "red",
     xlab = "", ylab = "", yaxt = "n", xaxt = "n",
     ylim = pseudo_log(c(0, 100))
)
lines(rList[2:9], pseudo_log(y_blue),  lwd = 4, col = "blue",      lty = "dashed")
lines(rList[2:9], pseudo_log(y_green), lwd = 4, col = "darkgreen", lty = "dotted")
lines(rList[2:9], pseudo_log(y_purp),  lwd = 4, col = "purple4",    lty = "longdash")

# --- Axes ---
axis(1, at = rList[2:9], las = 1, cex.axis = 2.3)
axis(2, at = pseudo_log(yticks),
     labels = c("0", expression(10^0), expression(10^1), expression(10^2)),
     las = 1, cex.axis = 2.3)
box()

dev.off()

#False positives Near

file = paste0("false_positives_final-k",k,"-l",l,"-rho1",rho1,".pdf")
setwd("~")
pdf(file,width = 8, height = 6) 


plot(rList[1:9], faux_positifsNear[1:9,3,1]/((1 - rList[1:9]/100)*n)*100,
     type = "l", lwd = 4,col = "red", 
     xlab = "", ylab = "",   # Pas de label
     yaxt = "n", xaxt = "n", # On masque les axes par défaut
     #log = "y",              # Échelle logarithmique Y
     ylim = c(0, 20)     # Plage Y adaptée à tes ticks log
)
# 
# # Autres courbes
lines(rList[1:9], faux_positifsNear[1:9,2,1]/((1 - rList[1:9]/100)*n)*100,lwd = 4,col = "blue", lty = "dashed")
lines(rList[1:9], faux_positifsNear[1:9,1,1]/((1 - rList[1:9]/100)*n)*100,lwd = 4, col = "darkgreen", lty = "dotted")
lines(rList[1:9], faux_positifsOracle[1:9,1]/((1 - rList[1:9]/100)*n)*100,lwd = 4, col = "purple", lty = "longdash")


x_ticks <- rList[2:length(rList)]
axis(2, at = seq(0, 100, by = 5), las = 1,cex.axis = 2.3)

axis(1, at = x_ticks,las = 1,cex.axis = 2.3)
dev.off()


k = k1val[3]; l = l1valup1[3]; rho1 = rho1val[3]

##Sigma error

setwd("~")

file <- paste0("erreurs_Sigma_final-k", k, "-l", l, "-rho1", rho1, ".pdf")

pdf(file, width = 8, height = 6)


plot(rList[1:9], moyenne_erreursSigmaMed[n, 1:9, 3],
     type = "l", lwd = 4, col = "red",
     xlab = "", ylab = "",
     yaxt = "n", xaxt = "n",
     log = "y",                       # ✅ active l’échelle log pour Y
     ylim = c(1e-1, 1e2)              # borne Y cohérente
)

# Autres courbes
lines(rList[1:9], moyenne_erreursSigmaMed[n, 1:9, 2], lwd = 4, col = "blue", lty = "dashed")
lines(rList[1:9], moyenne_erreursSigmaMed[n,1:9 ,1], lwd = 4, col = "darkgreen", lty = "dotted")

# Axe Y 
log_ticks <- 10^seq(-1, 2, by = 1)
axis(2, at = log_ticks,
     labels = parse(text = paste0("10^", -1:2)),
     las = 1, cex.axis = 2.1)

# Axe X 
x_ticks <- rList
axis(1, at = x_ticks, las = 1, cex.axis = 2.3)

box()  # bordure
dev.off()



#False negatives
file = paste0("false_negatives_final-k", k, "-l", l, "-rho1", rho1, ".pdf")
setwd("~")
pdf(file, width = 8, height = 6)

y_red   <- moyenne_faux_negatifsMed[2:9,3]/((rList[2:9]/100)*n)*100
y_blue  <- moyenne_faux_negatifsMed[2:9,2]/((rList[2:9]/100)*n)*100
y_green <- moyenne_faux_negatifsMed[2:9,1]/((rList[2:9]/100)*n)*100
y_purp  <- moyenne_faux_negatifsMedOracle[2:9,1]/((rList[2:9]/100)*n)*100


# --- Échelle des ticks ---
yticks <- c(0, 1, 10, 100)

# --- Plot principal ---
plot(rList[2:9], pseudo_log(y_red),
     type = "l", lwd = 4, col = "red",
     xlab = "", ylab = "", yaxt = "n", xaxt = "n",
     ylim = pseudo_log(c(0, 100))
)
lines(rList[2:9], pseudo_log(y_blue),  lwd = 4, col = "blue",      lty = "dashed")
lines(rList[2:9], pseudo_log(y_green), lwd = 4, col = "darkgreen", lty = "dotted")
lines(rList[2:9], pseudo_log(y_purp),  lwd = 4, col = "purple4",    lty = "longdash")

# --- Axes ---
axis(1, at = rList[2:9], las = 1, cex.axis = 2.3)
axis(2, at = pseudo_log(yticks),
     labels = c("0", expression(10^0), expression(10^1), expression(10^2)),
     las = 1, cex.axis = 2.3)
box()

dev.off()

#False positives Med

file = paste0("false_positives_final-k",k,"-l",l,"-rho1",rho1,".pdf")
setwd("~")
pdf(file,width = 8, height = 6) 


plot(rList[1:9], faux_positifsMed[1:9,3,1]/((1 - rList[1:9]/100)*n)*100,
     type = "l", lwd = 4,col = "red", 
     xlab = "", ylab = "",   # Pas de label
     yaxt = "n", xaxt = "n", # On masque les axes par défaut
     #log = "y",              # Échelle logarithmique Y
     ylim = c(0, 20)     # Plage Y adaptée à tes ticks log
)
# 
# # Autres courbes
lines(rList[1:9], faux_positifsMed[1:9,2,1]/((1 - rList[1:9]/100)*n)*100,lwd = 4,col = "blue", lty = "dashed")
lines(rList[1:9], faux_positifsMed[1:9,1,1]/((1 - rList[1:9]/100)*n)*100,lwd = 4, col = "darkgreen", lty = "dotted")
lines(rList[1:9], faux_positifsOracleMed[1:9,1]/((1 - rList[1:9]/100)*n)*100,lwd = 4, col = "purple", lty = "longdash")

#  
x_ticks <- rList[1:9]
axis(2, at = seq(0, 100, by = 5), las = 1,cex.axis = 2.1)

axis(1, at = x_ticks,las = 1,cex.axis = 2.3)
dev.off()

k = k1val[4]; l = l1valup1[4]; rho1 = rho1val[4]


file <- paste0("erreurs_Sigma_final-k", k, "-l", l, "-rho1", rho1, ".pdf")

pdf(file, width = 8, height = 6)

# Sigma estimation error

plot(rList[1:9], moyenne_erreursSigmaMed2[n, 1:9, 3],
     type = "l", lwd = 4, col = "red",
     xlab = "", ylab = "",
     yaxt = "n", xaxt = "n",
     log = "y",                       # ✅ active l’échelle log pour Y
     ylim = c(1e-1, 1e2)              # borne Y cohérente
)

# Autres courbes
lines(rList[1:9], moyenne_erreursSigmaMed2[n, 1:9, 2], lwd = 4, col = "blue", lty = "dashed")
lines(rList[1:9],moyenne_erreursSigmaMed2[n,1:9 , 1], lwd = 4, col = "darkgreen", lty = "dotted")

# Axe Y logarithmique lisible
log_ticks <- 10^seq(-1, 2, by = 1)
axis(2, at = log_ticks,
     labels = parse(text = paste0("10^", -1:2)),
     las = 1, cex.axis = 2.1)

# Axe X 
x_ticks <- rList
axis(1, at = x_ticks, las = 1, cex.axis = 2.3)

box()  # bordure
dev.off()


#False negatives

file = paste0("false_negatives_final-k", k, "-l", l, "-rho1", rho1, ".pdf")
setwd("~")
pdf(file, width = 8, height = 6)
# --- Données transformées ---

y_red   <- moyenne_faux_negatifsMed2[2:9,3]/((rList[2:9]/100)*n)*100
y_blue  <- moyenne_faux_negatifsMed2[2:9,2]/((rList[2:9]/100)*n)*100
y_green <- moyenne_faux_negatifsMed2[2:9,1]/((rList[2:9]/100)*n)*100
y_purp  <- moyenne_faux_negatifsMed2Oracle[2:9,1]/((rList[2:9]/100)*n)*100


# --- Échelle des ticks ---
yticks <- c(0, 1, 10, 100)

# --- Plot principal ---
plot(rList[2:9], pseudo_log(y_red),
     type = "l", lwd = 4, col = "red",
     xlab = "", ylab = "", yaxt = "n", xaxt = "n",
     ylim = pseudo_log(c(0, 100))
)
lines(rList[2:9], pseudo_log(y_blue),  lwd = 4, col = "blue",      lty = "dashed")
lines(rList[2:9], pseudo_log(y_green), lwd = 4, col = "darkgreen", lty = "dotted")
lines(rList[2:9], pseudo_log(y_purp),  lwd = 4, col = "purple4",    lty = "longdash")

# --- Axes ---
axis(1, at = rList[2:9], las = 1, cex.axis = 2.3)
axis(2, at = pseudo_log(yticks),
     labels = c("0", expression(10^0), expression(10^1), expression(10^2)),
     las = 1, cex.axis = 2.3)
box()

dev.off()





file = paste0("false_positives_final-k",k,"-l",l,"-rho1",rho1,".pdf")
setwd("~")
pdf(file,width = 8, height = 6) 


plot(rList[1:9], moyenne_faux_positifsMed2[1:9,3]/((1 - rList[1:9]/100)*n)*100,
     type = "l", lwd = 4,col = "red", 
     xlab = "", ylab = "",   # Pas de label
     yaxt = "n", xaxt = "n", # On masque les axes par défaut
     #log = "y",              # Échelle logarithmique Y
     ylim = c(0, 20)     # Plage Y adaptée à tes ticks log
)
# 
# # Autres courbes
lines(rList[1:9], moyenne_faux_positifsMed2[1:9,2]/((1 - rList[1:9]/100)*n)*100,lwd = 4,col = "blue", lty = "dashed")
lines(rList[1:9], moyenne_faux_positifsMed2[1:9,1]/((1 - rList[1:9]/100)*n)*100,lwd = 4, col = "darkgreen", lty = "dotted")
lines(rList[1:9], moyenne_faux_positifsMed2Oracle[1:9,1]/((1 - rList[1:9]/100)*n)*100,lwd = 4, col = "purple", lty = "longdash")

#  
# # Axe Y logarithmique
# log_ticks <- 10^seq(-1, 10, by = 1)
# axis(2, at = log_ticks,
#      labels = parse(text = paste0("10^", -1:10)),
#      las = 1)
x_ticks <- rList[1:9]
axis(2, at = seq(0, 100, by = 5), las = 1,cex.axis = 2.3)

axis(1, at = x_ticks,las = 1,cex.axis = 2.3)
dev.off()

k = k1val[5]; l = l1valup1[5]; rho1 = rho1val[5]

#Sigma error 

setwd("~")

file <- paste0("erreurs_Sigma_final-k", k, "-l", l, "-rho1", rho1, ".pdf")

pdf(file, width = 8, height = 6)

# Plot principal — axe Y en log (automatique)
plot(rList[1:9], moyenne_erreursSigmaMed3[n, 1:9, 3],
     type = "l", lwd = 4, col = "red",
     xlab = "", ylab = "",
     yaxt = "n", xaxt = "n",
     log = "y",                       # ✅ active l’échelle log pour Y
     ylim = c(1e-1, 1e2)              # borne Y cohérente
)

# Autres courbes
lines(rList[1:9], moyenne_erreursSigmaMed3[n,1:9 , 2], lwd = 4, col = "blue", lty = "dashed")
lines(rList[1:9], moyenne_erreursSigmaMed3[n,1:9 , 1], lwd = 4, col = "darkgreen", lty = "dotted")

# Axe Y logarithmique lisible
log_ticks <- 10^seq(-1, 2, by = 1)
axis(2, at = log_ticks,
     labels = parse(text = paste0("10^", -1:2)),
     las = 1, cex.axis = 2.1)

# Axe X (je suppose que tu veux une échelle régulière)
x_ticks <- rList[1:9]
axis(1, at = x_ticks, las = 1, cex.axis = 2.3)

box()  # bordure
dev.off()




#False negatives

file = paste0("false_negatives_final-k", k, "-l", l, "-rho1", rho1, ".pdf")
setwd("~/figures")
pdf(file, width = 8, height = 6)
# --- Données transformées ---
y_red   <- moyenne_faux_negatifsMed3[2:9,3]/((rList[2:9]/100)*n)*100
y_blue  <- moyenne_faux_negatifsMed3[2:9,2]/((rList[2:9]/100)*n)*100
y_green <- moyenne_faux_negatifsMed3[2:9,1]/((rList[2:9]/100)*n)*100
y_purp  <- moyenne_faux_negatifsMed3Oracle[2:9,1]/((rList[2:9]/100)*n)*100

# --- Échelle des ticks ---
yticks <- c(0, 1, 10, 100)

# --- Plot principal ---
plot(rList[2:9], pseudo_log(y_red),
     type = "l", lwd = 4, col = "red",
     xlab = "", ylab = "", yaxt = "n", xaxt = "n",
     ylim = pseudo_log(c(0, 100))
)
lines(rList[2:9], pseudo_log(y_blue),  lwd = 4, col = "blue",      lty = "dashed")
lines(rList[2:9], pseudo_log(y_green), lwd = 4, col = "darkgreen", lty = "dotted")
lines(rList[2:9], pseudo_log(y_purp),  lwd = 4, col = "purple4",    lty = "longdash")

# --- Axes ---
axis(1, at = rList[2:9], las = 1, cex.axis = 2.3)
axis(2, at = pseudo_log(yticks),
     labels = c("0", expression(10^0), expression(10^1), expression(10^2)),
     las = 1, cex.axis = 2.3)
box()

dev.off()

file = paste0("false_positives_final-k",k,"-l",l,"-rho1",rho1,".pdf")
setwd("~")
pdf(file,width = 8, height = 6) 


plot(rList[1:9], moyenne_faux_positifsMed3[1:9,3]/((1 - rList[1:9]/100)*n)*100,
     type = "l", lwd = 4,col = "red", 
     xlab = "", ylab = "",   # Pas de label
     yaxt = "n", xaxt = "n", # On masque les axes par défaut
     #log = "y",              # Échelle logarithmique Y
     ylim = c(0, 20)     # Plage Y adaptée à tes ticks log
)
# 
# # Autres courbes
lines(rList[1:9], moyenne_faux_positifsMed3[1:9,2]/((1 - rList[1:9]/100)*n)*100,lwd = 4,col = "blue", lty = "dashed")
lines(rList[1:9], moyenne_faux_positifsMed3[1:9,1]/((1 - rList[1:9]/100)*n)*100,lwd = 4, col = "darkgreen", lty = "dotted")
lines(rList[1:9], moyenne_faux_positifsMed3Oracle[1:9,1]/((1 - rList[1:9]/100)*n)*100,lwd = 4, col = "purple", lty = "longdash")

#  
x_ticks <- rList[1:9]
axis(2, at = seq(0, 100, by = 5), las = 1,cex.axis = 2.3)

axis(1, at = x_ticks,las = 1,cex.axis = 2.3)
dev.off()




########################Trajectoires####################
x_vals = 1:n

ind_rates = c(2,5,7,9)




for(r in ind_rates){
  rate = rList[r]
  file = paste0("erreurSigma-traj-k",k,"-l",l,"-rho1",rho1,"-rate",rate,".pdf")
  setwd("~")
  pdf(file,width = 8, height = 6) 
  
  
  plot(x_vals,moyenne_erreursSigmaMed3[,r,3],     type = "l", lwd = 4, col= "red",
       xlab = "", ylab = "",   # Pas de label
       yaxt = "n", xaxt = "n", # On masque les axes par défaut
       log = "y",              # Échelle logarithmique Y
       ylim = c(1e-1, 1e2)     # Plage Y adaptée à tes ticks log
  )
  lines(x_vals,moyenne_erreursSigmaMed3[,r,1],lty = "dotted",col = "darkgreen",lwd = 4)
  
  lines(x_vals,moyenne_erreursSigmaMed3[,r,2],lty = "dashed",col = "blue",lwd = 4)
  
  #lines(x_vals,taux_fn_cum[,r,4]*100,lty = "dotted",col = "orange",lwd = 4)
  # Axe Y logarithmique lisible
  log_ticks <- 10^seq(-1, 2, by = 1)
  axis(2, at = log_ticks,
       labels = parse(text = paste0("10^", -1:2)),
       las = 1, cex.axis = 2.3)
  
  axis(1, at = seq(1000, max(x_vals), by = 1000), las = 1, cex.axis =2.3)
  dev.off()
}


for(r in ind_rates){
  rate = rList[r]
  file = paste0("false_negatives-traj-k",k,"-l",l,"-rho1",rho1,"-rate",rate,".pdf")
  setwd("~")
  pdf(file,width = 8, height = 6) 
  
  # Transformation log10(1 + y)
  y_red   <- log10(1 + taux_fn_cum[, r, 3] * 100)
  y_green <- log10(1 + taux_fn_cum[, r, 1] * 100)
  y_blue  <- log10(1 + taux_fn_cum[, r, 2] * 100)
  y_purp  <- log10(1 + taux_fn_cum[, r, 4] * 100)
  
  # Tracé principal
  plot(x_vals, y_red,
       type = "l", lwd = 4, col = "red",
       xlab = "", ylab = "",
       yaxt = "n", xaxt = "n",
       ylim = log10(1 + c(0, 100))
  )
  lines(x_vals, y_green, lty = "dotted",  col = "darkgreen", lwd = 4)
  lines(x_vals, y_blue,  lty = "dashed",  col = "blue", lwd = 4)
  lines(x_vals, y_purp,  lty = "longdash", col = "purple4", lwd = 4)
  
  # Axe Y : puissances de 10 avec 0 au bon endroit
  yticks <- c(0, 1, 10, 100)
  axis(2, at = log10(1 + yticks),
       labels = c("0", expression(10^0), expression(10^1), expression(10^2)),
       las = 1, cex.axis = 1.8)
  
  # Axe X
  axis(1, at = seq(1000, max(x_vals), by = 1000), las = 1, cex.axis = 1.8)
  
  box()
  dev.off()
}


for(r in ind_rates){
  rate = rList[r]
  file = paste0("false_positives-traj-k",k,"-l",l,"-rho1",rho1,"-rate",rate,".pdf")
  setwd("~")
  pdf(file,width = 8, height = 6) 
  
  
  plot(x_vals,taux_fp_cum[,r,3]*100,     type = "l", lwd = 4, col= "red",
       xlab = "", ylab = "",   # Pas de label
       yaxt = "n", xaxt = "n", # On masque les axes par défaut
       #log = "y",              # Échelle logarithmique Y
       ylim = c(0, 10)     # Plage Y adaptée à tes ticks log
  )
  lines(x_vals,taux_fp_cum[,r,1]*100,lty = "dotted",col = "darkgreen",lwd = 4)
  
  lines(x_vals,taux_fp_cum[,r,2]*100,lty = "dashed",col = "blue",lwd = 4)
  
  lines(x_vals, taux_fp_cum[,r,4]*100, lty = "longdash", col = "purple4", lwd = 4)
  axis(2, at = seq(0, 10, by = 5), las = 1,cex.axis = 1.8)
  
  axis(1, at = seq(1000, max(x_vals), by = 1000), las = 1, cex.axis = 1.8)
  dev.off()
}



# ##############################################
# Temps calculs
# ################# #############################
setwd("~")
file = paste("boxplotn", n, "d", d, ".pdf", sep = "")

pdf(file = file, width = 8, height = 6)

# tracer le boxplot sans axe y
boxplot(temps[2,], temps[3,],
        col = c("blue", "red"),
        lwd = 2,
        las = 1,
        ylab = "",
        xlab = "",
        cex.axis = 2.2,
        cex.lab = 2.2,
        log = "y",                      # axe Y en log
        ylim = c(1e-1, 100),               # de 1 à 100
        names = c("online", "streaming"),
        axes = FALSE)                   # désactive les axes par défaut

# tracer l'axe X
axis(1, at = 1:2, labels = c("online", "streaming"), cex.axis = 2.2)

# Axe Y logarithmique (0.1, 1, 10, 100)
y_ticks <- c(0.1, 1, 10, 100)
axis(2, at = y_ticks,
     labels = parse(text = c("10^-1", "10^0", "10^1", "10^2")),
     las = 1, cex.axis = 2.1)


# cadre autour du graphe
box()

dev.off()


