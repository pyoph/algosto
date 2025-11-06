###################################################
#Affiche contamination scenarios
##################################################



afficherContaminationScenarios = function(k,l,rho1,contamination,rate,a,b){
  contParam = ParmsF1(m1, k, l, rho1)
  
    data = genererEchantillon(n,n,mu1 = mu0,mu2 = contParam$mu1,Sigma1 = Sigma0,Sigma2 = contParam$Sigma1,rate )
    
 Z = data$Z
 # Création d'un dataframe pour ggplot
 df <- data.frame(X1 = data$Z[,1], X2 = data$Z[,2], 
                  Label = factor(data$labelsVrais, levels = c(0, 1)))
 #Création du plot
 # Create the plot (English version)
 p <- ggplot(df, aes(x = X1, y = X2, color = Label, alpha = Label)) +
   geom_point(size = 3) +
   ylim(a,b) +
   xlim(a,b) +
   scale_color_manual(
     values = c("0" = "blue", "1" = "red"),
     labels = c("Inlier", "Outlier"),
     name = "Status"
   ) +
   scale_alpha_manual(
     values = c("0" = 0.5, "1" = 0.5),  
     guide = "none"  # Hide alpha from legend
   ) +
   labs(
     title = "",
     subtitle = paste("Rate:", rate, "%, k =", k, ", l =", l, ", rho1 =", rho1),
     x = "",
     y = ""
   ) +
   theme_minimal() +
   theme(
     legend.position = "bottom",
     plot.title = element_text(face = "bold")
   )
  return(p)
}

a = -5
b = 5

k =k1val[1];l=l1val[1];rho1 = rho1val[1];
rate  = 5
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

file = paste0("contaminScen_no_outlierr",rate,".pdf")

pdf(file, width = 8, height = 6)  # Ouvre un device PDF

# Créer le plot principal
plot(df$X1, df$X2,
     col = cols,
     pch = 19,        # points pleins
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

file = paste0("contaminScen","-k",k,"-l",l,"-rho1",rho1,"-r",rate,".pdf")

pdf(file, width = 8, height = 6)  # Ouvre un device PDF

# Créer le plot principal
plot(df$X1, df$X2,
     col = cols,
     pch = 19,        # points pleins
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

file = paste0("contaminScen","-k",k,"-l",l,"-rho1",rho1,"-r",rate,".pdf")

pdf(file, width = 8, height = 6)  # Ouvre un device PDF

# Créer le plot principal
plot(df$X1, df$X2,
     col = cols,
     pch = 19,        # points pleins
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

file = paste0("contaminScen","-k",k,"-l",l,"-rho1",rho1,"-r",rate,".pdf")

pdf(file, width = 8, height = 6)  # Ouvre un device PDF

# Créer le plot principal
plot(df$X1, df$X2,
     col = cols,
     pch = 19,        # points pleins
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

file = paste0("contaminScen","-k",k,"-l",l,"-rho1",rho1,"-r",rate,".pdf")

setwd("~")

pdf(file, width = 8, height = 6)  # Ouvre un device PDF

# Créer le plot principal
plot(df$X1, df$X2,
     col = cols,
     pch = 19,        # points pleins
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



#################KL 50####################

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

file = paste0("contaminScen","-k",k,"-l",l,"-rho1",rho1,"-r",rate,".pdf")

pdf(file, width = 8, height = 6)  # Ouvre un device PDF

# Créer le plot principal
plot(df$X1, df$X2,
     col = cols,
     pch = 19,        # points pleins
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


#################KL 100####################

k =k1val[6];l=l1valup1[6];rho1 = rho1val[6];

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

file = paste0("contaminScen","-k",k,"-l",l,"-rho1",rho1,"-r",rate,".pdf")

pdf(file, width = 8, height = 6)  # Ouvre un device PDF

# Créer le plot principal
plot(df$X1, df$X2,
     col = cols,
     pch = 19,        # points pleins
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


###########################Densités Khi2 pour différentes KL########################################

##############Densité Khi2 pour KL 1, 5,10,25################################################################
# 
# 
# # Domaine
# x_vals <- seq(0, 40, length.out = 1000)
# 
# # Densités
# dens_centrale <- dchisq(x_vals, df = d)
# ################KL 1#####################"
# 
# contParam = ParmsF1(m1, k1val[2], l1valup1[2], rho1val[2])
# mu_1 = contParam$mu1
# Sigma1 = contParam$Sigma1
# lambdaKL1 <- as.numeric(t(mu_1 - mu0) %*% solve(Sigma0) %*% (mu_1 - mu0))
# dens_noncentraleKL1 <- dchisq(x_vals, df = d, ncp = lambdaKL1)
# ############################KL 5#########################
# contParam = ParmsF1(m1, k1val[3], l1valup1[3], rho1val[3])
# mu_1 = contParam$mu1
# Sigma1 = contParam$Sigma1
# lambdaKL5 <- as.numeric(t(mu_1 - mu0) %*% solve(Sigma1) %*% (mu_1 - mu0))
# dens_noncentraleKL5 <- dchisq(x_vals, df = d, ncp = lambdaKL5)
# 
# ###########################KL 10###########################
# contParam = ParmsF1(m1, k1val[4], l1valup1[4], rho1val[4])
# mu_1 = contParam$mu1
# Sigma1 = contParam$Sigma1
# lambdaKL10 <- as.numeric(t(mu_1 - mu0) %*% solve(Sigma1) %*% (mu_1 - mu0))
# dens_noncentraleKL10 <- dchisq(x_vals, df = d, ncp = lambdaKL10)
# 
# 
# ###########################KL 25###########################
# contParam = ParmsF1(m1, k1val[5], l1valup1[5], rho1val[5])
# mu_1 = contParam$mu1
# Sigma1 = contParam$Sigma1
# lambdaKL25 <- as.numeric(t(mu_1 - mu0) %*% solve(Sigma1) %*% (mu_1 - mu0))
# dens_noncentraleKL25 <- dchisq(x_vals, df = d, ncp = lambdaKL25)
# 
# ##########################KL 50###############################
# contParam = ParmsF1(m1, k1val[6], l1valup1[6], rho1val[6])
# mu_1 = contParam$mu1
# Sigma1 = contParam$Sigma1
# lambdaKL50 <- as.numeric(t(mu_1 - mu0) %*% solve(Sigma1) %*% (mu_1 - mu0))
# dens_noncentraleKL50 <- dchisq(x_vals, df = d, ncp = lambdaKL50)
# ##########################KL 100###############################
# contParam = ParmsF1(m1, k1val[7], l1valup1[7], rho1val[7])
# mu_1 = contParam$mu1
# Sigma1 = contParam$Sigma1
# lambdaKL100 <- as.numeric(t(mu_1 - mu0) %*% solve(Sigma1) %*% (mu_1 - mu0))
# dens_noncentraleKL100 <- dchisq(x_vals, df = d, ncp = lambdaKL100)
# 
# setwd("~")
# 
# file = paste0("khi2density",".pdf")
# 
# pdf(file, width = 8, height = 6)  # Ouvre un device PDF
# 
# # Tracé sans titre ni légende
# plot(x_vals, dens_centrale, type = "l", lwd = 4, col = "blue",
#      ylim = c(0, max(dens_centrale)),
#      xlab = "", ylab = "", axes = FALSE)
# lines(x_vals, dens_noncentraleKL1,  col = "darkgreen", lwd = 4)
# lines(x_vals, dens_noncentraleKL5,  col = "orange",    lwd = 4)
# lines(x_vals, dens_noncentraleKL10, col = "red",       lwd = 4)
# lines(x_vals, dens_noncentraleKL25, col = "grey",   lwd = 4)
# lines(x_vals, dens_noncentraleKL50, col = "darkred",   lwd = 4)
# lines(x_vals, dens_noncentraleKL100, col = "black",   lwd = 4)
# 
# axis(1)
# axis(2)
# 
# dev.off()
# ###########################

#######################################################False positives and false negatives and Sigma error final########################################
k = k1val[2]; l = l1valup1[2]; rho1 = rho1val[2]


file = paste0("false_negatives_final-k", k, "-l", l, "-rho1", rho1, ".pdf")
setwd("~")
pdf(file, width = 8, height = 6)
# --- Données transformées ---
epsilon <- 0
y_red   <- moyenne_faux_negatifsMed[2:9,3]/((rList[2:9]/100)*n)*100
y_blue  <- moyenne_faux_negatifsMed[2:9,2]/((rList[2:9]/100)*n)*100
y_green <- moyenne_faux_negatifsMed[2:9,1]/((rList[2:9]/100)*n)*100
y_purp  <- moyenne_faux_negatifsMedOracle[2:9,1]/((rList[2:9]/100)*n)*100

# --- Transformation pseudo-log (linéaire entre 0 et 1, log10 au-delà) ---
pseudo_log <- function(y) {
  ifelse(y <= 1, y, 1 + log10(y))  # 0–1 linéaire, puis log
}

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
lines(rList[2:9], pseudo_log(y_purp),  lwd = 4, col = "purple",    lty = "longdash")

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
lines(rList[1:9], faux_positifsOracle[1:9,1]/((1 - rList[1:9]/100)*n)*100,lwd = 4, col = "orange", lty = "dotted")

#  
# # Axe Y logarithmique
# log_ticks <- 10^seq(-1, 10, by = 1)
# axis(2, at = log_ticks,
#      labels = parse(text = paste0("10^", -1:10)),
#      las = 1)
x_ticks <- rList[2:length(rList)]
axis(2, at = seq(0, 100, by = 5), las = 1,cex.axis = 2.3)

axis(1, at = x_ticks,las = 1,cex.axis = 2.3)
dev.off()



setwd("~")

file <- paste0("erreurs_Sigma_final-k", k, "-l", l, "-rho1", rho1, ".pdf")

pdf(file, width = 8, height = 6)

# Plot principal — axe Y en log (automatique)
plot(rList[1:9], erreursSigmaNear[n,1:9 , 3, 1],
     type = "l", lwd = 4, col = "red",
     xlab = "", ylab = "",
     yaxt = "n", xaxt = "n",
     log = "y",                       # ✅ active l’échelle log pour Y
     ylim = c(1e-1, 1e2)              # borne Y cohérente
)

# Autres courbes
lines(rList[1:9], erreursSigmaNear[n, 1:9, 2, 1], lwd = 4, col = "blue", lty = "dashed")
lines(rList[1:9], erreursSigmaNear[n, 1:9, 1, 1], lwd = 4, col = "darkgreen", lty = "dotted")

# Axe Y logarithmique lisible
log_ticks <- 10^seq(-1, 2, by = 1)
axis(2, at = log_ticks,
     labels = parse(text = paste0("10^", -1:2)),
     las = 1, cex.axis = 1.8)

# Axe X (je suppose que tu veux une échelle régulière)
x_ticks <- rList
axis(1, at = x_ticks, las = 1, cex.axis = 2.3)

box()  # bordure
dev.off()


k = k1val[3]; l = l1valup1[3]; rho1 = rho1val[3]


file = paste0("false_negatives_final-k", k, "-l", l, "-rho1", rho1, ".pdf")
setwd("~")
pdf(file, width = 8, height = 6)
# --- Données transformées ---
epsilon <- 0
y_red   <- moyenne_faux_negatifsMed[2:9,3]/((rList[2:9]/100)*n)*100
y_blue  <- moyenne_faux_negatifsMed[2:9,2]/((rList[2:9]/100)*n)*100
y_green <- moyenne_faux_negatifsMed[2:9,1]/((rList[2:9]/100)*n)*100
y_purp  <- moyenne_faux_negatifsMedOracle[2:9,1]/((rList[2:9]/100)*n)*100

# --- Transformation pseudo-log (linéaire entre 0 et 1, log10 au-delà) ---
pseudo_log <- function(y) {
  ifelse(y <= 1, y, 1 + log10(y))  # 0–1 linéaire, puis log
}

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
lines(rList[2:9], pseudo_log(y_purp),  lwd = 4, col = "purple",    lty = "longdash")

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
lines(rList[1:9], faux_positifsOracleMed[1:9,1]/((1 - rList[1:9]/100)*n)*100,lwd = 4, col = "orange", lty = "dotted")

#  
# # Axe Y logarithmique
# log_ticks <- 10^seq(-1, 10, by = 1)
# axis(2, at = log_ticks,
#      labels = parse(text = paste0("10^", -1:10)),
#      las = 1)
x_ticks <- rList[1:9]
axis(2, at = seq(0, 100, by = 5), las = 1,cex.axis = 2.1)

axis(1, at = x_ticks,las = 1,cex.axis = 2.3)
dev.off()


setwd("~")

file <- paste0("erreurs_Sigma_final-k", k, "-l", l, "-rho1", rho1, ".pdf")

pdf(file, width = 8, height = 6)

# Plot principal — axe Y en log (automatique)
plot(rList[1:9], erreursSigmaMed[n, 1:9, 3, 1],
     type = "l", lwd = 4, col = "red",
     xlab = "", ylab = "",
     yaxt = "n", xaxt = "n",
     log = "y",                       # ✅ active l’échelle log pour Y
     ylim = c(1e-1, 1e2)              # borne Y cohérente
)

# Autres courbes
lines(rList[1:9], erreursSigmaMed[n, 1:9, 2, 1], lwd = 4, col = "blue", lty = "dashed")
lines(rList[1:9], erreursSigmaMed[n,1:9 , 1, 1], lwd = 4, col = "darkgreen", lty = "dotted")

# Axe Y logarithmique lisible
log_ticks <- 10^seq(-1, 2, by = 1)
axis(2, at = log_ticks,
     labels = parse(text = paste0("10^", -1:2)),
     las = 1, cex.axis = 1.8)

# Axe X (je suppose que tu veux une échelle régulière)
x_ticks <- rList
axis(1, at = x_ticks, las = 1, cex.axis = 2.3)

box()  # bordure
dev.off()

k = k1val[4]; l = l1valup1[4]; rho1 = rho1val[4]

file = paste0("false_negatives_final-k", k, "-l", l, "-rho1", rho1, ".pdf")
setwd("~")
pdf(file, width = 8, height = 6)
# --- Données transformées ---
epsilon <- 0
y_red   <- moyenne_faux_negatifsMed2[2:9,3]/((rList[2:9]/100)*n)*100
y_blue  <- moyenne_faux_negatifsMed2[2:9,2]/((rList[2:9]/100)*n)*100
y_green <- moyenne_faux_negatifsMed2[2:9,1]/((rList[2:9]/100)*n)*100
y_purp  <- moyenne_faux_negatifsMed2Oracle[2:9,1]/((rList[2:9]/100)*n)*100

# --- Transformation pseudo-log (linéaire entre 0 et 1, log10 au-delà) ---
pseudo_log <- function(y) {
  ifelse(y <= 1, y, 1 + log10(y))  # 0–1 linéaire, puis log
}

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
lines(rList[2:9], pseudo_log(y_purp),  lwd = 4, col = "purple",    lty = "longdash")

# --- Axes ---
axis(1, at = rList[2:9], las = 1, cex.axis = 2.3)
axis(2, at = pseudo_log(yticks),
     labels = c("0", expression(10^0), expression(10^1), expression(10^2)),
     las = 1, cex.axis = 2.3)
box()

dev.off()





file = paste0("false_negatives_final-k", k, "-l", l, "-rho1", rho1, ".pdf")
setwd("~")
pdf(file, width = 8, height = 6)

# Données à tracer
y1 <- moyenne_faux_negatifsMed2[2:9,3]/((rList[2:9]/100)*n)*100
y2 <- moyenne_faux_negatifsMed2[2:9,2]/((rList[2:9]/100)*n)*100
y3 <- moyenne_faux_negatifsMed2[2:9,1]/((rList[2:9]/100)*n)*100
y4 <- moyenne_faux_negatifsMed2Oracle[2:9,1]/((rList[2:9]/100)*n)*100

# Transformation log10(y + 1) pour éviter les zéros
y1_log <- log10(y1 + 1)
y2_log <- log10(y2 + 1)
y3_log <- log10(y3 + 1)
y4_log <- log10(y4 + 1)

# Tracé principal
plot(rList[2:9], y1_log,
     type = "l", lwd = 4, col = "red",
     xlab = "", ylab = "",
     yaxt = "n", xaxt = "n",
     ylim = log10(c(1e-1, 1e2)) # échelle en puissance de 10
)
lines(rList[2:9], y2_log, lwd = 4, col = "blue", lty = "dashed")
lines(rList[2:9], y3_log, lwd = 4, col = "darkgreen", lty = "dotted")
lines(rList[2:9], y4_log, lwd = 4, col = "purple", lty = "longdash")

# Axe X
axis(1, at = rList[2:9], las = 1, cex.axis = 2.3)

# Axe Y avec puissances de 10
yticks_real <- 10^seq(-1, 2, by = 1)
axis(2, at = log10(yticks_real),
     labels = parse(text = paste0("10^", -1:2)),
     las = 1, cex.axis = 2.1)

box()
dev.off()


file = paste0("false_positives_final-k",k,"-l",l,"-rho1",rho1,".pdf")
setwd("~")
pdf(file,width = 8, height = 6) 


plot(rList[1:9], faux_positifsMed2[1:9,3,1]/((1 - rList[1:9]/100)*n)*100,
     type = "l", lwd = 4,col = "red", 
     xlab = "", ylab = "",   # Pas de label
     yaxt = "n", xaxt = "n", # On masque les axes par défaut
     #log = "y",              # Échelle logarithmique Y
     ylim = c(0, 20)     # Plage Y adaptée à tes ticks log
)
# 
# # Autres courbes
lines(rList[1:9], faux_positifsMed2[1:9,2,1]/((1 - rList[1:9]/100)*n)*100,lwd = 4,col = "blue", lty = "dashed")
lines(rList[1:9], faux_positifsMed2[1:9,1,1]/((1 - rList[1:9]/100)*n)*100,lwd = 4, col = "darkgreen", lty = "dotted")
lines(rList[1:9], faux_positifsOracleMed2[1:9,1]/((1 - rList[1:9]/100)*n)*100,lwd = 4, col = "orange", lty = "dotted")

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

setwd("~")

file <- paste0("erreurs_Sigma_final-k", k, "-l", l, "-rho1", rho1, ".pdf")

pdf(file, width = 8, height = 6)

# Plot principal — axe Y en log (automatique)
plot(rList[1:9], erreursSigmaMed2[n, 1:9, 3, 1],
     type = "l", lwd = 4, col = "red",
     xlab = "", ylab = "",
     yaxt = "n", xaxt = "n",
     log = "y",                       # ✅ active l’échelle log pour Y
     ylim = c(1e-1, 1e2)              # borne Y cohérente
)

# Autres courbes
lines(rList[1:9], erreursSigmaMed2[n, 1:9, 2, 1], lwd = 4, col = "blue", lty = "dashed")
lines(rList[1:9], erreursSigmaMed2[n,1:9 , 1, 1], lwd = 4, col = "darkgreen", lty = "dotted")

# Axe Y logarithmique lisible
log_ticks <- 10^seq(-1, 2, by = 1)
axis(2, at = log_ticks,
     labels = parse(text = paste0("10^", -1:2)),
     las = 1, cex.axis = 1.8)

# Axe X (je suppose que tu veux une échelle régulière)
x_ticks <- rList
axis(1, at = x_ticks, las = 1, cex.axis = 2.3)

box()  # bordure
dev.off()

k = k1val[5]; l = l1valup1[5]; rho1 = rho1val[5]
file = paste0("false_negatives_final-k", k, "-l", l, "-rho1", rho1, ".pdf")
setwd("~/figures")
pdf(file, width = 8, height = 6)
# --- Données transformées ---
epsilon <- 0
y_red   <- moyenne_faux_negatifsMed3[2:9,3]/((rList[2:9]/100)*n)*100
y_blue  <- moyenne_faux_negatifsMed3[2:9,2]/((rList[2:9]/100)*n)*100
y_green <- moyenne_faux_negatifsMed3[2:9,1]/((rList[2:9]/100)*n)*100
y_purp  <- moyenne_faux_negatifsMed3Oracle[2:9,1]/((rList[2:9]/100)*n)*100

# --- Transformation pseudo-log (linéaire entre 0 et 1, log10 au-delà) ---
pseudo_log <- function(y) {
  ifelse(y <= 1, y, 1 + log10(y))  # 0–1 linéaire, puis log
}

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
lines(rList[2:9], pseudo_log(y_purp),  lwd = 4, col = "purple",    lty = "longdash")

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

setwd("~")

plot(rList[1:9], faux_positifsMed3[1:9,3,1]/((1 - rList[1:9]/100)*n)*100,
     type = "l", lwd = 4,col = "red", 
     xlab = "", ylab = "",   # Pas de label
     yaxt = "n", xaxt = "n", # On masque les axes par défaut
     #log = "y",              # Échelle logarithmique Y
     ylim = c(0, 20)     # Plage Y adaptée à tes ticks log
)
# 
# # Autres courbes
lines(rList[1:9], faux_positifsMed3[1:9,2,1]/((1 - rList[1:9]/100)*n)*100,lwd = 4,col = "blue", lty = "dashed")
lines(rList[1:9], faux_positifsMed3[1:9,1,1]/((1 - rList[1:9]/100)*n)*100,lwd = 4, col = "darkgreen", lty = "dotted")
lines(rList[1:9], faux_positifsOracleMed3[1:9,1]/((1 - rList[1:9]/100)*n)*100,lwd = 4, col = "orange", lty = "dotted")

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




setwd("~")

file <- paste0("erreurs_Sigma_final-k", k, "-l", l, "-rho1", rho1, ".pdf")

pdf(file, width = 8, height = 6)

# Plot principal — axe Y en log (automatique)
plot(rList[1:9], erreursSigmaMed3[n, 1:9, 3, 1],
     type = "l", lwd = 4, col = "red",
     xlab = "", ylab = "",
     yaxt = "n", xaxt = "n",
     log = "y",                       # ✅ active l’échelle log pour Y
     ylim = c(1e-1, 1e2)              # borne Y cohérente
)

# Autres courbes
lines(rList[1:9], erreursSigmaMed3[n,1:9 , 2, 1], lwd = 4, col = "blue", lty = "dashed")
lines(rList[1:9], erreursSigmaMed3[n,1:9 , 1, 1], lwd = 4, col = "darkgreen", lty = "dotted")

# Axe Y logarithmique lisible
log_ticks <- 10^seq(-1, 2, by = 1)
axis(2, at = log_ticks,
     labels = parse(text = paste0("10^", -1:2)),
     las = 1, cex.axis = 1.8)

# Axe X (je suppose que tu veux une échelle régulière)
x_ticks <- rList[1:9]
axis(1, at = x_ticks, las = 1, cex.axis = 2.3)

box()  # bordure
dev.off()



# 
# k = k1val[6]; l = l1valup1[6];rho1 = rho1val[6]
# 
# file = paste0("false_negatives_final-k",k,"-l",l,"-rho1",rho1,".pdf")
# setwd("~")
# pdf(file,width = 8, height = 6) 
# 
# plot(rList[2:9], faux_negatifsMed4[2:9,3,1]/((rList[2:9]/100)*n)*100,
#      type = "l", lwd = 4,col = "red", 
#      xlab = "", ylab = "",   # Pas de label
#      yaxt = "n", xaxt = "n", # On masque les axes par défaut
#      #log = "y",              # Échelle logarithmique Y
#      ylim = c(0, 100)     # Plage Y adaptée à tes ticks log
# )
# # 
# # # Autres courbes
# lines(rList[2:9], faux_negatifsMed4[2:9,2,1]/((rList[2:9]/100)*n)*100,lwd = 4,col = "blue", lty = "dashed")
# lines(rList[2:9], faux_negatifsMed4[2:9,1,1]/((rList[2:9]/100)*n)*100,lwd = 4, col = "darkgreen", lty = "dotted")
# lines(rList[2:9], faux_negatifsOracleMed4[2:9,1]/((rList[2:9]/100)*n)*100,lwd = 4, col = "orange", lty = "dotted")
# 
# #  
# # # Axe Y logarithmique
# # log_ticks <- 10^seq(-1, 10, by = 1)
# # axis(2, at = log_ticks,
# #      labels = parse(text = paste0("10^", -1:10)),
# #      las = 1)
# x_ticks <- rList[2:9]
# axis(2, at = seq(0, 100, by = 5), las = 1,cex.axis = 2.3)
# 
# axis(1, at = x_ticks,las = 1,cex.axis = 2.3)
# dev.off()
# 
# 
# 
# file = paste0("false_positives_final-k",k,"-l",l,"-rho1",rho1,".pdf")
# setwd("~")
# pdf(file,width = 8, height = 6) 
# 
# setwd("~")
# 
# plot(rList[1:9], faux_positifsMed4[1:9,3,1]/((1 - rList[1:9]/100)*n)*100,
#      type = "l", lwd = 4,col = "red", 
#      xlab = "", ylab = "",   # Pas de label
#      yaxt = "n", xaxt = "n", # On masque les axes par défaut
#      #log = "y",              # Échelle logarithmique Y
#      ylim = c(0, 20)     # Plage Y adaptée à tes ticks log
# )
# # 
# # # Autres courbes
# lines(rList[1:9], faux_positifsMed4[1:9,2,1]/((1 - rList[1:9]/100)*n)*100,lwd = 4,col = "blue", lty = "dashed")
# lines(rList[1:9], faux_positifsMed4[1:9,1,1]/((1 - rList[1:9]/100)*n)*100,lwd = 4, col = "darkgreen", lty = "dotted")
# lines(rList[1:9], faux_positifsOracleMed4[1:9,1]/((1 - rList[1:9]/100)*n)*100,lwd = 4, col = "orange", lty = "dotted")
# 
# #  
# # # Axe Y logarithmique
# # log_ticks <- 10^seq(-1, 10, by = 1)
# # axis(2, at = log_ticks,
# #      labels = parse(text = paste0("10^", -1:10)),
# #      las = 1)
# x_ticks <- rList[1:9]
# axis(2, at = seq(0, 100, by = 5), las = 1,cex.axis = 2.3)
# 
# axis(1, at = x_ticks,las = 1,cex.axis = 2.3)
# dev.off()
# 
# 
# 
# file = paste0("erreurs_Sigma_final-k",k,"-l",l,"-rho1",rho1,".pdf")
# setwd("~")
# pdf(file,width = 8, height = 6) 
# 
# 
# # Plot principal — axe Y en log (automatique)
# plot(rList[1:9], erreursSigmaMed4[n,1:9 , 3, 1],
#      type = "l", lwd = 4, col = "red",
#      xlab = "", ylab = "",
#      yaxt = "n", xaxt = "n",
#      log = "y",                       # ✅ active l’échelle log pour Y
#      ylim = c(1e-1, 1e2)              # borne Y cohérente
# )
# 
# # Autres courbes
# lines(rList[1:9], erreursSigmaMed4[n, 1:9, 2, 1], lwd = 4, col = "blue", lty = "dashed")
# lines(rList[1:9], erreursSigmaMed4[n, 1:9, 1, 1], lwd = 4, col = "darkgreen", lty = "dotted")
# 
# # Axe Y logarithmique lisible
# log_ticks <- 10^seq(-1, 2, by = 1)
# axis(2, at = log_ticks,
#      labels = parse(text = paste0("10^", -1:2)),
#      las = 1, cex.axis = 1.8)
# 
# # Axe X (je suppose que tu veux une échelle régulière)
# x_ticks <- rList
# axis(1, at = x_ticks, las = 1, cex.axis = 2.3)
# 
# box()  # bordure
# dev.off()
# 
# 
# 
# k = k1val[7]; l = l1valup1[7];rho1 = rho1val[7]
# 
# file = paste0("false_negatives_final-k",k,"-l",l,"-rho1",rho1,".pdf")
# setwd("~")
# pdf(file,width = 8, height = 6) 
# 
# plot(rList[2:9], faux_negatifsMed5[2:9,3,1]/((rList[2:9]/100)*n)*100,
#      type = "l", lwd = 4,col = "red", 
#      xlab = "", ylab = "",   # Pas de label
#      yaxt = "n", xaxt = "n", # On masque les axes par défaut
#      #log = "y",              # Échelle logarithmique Y
#      ylim = c(0, 100)     # Plage Y adaptée à tes ticks log
# )
# # 
# # # Autres courbes
# lines(rList[2:9], faux_negatifsMed5[2:9,2,1]/((rList[2:9]/100)*n)*100,lwd = 4,col = "blue", lty = "dashed")
# lines(rList[2:9], faux_negatifsMed5[2:9,1,1]/((rList[2:9]/100)*n)*100,lwd = 4, col = "darkgreen", lty = "dotted")
# lines(rList[2:9], faux_negatifsOracleMed5[2:9,1]/((rList[2:9]/100)*n)*100,lwd = 4, col = "orange", lty = "dotted")
# 
# #  
# # # Axe Y logarithmique
# # log_ticks <- 10^seq(-1, 10, by = 1)
# # axis(2, at = log_ticks,
# #      labels = parse(text = paste0("10^", -1:10)),
# #      las = 1)
# x_ticks <- rList[2:9]
# axis(2, at = seq(0, 100, by = 5), las = 1,cex.axis = 2.3)
# 
# axis(1, at = x_ticks,las = 1,cex.axis = 2.3)
# dev.off()
# 
# 
# file = paste0("erreurs_Sigma_final-k",k,"-l",l,"-rho1",rho1,".pdf")
# setwd("~")
# pdf(file,width = 8, height = 6) 
# 
# 
# # Plot principal — axe Y en log (automatique)
# plot(rList[1:9], erreursSigmaMed5[n,1:9 , 3, 1],
#      type = "l", lwd = 4, col = "red",
#      xlab = "", ylab = "",
#      yaxt = "n", xaxt = "n",
#      log = "y",                       
#      ylim = c(1e-1, 1e2)              
# )
# 
# # Autres courbes
# lines(rList[1:9], erreursSigmaMed5[n,1:9 , 2, 1], lwd = 4, col = "blue", lty = "dashed")
# lines(rList[1:9], erreursSigmaMed5[n, 1:9, 1, 1], lwd = 4, col = "darkgreen", lty = "dotted")
# 
# # Axe Y logarithmique lisible
# log_ticks <- 10^seq(-1, 2, by = 1)
# axis(2, at = log_ticks,
#      labels = parse(text = paste0("10^", -1:2)),
#      las = 1, cex.axis = 1.8)
# 
# # Axe X (je suppose que tu veux une échelle régulière)
# x_ticks <- rList
# axis(1, at = x_ticks, las = 1, cex.axis = 2.3)
# 
# box()  # bordure
# dev.off()
# 
# file = paste0("false_positives_final-k",k,"-l",l,"-rho1",rho1,".pdf")
# setwd("~")
# pdf(file,width = 8, height = 6) 
# 
# setwd("~")
# 
# plot(rList[1:9], faux_positifsMed5[1:9,3,1]/((1 - rList[1:9]/100)*n)*100,
#      type = "l", lwd = 4,col = "red", 
#      xlab = "", ylab = "",   # Pas de label
#      yaxt = "n", xaxt = "n", # On masque les axes par défaut
#      #log = "y",              # Échelle logarithmique Y
#      ylim = c(0, 20)     # Plage Y adaptée à tes ticks log
# )
# # 
# # # Autres courbes
# lines(rList[1:9], faux_positifsMed5[1:9,2,1]/((1 - rList[1:9]/100)*n)*100,lwd = 4,col = "blue", lty = "dashed")
# lines(rList[1:9], faux_positifsMed5[1:9,1,1]/((1 - rList[1:9]/100)*n)*100,lwd = 4, col = "darkgreen", lty = "dotted")
# lines(rList[1:9], faux_positifsOracleMed5[1:9,1]/((1 - rList[1:9]/100)*n)*100,lwd = 4, col = "orange", lty = "dotted")
# 
# #  
# # # Axe Y logarithmique
# # log_ticks <- 10^seq(-1, 10, by = 1)
# # axis(2, at = log_ticks,
# #      labels = parse(text = paste0("10^", -1:10)),
# #      las = 1)
# x_ticks <- rList[1:9]
# axis(2, at = seq(0, 100, by = 5), las = 1,cex.axis = 2.3)
# 
# axis(1, at = x_ticks,las = 1,cex.axis = 2.3)
# dev.off()
# 
# 
# 
# 
# 
# ##################################Cumulative outlier detection###################################
# 
# ########################Cumulative outlier detection###################
# 
# 
# resFar5Naive = cumulativeOutlierDetection(labelsVraisFar[,2],majority_vote_Far[,2,1],5)
# resFar5Online = cumulativeOutlierDetection(labelsVraisFar[,2],majority_vote_Far[,2,2],5)
# resFar5Streaming = cumulativeOutlierDetection(labelsVraisFar[,2],majority_vote_Far[,2,3],5)
# 
# # Données supposées : n correspond à la taille de resFar5Naive$taux_outliers_detectes_vrais
# x_vals <- 1:n
# 
# setwd("~")
# 
# file = paste0("cumulativeOutlierDetFarr5",".png")
# png(file, width = 1800, height = 1200, res = 200)
# # --- Plot principal (échelle log sur X)
# plot(x_vals, resFar5Naive$taux_outliers_detectes_vrais,
#      type = "l", lwd = 4, lty = "dotted",
#      xlab = "", ylab = "",
#      yaxt = "n", xaxt = "n",
#      ylim = c(0, 100),
#      col = "darkgreen",
#      log = "x"        # <=== Échelle logarithmique sur l’axe X
# )
# 
# # --- Autres lignes
# lines(x_vals, resFar5Online$taux_outliers_detectes_vrais,
#       lwd = 4, col = adjustcolor("blue",alpha.f = 0.5), lty = "dashed")
# lines(x_vals, resFar5Streaming$taux_outliers_detectes_vrais,
#       lwd = 4, col = adjustcolor("red",alpha.f = 0.8), lty = "solid")
# 
# # Ligne invisible (orange)
# lines(x_vals, resFar5Naive$taux_outliers_vrais,
#       lwd = 4, col = adjustcolor("orange", alpha.f = 0), type = "l")
# 
# # --- Axes manuels
# 
# # Axe X (logarithmique)
# log_ticks_x <- 10^seq(0, floor(log10(n)), by = 1)   # 1, 10, 100, 1000, ...
# axis(1, at = log_ticks_x,
#      labels = parse(text = paste0("10^", seq(0, floor(log10(n))))),
#      las = 1, cex.axis = 2)
# 
# # Axe Y (linéaire)
# axis(2, at = seq(0, 100, by = 10), las = 1, cex.axis = 2.5)
# dev.off() 
# 
# setwd("~")
# 
# file = paste0("cumulativeOutlierDetFarr20",".png")
# png(file, width = 1800, height = 1200, res = 200)
# resFar20Naive = cumulativeOutlierDetection(labelsVraisFar[,5],majority_vote_Far[,5,1],20)
# resFar20Online = cumulativeOutlierDetection(labelsVraisFar[,5],majority_vote_Far[,5,2],20)
# resFar20Streaming = cumulativeOutlierDetection(labelsVraisFar[,5],majority_vote_Far[,5,3],20)
# x_vals <- 1:n
# 
# # --- Plot principal (échelle log sur X)
# plot(x_vals, resFar20Naive$taux_outliers_detectes_vrais,
#      type = "l", lwd = 4, lty = "dotted",
#      xlab = "", ylab = "",
#      yaxt = "n", xaxt = "n",
#      ylim = c(0, 100),
#      col = "darkgreen",
#      log = "x"        # <=== Échelle logarithmique sur l’axe X
# )
# 
# # --- Autres lignes
# lines(x_vals, resFar20Online$taux_outliers_detectes_vrais,
#       lwd = 4, col = adjustcolor("blue",alpha.f = 0.5), lty = "dashed")
# lines(x_vals, resFar20Streaming$taux_outliers_detectes_vrais,
#       lwd = 4, col = adjustcolor("red",alpha.f = 0.8), lty = "solid")
# 
# # Ligne invisible (orange)
# lines(x_vals, resFar20Naive$taux_outliers_vrais,
#       lwd = 4, col = adjustcolor("orange", alpha.f = 0), type = "l")
# 
# # --- Axes manuels
# 
# # Axe X (logarithmique)
# log_ticks_x <- 10^seq(0, floor(log10(n)), by = 1)   # 1, 10, 100, 1000, ...
# axis(1, at = log_ticks_x,
#      labels = parse(text = paste0("10^", seq(0, floor(log10(n))))),
#      las = 1, cex.axis = 2)
# 
# # Axe Y (linéaire)
# axis(2, at = seq(0, 100, by = 10), las = 1, cex.axis = 2.5)
# 
# dev.off()
# 
# setwd("~")
# 
# file = paste0("cumulativeOutlierDetFarr30",".png")
# png(file, width = 1800, height = 1200, res = 200)
# 
# resFar30Naive = cumulativeOutlierDetection(labelsVraisFar[,7],majority_vote_Far[,7,1],30)
# resFar30Online = cumulativeOutlierDetection(labelsVraisFar[,7],majority_vote_Far[,7,2],30)
# resFar30Streaming = cumulativeOutlierDetection(labelsVraisFar[,7],majority_vote_Far[,7,3],30)
# 
# # --- Plot principal (échelle log sur X)
# plot(x_vals, resFar30Naive$taux_outliers_detectes_vrais,
#      type = "l", lwd = 4, lty = "dotted",
#      xlab = "", ylab = "",
#      yaxt = "n", xaxt = "n",
#      ylim = c(0, 100),
#      col = "darkgreen",
#      log = "x"        # <=== Échelle logarithmique sur l’axe X
# )
# 
# # --- Autres lignes
# lines(x_vals, resFar30Online$taux_outliers_detectes_vrais,
#       lwd = 4, col = adjustcolor("blue",alpha.f = 0.5), lty = "dashed")
# lines(x_vals, resFar30Streaming$taux_outliers_detectes_vrais,
#       lwd = 4, col = adjustcolor("red",alpha.f = 0.8), lty = "solid")
# 
# # Ligne invisible (orange)
# lines(x_vals, resFar30Naive$taux_outliers_vrais,
#       lwd = 4, col = adjustcolor("orange", alpha.f = 0.7), type = "l")
# 
# # --- Axes manuels
# 
# # Axe X (logarithmique)
# log_ticks_x <- 10^seq(0, floor(log10(n)), by = 1)   # 1, 10, 100, 1000, ...
# axis(1, at = log_ticks_x,
#      labels = parse(text = paste0("10^", seq(0, floor(log10(n))))),
#      las = 1, cex.axis = 2)
# 
# # Axe Y (linéaire)
# axis(2, at = seq(0, 100, by = 10), las = 1, cex.axis = 2.5)
# 
# dev.off()
# 
# resFar40Naive = cumulativeOutlierDetection(labelsVraisFar[,11],majority_vote_Far[,11,1],50)
# resFar40Online = cumulativeOutlierDetection(labelsVraisFar[,11],majority_vote_Far[,11,2],50)
# resFar40Streaming = cumulativeOutlierDetection(labelsVraisFar[,11],majority_vote_Far[,11,3],50)
# 
# setwd("~")
# 
# file = paste0("cumulativeOutlierDetFarr40",".png")
# png(file, width = 1800, height = 1200, res = 200)
# 
# 
# # --- Plot principal (échelle log sur X)
# 
# 
# 
# plot(x_vals, resFar40Naive$taux_outliers_detectes_vrais,
#      type = "l", lwd = 4, lty = "dotted",
#      xlab = "", ylab = "",
#      yaxt = "n", xaxt = "n",
#      ylim = c(0, 100),
#      col = "darkgreen",
#      log = "x"        # <=== Échelle logarithmique sur l’axe X
# )
# 
# # --- Autres lignes
# lines(x_vals, resFar40Online$taux_outliers_detectes_vrais,
#       lwd = 4, col = adjustcolor("blue",alpha.f = 0.5), lty = "dashed")
# lines(x_vals, resFar40Streaming$taux_outliers_detectes_vrais,
#       lwd = 4, col = adjustcolor("red",alpha.f = 0.8), lty = "solid")
# 
# # Ligne invisible (orange)
# lines(x_vals, resFar40Naive$taux_outliers_vrais,
#       lwd = 4, col = adjustcolor("orange", alpha.f = 0.7), type = "l")
# 
# # --- Axes manuels
# 
# # Axe X (logarithmique)
# log_ticks_x <- 10^seq(0, floor(log10(n)), by = 1)   # 1, 10, 100, 1000, ...
# axis(1, at = log_ticks_x,
#      labels = parse(text = paste0("10^", seq(0, floor(log10(n))))),
#      las = 1, cex.axis = 2)
# 
# # Axe Y (linéaire)
# axis(2, at = seq(0, 100, by = 10), las = 1, cex.axis = 2.5)
# 
# dev.off()
# 
# 
# # ##############################################
# # Temps calculs
# # ################# #############################
# 
# temps_calculTout = res$temps_calcul
# 
# 
# # Méthodes et indices souhaités
# methodes <- c("sampleCovOnline", "samplecovTrimmed", "sampleCovOffline", "comedianeOffline",
#               "comedianeOfflineShrinkage", "OGK", "FASTMCD", "offline", "online", "streaming")
# 
# # Indices à garder : 1, 2, 6 à 10
# indices_gardes <- c(1, 6, 7, 8, 9, 10)
# 
# # Supposons que taux_index est défini
# taux_index <- 3  # par exemple
# 
# # Extraction des données pour ce taux et méthodes sélectionnées
# temps_sel <- temps_calculTout[taux_index, indices_gardes, ]  # dims : méthodes sélectionnées x runs
# 
# # Transformation en data frame long
# df_temps <- melt(temps_sel, varnames = c("MethodeIndex", "Run"), value.name = "Temps")
# 
# # Remplacement par les noms des méthodes sélectionnées
# df_temps$Methode <- factor(df_temps$MethodeIndex, 
#                            levels = 1:length(indices_gardes), 
#                            labels = methodes[indices_gardes])
# 
# # Plot boxplot
# ggplot(df_temps, aes(x = Methode, y = Temps)) +
#   geom_boxplot(fill = "lightblue", outlier.color = "red", outlier.shape = 1) +
#   labs(
#     title = "Boxplot of computation times",
#     x = "Method",
#     y = "Time (seconds)"
#   ) +
#   scale_y_log10() +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
