##########################################
#Erreurs KL k,l rho1 false positives, false negatives
#############################################
setwd("~")
pdf("frobnormerrorklkd10.pdf", width = 7, height = 5)

#######################Plot erreurs Sigma k#############################
plot(rList, erreursSigmak[,4,3],
     type = "l", lwd = 4, col = "red",
     log = "y",              # axe Y en log10
     xlab = "", ylab = "",   # pas de noms d'axes
     yaxt = "n", xaxt = "n",
     ylim = c(1e-1, 1e2)) # on redessine les axes manuellement

# Axe Y en log10 avec labels 10^k
# ylim <- range(erreursSigmak[,4,3], na.rm = TRUE)
# ylim[1] <- max(ylim[1], min(erreursSigmaRho1[erreursSigmalsup1[,4,3] > 0,4,3], na.rm = TRUE)) # > 0
#############KL 100############################
lines(rList,erreursSigmak[,4,2],col = "red",lwd = 4,lty = "dashed")
lines(rList,erreursSigmak[,4,1],col = "red",lwd = 4,lty = "dotted")
#############KL 10############################
lines(rList,erreursSigmak[,3,3],col = "orange",lwd = 4)
lines(rList,erreursSigmak[,3,2],col = "orange",lty = "dashed",lwd = 4)
lines(rList,erreursSigmak[,3,1],col = "orange",lty = "dotted",lwd = 4)
#############KL 1############################
lines(rList,erreursSigmak[,2,3],col = "darkgreen",lwd = 4)
lines(rList,erreursSigmak[,2,2],col = "darkgreen",lty = "dashed",lwd = 4)
lines(rList,erreursSigmak[,2,1],col = "darkgreen",lty = "dotted",lwd = 4)
#############KL 0############################
lines(rList,erreursSigmak[,1,3],col = "blue",lwd = 4)
lines(rList,erreursSigmak[,1,2],col = "blue",lty = "dashed",lwd = 4)
lines(rList,erreursSigmak[,1,1],col = "blue",lty = "dotted",lwd = 4)


# log_ticks <- 10^seq(floor(log10(ylim[1])), ceiling(log10(ylim[2])), by = 1)
# axis(2, at = log_ticks,
#      labels = parse(text = paste0("10^", seq(floor(log10(ylim[1])),
#                                              ceiling(log10(ylim[2])), 1))))


log_ticks <- 10^seq(-1, 2, by = 1)  # ticks de 10^-2 à 10^5
axis(2, at = log_ticks,
     labels = parse(text = paste0("10^", -1:2)),cex.axis = 2.3)
# Axe X tous les 50
axis(1, at = seq(min(rList), max(rList), by = 5),cex.axis= 2.3)
dev.off()


setwd("~")
pdf("frobnormerrorkllsup1d10.pdf", width = 7, height = 5)

#######################Plot erreurs Sigma lsup 1#############################
plot(rList, erreursSigmalsup1[,4,3],
     type = "l", lwd = 4, col = "red",
     log = "y",              # axe Y en log10
     xlab = "", ylab = "",   # pas de noms d'axes
     yaxt = "n", xaxt = "n",
     ylim = c(1e-1, 1e2)) # on redessine les axes manuellement

# Axe Y en log10 avec labels 10^k
# ylim <- range(erreursSigmak[,4,3], na.rm = TRUE)
# ylim[1] <- max(ylim[1], min(erreursSigmaRho1[erreursSigmalsup1[,4,3] > 0,4,3], na.rm = TRUE)) # > 0
#############KL 100############################
lines(rList,erreursSigmalsup1[,4,2],col = "red",lwd = 4,lty = "dashed")
lines(rList,erreursSigmalsup1[,4,1],col = "red",lwd = 4,lty = "dotted")
#############KL 10############################
lines(rList,erreursSigmalsup1[,3,3],col = "orange",lwd = 4)
lines(rList,erreursSigmalsup1[,3,2],col = "orange",lty = "dashed",lwd = 4)
lines(rList,erreursSigmalsup1[,3,1],col = "orange",lty = "dotted",lwd = 4)
#############KL 1############################
lines(rList,erreursSigmalsup1[,2,3],col = "darkgreen",lwd = 4)
lines(rList,erreursSigmalsup1[,2,2],col = "darkgreen",lty = "dashed",lwd = 4)
lines(rList,erreursSigmalsup1[,2,1],col = "darkgreen",lty = "dotted",lwd = 4)
#############KL 0############################
lines(rList,erreursSigmalsup1[,1,3],col = "blue",lwd = 4)
lines(rList,erreursSigmalsup1[,1,2],col = "blue",lty = "dashed",lwd = 4)
lines(rList,erreursSigmalsup1[,1,1],col = "blue",lty = "dotted",lwd = 4)


# log_ticks <- 10^seq(floor(log10(ylim[1])), ceiling(log10(ylim[2])), by = 1)
# axis(2, at = log_ticks,
#      labels = parse(text = paste0("10^", seq(floor(log10(ylim[1])),
#                                              ceiling(log10(ylim[2])), 1))))


log_ticks <- 10^seq(-1, 2, by = 1)  # ticks de 10^-2 à 10^5
axis(2, at = log_ticks,
     labels = parse(text = paste0("10^", -1:2)),cex.axis = 2.3)
# Axe X tous les 50
axis(1, at = seq(min(rList), max(rList), by = 5),cex.axis= 2.3)

dev.off()

setwd("~")
pdf("frobnormerrorklrho1d10.pdf", width = 7, height = 5)


#######################Plot erreurs Sigma Rho1#############################
plot(rList, erreursSigmaRho1[,4,3],
     type = "l", lwd = 4, col = "red",
     log = "y",              # axe Y en log10
     xlab = "", ylab = "",   # pas de noms d'axes
     yaxt = "n", xaxt = "n",
     ylim = c(1e-1, 1e2)) # on redessine les axes manuellement

# Axe Y en log10 avec labels 10^k
# ylim <- range(erreursSigmak[,4,3], na.rm = TRUE)
# ylim[1] <- max(ylim[1], min(erreursSigmaRho1[erreursSigmalsup1[,4,3] > 0,4,3], na.rm = TRUE)) # > 0
#############KL 100############################
lines(rList,erreursSigmaRho1[,4,2],col = "red",lwd = 4,lty = "dashed")
lines(rList,erreursSigmaRho1[,4,1],col = "red",lwd = 4,lty = "dotted")
#############KL 10############################
lines(rList,erreursSigmaRho1[,3,3],col = "orange",lwd = 4)
lines(rList,erreursSigmaRho1[,3,2],col = "orange",lty = "dashed",lwd = 4)
lines(rList,erreursSigmaRho1[,3,1],col = "orange",lty = "dotted",lwd = 4)
#############KL 1############################
lines(rList,erreursSigmaRho1[,2,3],col = "darkgreen",lwd = 4)
lines(rList,erreursSigmaRho1[,2,2],col = "darkgreen",lty = "dashed",lwd = 4)
lines(rList,erreursSigmaRho1[,2,1],col = "darkgreen",lty = "dotted",lwd = 4)
#############KL 0############################
lines(rList,erreursSigmaRho1[,1,3],col = "blue",lwd = 4)
lines(rList,erreursSigmaRho1[,1,2],col = "blue",lty = "dashed",lwd = 4)
lines(rList,erreursSigmaRho1[,1,1],col = "blue",lty = "dotted",lwd = 4)


# log_ticks <- 10^seq(floor(log10(ylim[1])), ceiling(log10(ylim[2])), by = 1)
# axis(2, at = log_ticks,
#      labels = parse(text = paste0("10^", seq(floor(log10(ylim[1])),
#                                              ceiling(log10(ylim[2])), 1))))


log_ticks <- 10^seq(-1, 2, by = 1)  # ticks de 10^-2 à 10^5
axis(2, at = log_ticks,
     labels = parse(text = paste0("10^", -1:2)),cex.axis = 2.3)
# Axe X tous les 50
axis(1, at = seq(min(rList), max(rList), by = 5),cex.axis= 2.3)

dev.off()

setwd("~")

pdf("falsenegativesklkd10.pdf", width = 7, height = 5)

#######################Plot faux négatifs k#############################
plot(rList[2:length(rList)], faux_negatifsRatek[(2:length(rList)),4,3],
     type = "l", lwd = 4, col = "red",
     #log = "y",              # axe Y en log10
     xlab = "", ylab = "",   # pas de noms d'axes
     yaxt = "n", xaxt = "n",
     #ylim = c(0, 1e2)) # on redessine les axes manuellement
     ylim = c(0, 1e2)
)
# Axe Y en log10 avec labels 10^k
# ylim <- range(faux_negatifsRatek[,4,3], na.rm = TRUE)
# ylim[1] <- max(ylim[1], min(faux_negatifsRatek[faux_negatifsRatek[,4,3] > 0,4,3], na.rm = TRUE)) # > 0
#############KL 100############################
lines(rList[2:length(rList)],faux_negatifsRateRho1[(2:length(rList)),4,2],col = "red",lwd = 4,lty = "dashed")
lines(rList[2:length(rList)],faux_negatifsRateRho1[(2:length(rList)),4,1],col = "red",lwd = 4,lty = "dotted")
#############KL 10############################
lines(rList[2:length(rList)],faux_negatifsRateRho1[(2:length(rList)),3,3],lwd = 4,col = "orange")
lines(rList[2:length(rList)],faux_negatifsRateRho1[(2:length(rList)),3,2],lwd = 4,col = "orange",lty = "dashed")
lines(rList[2:length(rList)],faux_negatifsRateRho1[(2:length(rList)),3,1],lwd = 4,col = "orange",lty = "dotted")
#############KL 1############################
lines(rList[2:length(rList)],faux_negatifsRateRho1[(2:length(rList)),2,3],lwd = 4,col = "darkgreen")
lines(rList[2:length(rList)],faux_negatifsRateRho1[(2:length(rList)),2,2],lwd = 4,col = "darkgreen",lty = "dashed")
lines(rList[2:length(rList)],faux_negatifsRateRho1[(2:length(rList)),2,1],lwd = 4,col = "darkgreen",lty = "dotted")
#############KL 0############################
lines(rList[2:length(rList)],faux_negatifsRateRho1[(2:length(rList)),1,3],lwd = 4,col = "blue")
lines(rList[2:length(rList)],faux_negatifsRateRho1[(2:length(rList)),1,2],,lwd = 4,col = "blue",lty = "dashed")
lines(rList[2:length(rList)],faux_negatifsRateRho1[(2:length(rList)),1,1],lwd = 4,col = "blue",lty = "dotted")


# log_ticks <- 10^seq(floor(log10(ylim[1])), ceiling(log10(ylim[2])), by = 1)
# axis(2, at = log_ticks,
#      labels = parse(text = paste0("10^", seq(floor(log10(ylim[1])),
#                                              ceiling(log10(ylim[2])), 1))))


# log_ticks <- 10^seq(-2, 2, by = 1)  # ticks de 10^-2 à 10^5
# axis(2, at = log_ticks,
#      labels = parse(text = paste0("10^", -2:2)))
# # Axe X tous les 50
axis(1, at = seq(5, max(rList), by = 5),cex.axis = 2.3)



# Axe X tous les 50

axis(2, at = seq(0, 100, by = 5),las = 1,cex.axis = 1.8)


dev.off()



setwd("~")

pdf("falsenegativeskllsup1d10.pdf", width = 7, height = 5)

#######################Plot faux négatifs l >= 1#############################
plot(rList[2:length(rList)], faux_negatifsRatelsup1[(2:length(rList)),4,3],
     type = "l", lwd = 4, col = "red",
     #log = "y",              # axe Y en log10
     xlab = "", ylab = "",   # pas de noms d'axes
     yaxt = "n", xaxt = "n",
     #ylim = c(0, 1e2)) # on redessine les axes manuellement
     ylim = c(0, 1e2)
)
# Axe Y en log10 avec labels 10^k
# ylim <- range(faux_negatifsRatek[,4,3], na.rm = TRUE)
# ylim[1] <- max(ylim[1], min(faux_negatifsRatek[faux_negatifsRatek[,4,3] > 0,4,3], na.rm = TRUE)) # > 0
#############KL 100############################
lines(rList[2:length(rList)],faux_negatifsRatelsup1[(2:length(rList)),4,2],col = "red",lwd = 4,lty = "dashed")
lines(rList[2:length(rList)],faux_negatifsRatelsup1[(2:length(rList)),4,1],col = "red",lwd = 4,lty = "dotted")
#############KL 10############################
lines(rList[2:length(rList)],faux_negatifsRatelsup1[(2:length(rList)),3,3],lwd = 4,col = "orange")
lines(rList[2:length(rList)],faux_negatifsRatelsup1[(2:length(rList)),3,2],lwd = 4,col = "orange",lty = "dashed")
lines(rList[2:length(rList)],faux_negatifsRatelsup1[(2:length(rList)),3,1],lwd = 4,col = "orange",lty = "dotted")
#############KL 1############################
lines(rList[2:length(rList)],faux_negatifsRatelsup1[(2:length(rList)),2,3],lwd = 4,col = "darkgreen")
lines(rList[2:length(rList)],faux_negatifsRatelsup1[(2:length(rList)),2,2],lwd = 4,col = "darkgreen",lty = "dashed")
lines(rList[2:length(rList)],faux_negatifsRatelsup1[(2:length(rList)),2,1],lwd = 4,col = "darkgreen",lty = "dotted")
#############KL 0############################
lines(rList[2:length(rList)],faux_negatifsRatelsup1[(2:length(rList)),1,3],lwd = 4,col = "blue")
lines(rList[2:length(rList)],faux_negatifsRatelsup1[(2:length(rList)),1,2],,lwd = 4,col = "blue",lty = "dashed")
lines(rList[2:length(rList)],faux_negatifsRatelsup1[(2:length(rList)),1,1],lwd = 4,col = "blue",lty = "dotted")


# log_ticks <- 10^seq(floor(log10(ylim[1])), ceiling(log10(ylim[2])), by = 1)
# axis(2, at = log_ticks,
#      labels = parse(text = paste0("10^", seq(floor(log10(ylim[1])),
#                                              ceiling(log10(ylim[2])), 1))))


# log_ticks <- 10^seq(-2, 2, by = 1)  # ticks de 10^-2 à 10^5
# axis(2, at = log_ticks,
#      labels = parse(text = paste0("10^", -2:2)))
# # Axe X tous les 50
axis(1, at = seq(5, max(rList), by = 5),cex.axis = 2.3)



# Axe X tous les 50

axis(2, at = seq(0, 100, by = 5),las = 1,cex.axis = 1.8)


dev.off()


setwd("~")

pdf("falsenegativesklrho1d10.pdf", width = 7, height = 5)

#######################Plot faux négatifs l >= 1#############################
plot(rList[2:length(rList)], faux_negatifsRateRho1[(2:length(rList)),4,3],
     type = "l", lwd = 4, col = "red",
     #log = "y",              # axe Y en log10
     xlab = "", ylab = "",   # pas de noms d'axes
     yaxt = "n", xaxt = "n",
     #ylim = c(0, 1e2)) # on redessine les axes manuellement
     ylim = c(0, 1e2)
)
# Axe Y en log10 avec labels 10^k
# ylim <- range(faux_negatifsRatek[,4,3], na.rm = TRUE)
# ylim[1] <- max(ylim[1], min(faux_negatifsRatek[faux_negatifsRatek[,4,3] > 0,4,3], na.rm = TRUE)) # > 0
#############KL 100############################
lines(rList[2:length(rList)],faux_negatifsRateRho1[(2:length(rList)),4,2],col = "red",lwd = 4,lty = "dashed")
lines(rList[2:length(rList)],faux_negatifsRateRho1[(2:length(rList)),4,1],col = "red",lwd = 4,lty = "dotted")
#############KL 10############################
lines(rList[2:length(rList)],faux_negatifsRateRho1[(2:length(rList)),3,3],lwd = 4,col = "orange")
lines(rList[2:length(rList)],faux_negatifsRateRho1[(2:length(rList)),3,2],lwd = 4,col = "orange",lty = "dashed")
lines(rList[2:length(rList)],faux_negatifsRateRho1[(2:length(rList)),3,1],lwd = 4,col = "orange",lty = "dotted")
#############KL 1############################
lines(rList[2:length(rList)],faux_negatifsRateRho1[(2:length(rList)),2,3],lwd = 4,col = "darkgreen")
lines(rList[2:length(rList)],faux_negatifsRateRho1[(2:length(rList)),2,2],lwd = 4,col = "darkgreen",lty = "dashed")
lines(rList[2:length(rList)],faux_negatifsRateRho1[(2:length(rList)),2,1],lwd = 4,col = "darkgreen",lty = "dotted")
#############KL 0############################
lines(rList[2:length(rList)],faux_negatifsRateRho1[(2:length(rList)),1,3],lwd = 4,col = "blue")
lines(rList[2:length(rList)],faux_negatifsRateRho1[(2:length(rList)),1,2],,lwd = 4,col = "blue",lty = "dashed")
lines(rList[2:length(rList)],faux_negatifsRateRho1[(2:length(rList)),1,1],lwd = 4,col = "blue",lty = "dotted")


# log_ticks <- 10^seq(floor(log10(ylim[1])), ceiling(log10(ylim[2])), by = 1)
# axis(2, at = log_ticks,
#      labels = parse(text = paste0("10^", seq(floor(log10(ylim[1])),
#                                              ceiling(log10(ylim[2])), 1))))


# log_ticks <- 10^seq(-2, 2, by = 1)  # ticks de 10^-2 à 10^5
# axis(2, at = log_ticks,
#      labels = parse(text = paste0("10^", -2:2)))
# # Axe X tous les 50
axis(1, at = seq(5, max(rList), by = 5),cex.axis = 2.3)



# Axe X tous les 50

axis(2, at = seq(0, 100, by = 5),las = 1,cex.axis = 1.8)


dev.off()


setwd("~")

pdf("falsepositivesklkd10.pdf", width = 7, height = 5)


#######################Plot faux positifs k#############################
plot(rList, faux_positifsRatek[,4,3],
     type = "l", lwd = 4, col = "red",
     #log = "y",              # axe Y en log10
     xlab = "", ylab = "",   # pas de noms d'axes
     yaxt = "n", xaxt = "n",
     #ylim = c(0, 1e2)) # on redessine les axes manuellement
     ylim = c(0, 20)
)
# Axe Y en log10 avec labels 10^k
# ylim <- range(faux_negatifsRatek[,4,3], na.rm = TRUE)
# ylim[1] <- max(ylim[1], min(faux_negatifsRatek[faux_negatifsRatek[,4,3] > 0,4,3], na.rm = TRUE)) # > 0
#############KL 100############################
lines(rList,faux_positifsRatek[,4,2],lwd = 4,col = "red",lty = "dashed")
lines(rList,faux_positifsRatek[,4,1],lwd = 4,col = "red",lty = "dotted")
#############KL 10############################
lines(rList,faux_positifsRatek[,3,3],lwd = 4,col = "orange")
lines(rList,faux_positifsRatek[,3,2],lwd = 4,col = "orange",lty = "dashed")
lines(rList,faux_positifsRatek[,3,1],lwd = 4,col = "orange",lty = "dotted")
#############KL 1############################
lines(rList,faux_positifsRatek[,2,3],lwd = 4,col = "darkgreen")
lines(rList,faux_positifsRatek[,2,2],lwd =4,col = "darkgreen",lty = "dashed")
lines(rList,faux_positifsRatek[,2,1],lwd = 4,col = "darkgreen",lty = "dotted")
#############KL 0############################
lines(rList,faux_positifsRatek[,1,3],lwd = 4,col = "blue")
lines(rList,faux_positifsRatek[,1,2],lwd = 4,col = "blue",lty = "dashed")
lines(rList,faux_positifsRatek[,1,1],lwd = 4,col = "blue",lty = "dotted")


# log_ticks <- 10^seq(floor(log10(ylim[1])), ceiling(log10(ylim[2])), by = 1)
# axis(2, at = log_ticks,
#      labels = parse(text = paste0("10^", seq(floor(log10(ylim[1])),
#                                              ceiling(log10(ylim[2])), 1))))


# log_ticks <- 10^seq(-2, 2, by = 1)  # ticks de 10^-2 à 10^5
# axis(2, at = log_ticks,
#      labels = parse(text = paste0("10^", -2:2)))
# # Axe X tous les 50
axis(1, at = seq(0, max(rList), by = 5),cex.axis = 2.3)



# Axe X tous les 50

axis(2, at = seq(0, 100, by = 5),las = 1,cex.axis = 2.3)


dev.off()



setwd("~")

pdf("falsepositiveskllsup1d10.pdf", width = 7, height = 5)


#######################Plot faux positifs k#############################
plot(rList, faux_positifsRatelsup1[,4,3],
     type = "l", lwd = 4, col = "red",
     #log = "y",              # axe Y en log10
     xlab = "", ylab = "",   # pas de noms d'axes
     yaxt = "n", xaxt = "n",
     #ylim = c(0, 1e2)) # on redessine les axes manuellement
     ylim = c(0, 20)
)
# Axe Y en log10 avec labels 10^k
# ylim <- range(faux_negatifsRatek[,4,3], na.rm = TRUE)
# ylim[1] <- max(ylim[1], min(faux_negatifsRatek[faux_negatifsRatek[,4,3] > 0,4,3], na.rm = TRUE)) # > 0
#############KL 100############################
lines(rList,faux_positifsRatelsup1[,4,2],lwd = 4,col = "red",lty = "dashed")
lines(rList,faux_positifsRatelsup1[,4,1],lwd = 4,col = "red",lty = "dotted")
#############KL 10############################
lines(rList,faux_positifsRatelsup1[,3,3],lwd = 4,col = "orange")
lines(rList,faux_positifsRatelsup1[,3,2],lwd = 4,col = "orange",lty = "dashed")
lines(rList,faux_positifsRatelsup1[,3,1],lwd = 4,col = "orange",lty = "dotted")
#############KL 1############################
lines(rList,faux_positifsRatelsup1[,2,3],lwd = 4,col = "darkgreen")
lines(rList,faux_positifsRatelsup1[,2,2],lwd =4,col = "darkgreen",lty = "dashed")
lines(rList,faux_positifsRatelsup1[,2,1],lwd = 4,col = "darkgreen",lty = "dotted")
#############KL 0############################
lines(rList,faux_positifsRatelsup1[,1,3],lwd = 4,col = "blue")
lines(rList,faux_positifsRatelsup1[,1,2],lwd = 4,col = "blue",lty = "dashed")
lines(rList,faux_positifsRatelsup1[,1,1],lwd = 4,col = "blue",lty = "dotted")


# log_ticks <- 10^seq(floor(log10(ylim[1])), ceiling(log10(ylim[2])), by = 1)
# axis(2, at = log_ticks,
#      labels = parse(text = paste0("10^", seq(floor(log10(ylim[1])),
#                                              ceiling(log10(ylim[2])), 1))))


# log_ticks <- 10^seq(-2, 2, by = 1)  # ticks de 10^-2 à 10^5
# axis(2, at = log_ticks,
#      labels = parse(text = paste0("10^", -2:2)))
# # Axe X tous les 50
axis(1, at = seq(0, max(rList), by = 5),cex.axis = 2.3)



# Axe X tous les 50

axis(2, at = seq(0, 100, by = 5),las = 1,cex.axis = 2.3)


dev.off()



setwd("~")

pdf("falsepositivesklrho1d10.pdf", width = 7, height = 5)


#######################Plot faux positifs k#############################
plot(rList, faux_positifsRateRho1[,4,3],
     type = "l", lwd = 4, col = "red",
     #log = "y",              # axe Y en log10
     xlab = "", ylab = "",   # pas de noms d'axes
     yaxt = "n", xaxt = "n",
     #ylim = c(0, 1e2)) # on redessine les axes manuellement
     ylim = c(0, 20)
)
# Axe Y en log10 avec labels 10^k
# ylim <- range(faux_negatifsRatek[,4,3], na.rm = TRUE)
# ylim[1] <- max(ylim[1], min(faux_negatifsRatek[faux_negatifsRatek[,4,3] > 0,4,3], na.rm = TRUE)) # > 0
#############KL 100############################
lines(rList,faux_positifsRateRho1[,4,2],lwd = 4,col = "red",lty = "dashed")
lines(rList,faux_positifsRateRho1[,4,1],lwd = 4,col = "red",lty = "dotted")
#############KL 10############################
lines(rList,faux_positifsRateRho1[,3,3],lwd = 4,col = "orange")
lines(rList,faux_positifsRateRho1[,3,2],lwd = 4,col = "orange",lty = "dashed")
lines(rList,faux_positifsRateRho1[,3,1],lwd = 4,col = "orange",lty = "dotted")
#############KL 1############################
lines(rList,faux_positifsRateRho1[,2,3],lwd = 4,col = "darkgreen")
lines(rList,faux_positifsRateRho1[,2,2],lwd =4,col = "darkgreen",lty = "dashed")
lines(rList,faux_positifsRateRho1[,2,1],lwd = 4,col = "darkgreen",lty = "dotted")
#############KL 0############################
lines(rList,faux_positifsRateRho1[,1,3],lwd = 4,col = "blue")
lines(rList,faux_positifsRateRho1[,1,2],lwd = 4,col = "blue",lty = "dashed")
lines(rList,faux_positifsRateRho1[,1,1],lwd = 4,col = "blue",lty = "dotted")


# log_ticks <- 10^seq(floor(log10(ylim[1])), ceiling(log10(ylim[2])), by = 1)
# axis(2, at = log_ticks,
#      labels = parse(text = paste0("10^", seq(floor(log10(ylim[1])),
#                                              ceiling(log10(ylim[2])), 1))))


# log_ticks <- 10^seq(-2, 2, by = 1)  # ticks de 10^-2 à 10^5
# axis(2, at = log_ticks,
#      labels = parse(text = paste0("10^", -2:2)))
# # Axe X tous les 50
axis(1, at = seq(0, max(rList), by = 5),cex.axis = 2.3)



# Axe X tous les 50

axis(2, at = seq(0, 100, by = 5),las = 1,cex.axis = 2.3)


dev.off()


#######################Plot erreurs Sigma l inf 1#############################
plot(rList, erreursSigmalinf1[,4,3],
     type = "l", lwd = 2, col = "red",
     log = "y",              # axe Y en log10
     xlab = "", ylab = "",   # pas de noms d'axes
     yaxt = "n", xaxt = "n",
     ylim = c(1e-2, 1e9)) # on redessine les axes manuellement
# Axe Y en log10 avec labels 10^k
# ylim <- range(erreursSigmak[,4,3], na.rm = TRUE)
# ylim[1] <- max(ylim[1], min(erreursSigmak[erreursSigmak[,4,3] > 0,4,3], na.rm = TRUE)) # > 0
#############KL 100############################
lines(rList,erreursSigmalinf1[,4,2],col = "red",lty = "dashed")
lines(rList,erreursSigmalinf1[,4,1],col = "red",lty = "dotted")
#############KL 10############################
lines(rList,erreursSigmalinf1[,3,3],col = "orange")
lines(rList,erreursSigmalinf1[,3,2],col = "orange",lty = "dashed")
lines(rList,erreursSigmalinf1[,3,1],col = "orange",lty = "dotted")
#############KL 1############################
lines(rList,erreursSigmalinf1[,2,3],col = "green")
lines(rList,erreursSigmalinf1[,2,2],col = "green",lty = "dashed")
lines(rList,erreursSigmalinf1[,2,1],col = "green",lty = "dotted")
#############KL 0############################
lines(rList,erreursSigmalinf1[,1,3],col = "blue")
lines(rList,erreursSigmalinf1[,1,2],col = "blue",lty = "dashed")
lines(rList,erreursSigmalinf1[,1,1],col = "blue",lty = "dotted")



log_ticks <- 10^seq(-2, 9, by = 1)  # ticks de 10^-2 à 10^5
axis(2, at = log_ticks,
     labels = parse(text = paste0("10^", -2:9)))
# Axe X tous les 50
axis(1, at = seq(min(rList), max(rList), by = 5))



#######################Plot faux positifs l <= 1#############################
plot(rList, faux_positifsRatelinf1[,4,3],
     type = "l", lwd = 2, col = "red",
     #log = "y",              # axe Y en log10
     xlab = "", ylab = "",   # pas de noms d'axes
     yaxt = "n", xaxt = "n"
     #ylim = c(0, 1e2)) # on redessine les axes manuellement
     # ylim = c(0, 1e2)
)
# Axe Y en log10 avec labels 10^k
# ylim <- range(faux_negatifsRatek[,4,3], na.rm = TRUE)
# ylim[1] <- max(ylim[1], min(faux_negatifsRatek[faux_negatifsRatek[,4,3] > 0,4,3], na.rm = TRUE)) # > 0
#############KL 100############################
lines(rList,faux_positifsRatelinf1[,4,2],col = "red",lty = "dashed")
lines(rList,faux_positifsRatelinf1[,4,1],col = "red",lty = "dotted")
#############KL 10############################
lines(rList,faux_positifsRatelinf1[,3,3],col = "orange")
lines(rList,faux_positifsRatelinf1[,3,2],col = "orange",lty = "dashed")
lines(rList,faux_positifsRatelinf1[,3,1],col = "orange",lty = "dotted")
#############KL 1############################
lines(rList,faux_positifsRatelinf1[,2,3],col = "green")
lines(rList,faux_positifsRatelinf1[,2,2],col = "green",lty = "dashed")
lines(rList,faux_positifsRatelinf1[,2,1],col = "green",lty = "dotted")
#############KL 0############################
lines(rList,faux_positifsRatelinf1[,1,3],col = "blue")
lines(rList,faux_positifsRatelinf1[,1,2],col = "blue",lty = "dashed")
lines(rList,faux_positifsRatelinf1[,1,1],col = "blue",lty = "dotted")


# log_ticks <- 10^seq(floor(log10(ylim[1])), ceiling(log10(ylim[2])), by = 1)
# axis(2, at = log_ticks,
#      labels = parse(text = paste0("10^", seq(floor(log10(ylim[1])),
#                                              ceiling(log10(ylim[2])), 1))))


# log_ticks <- 10^seq(-2, 2, by = 1)  # ticks de 10^-2 à 10^5
# axis(2, at = log_ticks,
#      labels = parse(text = paste0("10^", -2:2)))
# # Axe X tous les 50
axis(1, at = seq(5, max(rList), by = 5))




axis(2, at = seq(0, 100, by = 5),las = 1)





#######################Plot faux négatifs l <= 1#############################
plot(rList[2:length(rList)], faux_negatifsRatelinf1[(2:length(rList)),4,3],
     type = "l", lwd = 2, col = "red",
     #log = "y",              # axe Y en log10
     xlab = "", ylab = "",   # pas de noms d'axes
     yaxt = "n", xaxt = "n",
     #ylim = c(0, 1e2)) # on redessine les axes manuellement
     ylim = c(0, 1e2)
)
# Axe Y en log10 avec labels 10^k
# ylim <- range(faux_negatifsRatek[,4,3], na.rm = TRUE)
# ylim[1] <- max(ylim[1], min(faux_negatifsRatek[faux_negatifsRatek[,4,3] > 0,4,3], na.rm = TRUE)) # > 0
#############KL 100############################
lines(rList[2:length(rList)],faux_negatifsRatelinf1[(2:length(rList)),4,2],col = "red",lty = "dashed")
lines(rList[2:length(rList)],faux_negatifsRatelinf1[(2:length(rList)),4,1],col = "red",lty = "dotted")
#############KL 10############################
lines(rList[2:length(rList)],faux_negatifsRatelinf1[(2:length(rList)),3,3],col = "orange")
lines(rList[2:length(rList)],faux_negatifsRatelinf1[(2:length(rList)),3,2],col = "orange",lty = "dashed")
lines(rList[2:length(rList)],faux_negatifsRatelinf1[(2:length(rList)),3,1],col = "orange",lty = "dotted")
#############KL 1############################
lines(rList[2:length(rList)],faux_negatifsRatelinf1[(2:length(rList)),2,3],col = "green")
lines(rList[2:length(rList)],faux_negatifsRatelinf1[(2:length(rList)),2,2],col = "green",lty = "dashed")
lines(rList[2:length(rList)],faux_negatifsRatelinf1[(2:length(rList)),2,1],col = "green",lty = "dotted")
#############KL 0############################
lines(rList[2:length(rList)],faux_negatifsRatelinf1[(2:length(rList)),1,3],col = "blue")
lines(rList[2:length(rList)],faux_negatifsRatelinf1[(2:length(rList)),1,2],col = "blue",lty = "dashed")
lines(rList[2:length(rList)],faux_negatifsRatelinf1[(2:length(rList)),1,1],col = "blue",lty = "dotted")


# log_ticks <- 10^seq(floor(log10(ylim[1])), ceiling(log10(ylim[2])), by = 1)
# axis(2, at = log_ticks,
#      labels = parse(text = paste0("10^", seq(floor(log10(ylim[1])),
#                                              ceiling(log10(ylim[2])), 1))))


# log_ticks <- 10^seq(-2, 2, by = 1)  # ticks de 10^-2 à 10^5
# axis(2, at = log_ticks,
#      labels = parse(text = paste0("10^", -2:2)))
# # Axe X tous les 5
axis(1, at = seq(5, max(rList), by = 5))


axis(2, at = seq(0, 100, by = 5),las = 1)





#######################Plot erreurs Sigma l sup 1#############################
plot(rList, erreursSigmalsup1[,4,3],
     type = "l", lwd = 2, col = "red",
     log = "y",              # axe Y en log10
     xlab = "", ylab = "",   # pas de noms d'axes
     yaxt = "n", xaxt = "n",
     ylim = c(1e-2, 1e9)) # on redessine les axes manuellement
# # Axe Y en log10 avec labels 10^k
# ylim <- range(erreursSigmalsup1[,4,1], na.rm = TRUE)
# ylim[1] <- max(ylim[1], min(erreursSigmak[erreursSigmalsup1[,4,1] > 0,4,3], na.rm = TRUE)) # > 0
#############KL 100############################
lines(rList,erreursSigmalsup1[,4,2],col = "red",lty = "dashed")
lines(rList,erreursSigmalsup1[,4,1],col = "red",lty = "dotted")
#############KL 10############################
lines(rList,erreursSigmalsup1[,3,3],col = "orange")
lines(rList,erreursSigmalsup1[,3,2],col = "orange",lty = "dashed")
lines(rList,erreursSigmalsup1[,3,1],col = "orange",lty = "dotted")
#############KL 1############################
lines(rList,erreursSigmalsup1[,2,3],col = "green")
lines(rList,erreursSigmalsup1[,2,2],col = "green",lty = "dashed")
lines(rList,erreursSigmalsup1[,2,1],col = "green",lty = "dotted")
#############KL 0############################
lines(rList,erreursSigmalsup1[,1,3],col = "blue")
lines(rList,erreursSigmalsup1[,1,2],col = "blue",lty = "dashed")
lines(rList,erreursSigmalsup1[,1,1],col = "blue",lty = "dotted")

log_ticks <- 10^seq(-2, 9, by = 1)  # ticks de 10^-2 à 10^5
axis(2, at = log_ticks,
     labels = parse(text = paste0("10^", -2:9)))
# Axe X tous les 50
axis(1, at = seq(min(rList), max(rList), by = 5))



#######################Plot faux négatifs lsup 1#############################
plot(rList[2:length(rList)], faux_negatifsRatelsup1[(2:length(rList)),4,3],
     type = "l", lwd = 2, col = "red",
     #log = "y",              # axe Y en log10
     xlab = "", ylab = "",   # pas de noms d'axes
     yaxt = "n", xaxt = "n",
     #ylim = c(0, 1e2)) # on redessine les axes manuellement
     ylim = c(0, 1e2)
)
# Axe Y en log10 avec labels 10^k
# ylim <- range(faux_negatifsRatek[,4,3], na.rm = TRUE)
# ylim[1] <- max(ylim[1], min(faux_negatifsRatek[faux_negatifsRatek[,4,3] > 0,4,3], na.rm = TRUE)) # > 0
#############KL 100############################
lines(rList[2:length(rList)],faux_negatifsRatelsup1[(2:length(rList)),4,2],col = "red",lty = "dashed")
lines(rList[2:length(rList)],faux_negatifsRatelsup1[(2:length(rList)),4,1],col = "red",lty = "dotted")
#############KL 10############################
lines(rList[2:length(rList)],faux_negatifsRatelsup1[(2:length(rList)),3,3],col = "orange")
lines(rList[2:length(rList)],faux_negatifsRatelsup1[(2:length(rList)),3,2],col = "orange",lty = "dashed")
lines(rList[2:length(rList)],faux_negatifsRatelsup1[(2:length(rList)),3,1],col = "orange",lty = "dotted")
#############KL 1############################
lines(rList[2:length(rList)],faux_negatifsRatelsup1[(2:length(rList)),2,3],col = "green")
lines(rList[2:length(rList)],faux_negatifsRatelsup1[(2:length(rList)),2,2],col = "green",lty = "dashed")
lines(rList[2:length(rList)],faux_negatifsRatelsup1[(2:length(rList)),2,1],col = "green",lty = "dotted")
#############KL 0############################
lines(rList[2:length(rList)],faux_negatifsRatelsup1[(2:length(rList)),1,3],col = "blue")
lines(rList[2:length(rList)],faux_negatifsRatelsup1[(2:length(rList)),1,2],col = "blue",lty = "dashed")
lines(rList[2:length(rList)],faux_negatifsRatelsup1[(2:length(rList)),1,1],col = "blue",lty = "dotted")


# log_ticks <- 10^seq(floor(log10(ylim[1])), ceiling(log10(ylim[2])), by = 1)
# axis(2, at = log_ticks,
#      labels = parse(text = paste0("10^", seq(floor(log10(ylim[1])),
#                                              ceiling(log10(ylim[2])), 1))))


# log_ticks <- 10^seq(-2, 2, by = 1)  # ticks de 10^-2 à 10^5
# axis(2, at = log_ticks,
#      labels = parse(text = paste0("10^", -2:2)))
# # Axe X tous les 5
axis(1, at = seq(5, max(rList), by = 5))


axis(2, at = seq(0, 100, by = 5),las = 1)


#######################Plot faux positifs l >= 1#############################
plot(rList, faux_positifsRatelsup1[,4,3],
     type = "l", lwd = 2, col = "red",
     #log = "y",              # axe Y en log10
     xlab = "", ylab = "",   # pas de noms d'axes
     yaxt = "n", xaxt = "n"
     #ylim = c(0, 1e2)) # on redessine les axes manuellement
     # ylim = c(0, 1e2)
)
# Axe Y en log10 avec labels 10^k
# ylim <- range(faux_negatifsRatek[,4,3], na.rm = TRUE)
# ylim[1] <- max(ylim[1], min(faux_negatifsRatek[faux_negatifsRatek[,4,3] > 0,4,3], na.rm = TRUE)) # > 0
#############KL 100############################
lines(rList,faux_positifsRatelsup1[,4,2],col = "red",lty = "dashed")
lines(rList,faux_positifsRatelsup1[,4,1],col = "red",lty = "dotted")
#############KL 10############################
lines(rList,faux_positifsRatelsup1[,3,3],col = "orange")
lines(rList,faux_positifsRatelsup1[,3,2],col = "orange",lty = "dashed")
lines(rList,faux_positifsRatelsup1[,3,1],col = "orange",lty = "dotted")
#############KL 1############################
lines(rList,faux_positifsRatelsup1[,2,3],col = "green")
lines(rList,faux_positifsRatelsup1[,2,2],col = "green",lty = "dashed")
lines(rList,faux_positifsRatelsup1[,2,1],col = "green",lty = "dotted")
#############KL 0############################
lines(rList,faux_positifsRatelsup1[,1,3],col = "blue")
lines(rList,faux_positifsRatelsup1[,1,2],col = "blue",lty = "dashed")
lines(rList,faux_positifsRatelsup1[,1,1],col = "blue",lty = "dotted")


# log_ticks <- 10^seq(floor(log10(ylim[1])), ceiling(log10(ylim[2])), by = 1)
# axis(2, at = log_ticks,
#      labels = parse(text = paste0("10^", seq(floor(log10(ylim[1])),
#                                              ceiling(log10(ylim[2])), 1))))


# log_ticks <- 10^seq(-2, 2, by = 1)  # ticks de 10^-2 à 10^5
# axis(2, at = log_ticks,
#      labels = parse(text = paste0("10^", -2:2)))
# # Axe X tous les 50
axis(1, at = seq(5, max(rList), by = 5))

axis(2, at = seq(0, 8, by = 1),las = 1)


#######################Plot erreurs rho1 >= 0.3#############################
plot(rList, erreursSigmaRho1[,4,3],
     type = "l", lwd = 4, col = "red",
     log = "y",              # axe Y en log10
     xlab = "", ylab = "",   # pas de noms d'axes
     yaxt = "n", xaxt = "n", # on redessine les axes manuellement
     ylim = c(1e-2, 1e9)) 
# # Axe Y en log10 avec labels 10^k
# ylim <- range(erreursSigmalsup1[,4,1], na.rm = TRUE)
# ylim[1] <- max(ylim[1], min(erreursSigmak[erreursSigmalsup1[,4,1] > 0,4,3], na.rm = TRUE)) # > 0
#############KL 100############################
lines(rList,erreursSigmaRho1[,4,2],lwd = 4,col = "red",lty = "dashed")
lines(rList,erreursSigmaRho1[,4,1],lwd = 4,col = "red",lty = "dotted")
#############KL 10############################
lines(rList,erreursSigmaRho1[,3,3],lwd = 4,col = "orange")
lines(rList,erreursSigmaRho1[,3,2],lwd = 4,col = "orange",lty ="dashed")
lines(rList,erreursSigmaRho1[,3,1],lwd = 4,col = "orange",lty ="dotted")
#############KL 1############################
lines(rList,erreursSigmaRho1[,2,3],lwd = 4,col = "green")
lines(rList,erreursSigmaRho1[,2,2],lwd = 4,col = "green",lty ="dashed")
lines(rList,erreursSigmaRho1[,2,1],lwd = 4,col = "green",lty ="dotted")
#############KL 1############################
lines(rList,erreursSigmaRho1[,1,3],lwd = 4,col = "blue")
lines(rList,erreursSigmaRho1[,1,2],lwd = 4,col = "blue",lty ="dashed")
lines(rList,erreursSigmaRho1[,1,1],lwd = 4,col = "blue",lty ="dotted")


log_ticks <- 10^seq(-2, 9, by = 1)  # ticks de 10^-2 à 10^5
axis(2, at = log_ticks,
     labels = parse(text = paste0("10^", -2:9)),cex.axis = 1.8)
# Axe X tous les 50
axis(1, at = seq(min(rList), max(rList), by = 5),cex.axis = 1.8)

#######################Plot faux positifs rho1 >= 0.3 #############################
plot(rList, faux_positifsRateRho1[,4,3],
     type = "l", lwd = 4, col = "red",
     #log = "y",              # axe Y en log10
     xlab = "", ylab = "",   # pas de noms d'axes
     yaxt = "n", xaxt = "n"
     #ylim = c(0, 1e2)) # on redessine les axes manuellement
     # ylim = c(0, 1e2)
)
# Axe Y en log10 avec labels 10^k
# ylim <- range(faux_negatifsRatek[,4,3], na.rm = TRUE)
# ylim[1] <- max(ylim[1], min(faux_negatifsRatek[faux_negatifsRatek[,4,3] > 0,4,3], na.rm = TRUE)) # > 0
#############KL 100############################
lines(rList,faux_positifsRateRho1[,4,2],col = "red",lwd = 4,lty = "dashed")
lines(rList,faux_positifsRateRho1[,4,1],col = "red",lwd = 4,lty = "dotted")
#############KL 10############################
lines(rList,faux_positifsRateRho1[,3,3],lwd = 4,col = "orange")
lines(rList,faux_positifsRateRho1[,3,2],lwd = 4,col = "orange",lty = "dashed")
lines(rList,faux_positifsRateRho1[,3,1],lwd = 4,col = "orange",lty = "dotted")
#############KL 1############################
lines(rList,faux_positifsRateRho1[,2,3],lwd = 4,col = "green")
lines(rList,faux_positifsRateRho1[,2,2],lwd = 4,col = "green",lty = "dashed")
lines(rList,faux_positifsRateRho1[,2,1],lwd = 4,col = "green",lty = "dotted")
#############KL 0############################
lines(rList,faux_positifsRateRho1[,1,3],lwd = 4,col = "blue")
lines(rList,faux_positifsRateRho1[,1,2],lwd = 4,col = "blue",lty = "dashed")
lines(rList,faux_positifsRateRho1[,1,1],lwd = 4,col = "blue",lty = "dotted")


# log_ticks <- 10^seq(floor(log10(ylim[1])), ceiling(log10(ylim[2])), by = 1)
# axis(2, at = log_ticks,
#      labels = parse(text = paste0("10^", seq(floor(log10(ylim[1])),
#                                              ceiling(log10(ylim[2])), 1))))


# log_ticks <- 10^seq(-2, 2, by = 1)  # ticks de 10^-2 à 10^5
# axis(2, at = log_ticks,
#      labels = parse(text = paste0("10^", -2:2)))
# # Axe X tous les 50
axis(1, at = seq(5, max(rList), by = 5),cex.axis = 1.8)

axis(2, at = seq(0, 100, by = 5),las = 1,cex.axis = 1.8)






#######################Plot faux négatifs rho1 <= 0.3#############################
plot(rList[2:length(rList)], faux_negatifsRateRho1[(2:length(rList)),4,3],
     type = "l", lwd = 2, col = "red",
     #log = "y",              # axe Y en log10
     xlab = "", ylab = "",   # pas de noms d'axes
     yaxt = "n", xaxt = "n",
     #ylim = c(0, 1e2)) # on redessine les axes manuellement
     ylim = c(0, 1e2)
)
# Axe Y en log10 avec labels 10^k
# ylim <- range(faux_negatifsRatek[,4,3], na.rm = TRUE)
# ylim[1] <- max(ylim[1], min(faux_negatifsRatek[faux_negatifsRatek[,4,3] > 0,4,3], na.rm = TRUE)) # > 0
#############KL 100############################
lines(rList[2:length(rList)],faux_negatifsRateRho1[(2:length(rList)),4,2],col = "red",lty = "dashed")
lines(rList[2:length(rList)],faux_negatifsRateRho1[(2:length(rList)),4,1],col = "red",lty = "dotted")
#############KL 10############################
lines(rList[2:length(rList)],faux_negatifsRateRho1[(2:length(rList)),3,3],col = "orange")
lines(rList[2:length(rList)],faux_negatifsRateRho1[(2:length(rList)),3,2],col = "orange",lty = "dashed")
lines(rList[2:length(rList)],faux_negatifsRateRho1[(2:length(rList)),3,1],col = "orange",lty = "dotted")
#############KL 1############################
lines(rList[2:length(rList)],faux_negatifsRateRho1[(2:length(rList)),2,3],col = "green")
lines(rList[2:length(rList)],faux_negatifsRateRho1[(2:length(rList)),2,2],col = "green",lty = "dashed")
lines(rList[2:length(rList)],faux_negatifsRateRho1[(2:length(rList)),2,1],col = "green",lty = "dotted")
#############KL 0############################
lines(rList[2:length(rList)],faux_negatifsRateRho1[(2:length(rList)),1,3],col = "blue")
lines(rList[2:length(rList)],faux_negatifsRateRho1[(2:length(rList)),1,2],col = "blue",lty = "dashed")
lines(rList[2:length(rList)],faux_negatifsRateRho1[(2:length(rList)),1,1],col = "blue",lty = "dotted")


# log_ticks <- 10^seq(floor(log10(ylim[1])), ceiling(log10(ylim[2])), by = 1)
# axis(2, at = log_ticks,
#      labels = parse(text = paste0("10^", seq(floor(log10(ylim[1])),
#                                              ceiling(log10(ylim[2])), 1))))


# log_ticks <- 10^seq(-2, 2, by = 1)  # ticks de 10^-2 à 10^5
# axis(2, at = log_ticks,
#      labels = parse(text = paste0("10^", -2:2)))
# # Axe X tous les 5
axis(1, at = seq(5, max(rList), by = 5))


axis(2, at = seq(0, 100, by = 5),las = 1)




#######################Plot erreurs rho1 <= 0.3#############################
plot(rList, erreursSigmaRho1Neg[,4,3],
     type = "l", lwd = 2, col = "red",
     log = "y",              # axe Y en log10
     xlab = "", ylab = "",   # pas de noms d'axes
     yaxt = "n", xaxt = "n",
     ylim = c(1e-2, 1e9)) # on redessine les axes manuellement
# # Axe Y en log10 avec labels 10^k
# ylim <- range(erreursSigmalsup1[,4,1], na.rm = TRUE)
# ylim[1] <- max(ylim[1], min(erreursSigmak[erreursSigmalsup1[,4,1] > 0,4,3], na.rm = TRUE)) # > 0
#############KL 100############################
lines(rList,erreursSigmaRho1Neg[,4,2],col = "red",lty = "dashed")
lines(rList,erreursSigmaRho1Neg[,4,1],col = "red",lty = "dotted")
#############KL 10############################
lines(rList,erreursSigmaRho1Neg[,3,3],col = "orange")
lines(rList,erreursSigmaRho1Neg[,3,2],col = "orange",lty ="dashed")
lines(rList,erreursSigmaRho1Neg[,3,1],col = "orange",lty ="dotted")
#############KL 1############################
lines(rList,erreursSigmaRho1Neg[,2,3],col = "green")
lines(rList,erreursSigmaRho1Neg[,2,2],col = "green",lty ="dashed")
lines(rList,erreursSigmaRho1Neg[,2,1],col = "green",lty ="dotted")
#############KL 1############################
lines(rList,erreursSigmaRho1Neg[,1,3],col = "blue")
lines(rList,erreursSigmaRho1Neg[,1,2],col = "blue",lty ="dashed")
lines(rList,erreursSigmaRho1Neg[,1,1],col = "blue",lty ="dotted")


log_ticks <- 10^seq(-2, 9, by = 1)  
axis(2, at = log_ticks,
     labels = parse(text = paste0("10^", -2:9)))
# Axe X tous les 50
axis(1, at = seq(min(rList), max(rList), by = 5))



#######################Plot faux négatifs rho1 <= 0.3#############################
plot(rList[2:length(rList)], faux_negatifsRateRho1Neg[(2:length(rList)),4,3],
     type = "l", lwd = 2, col = "red",
     #log = "y",              # axe Y en log10
     xlab = "", ylab = "",   # pas de noms d'axes
     yaxt = "n", xaxt = "n",
     #ylim = c(0, 1e2)) # on redessine les axes manuellement
     ylim = c(0, 1e2)
)
# Axe Y en log10 avec labels 10^k
# ylim <- range(faux_negatifsRatek[,4,3], na.rm = TRUE)
# ylim[1] <- max(ylim[1], min(faux_negatifsRatek[faux_negatifsRatek[,4,3] > 0,4,3], na.rm = TRUE)) # > 0
#############KL 100############################
lines(rList[2:length(rList)],faux_negatifsRateRho1Neg[(2:length(rList)),4,2],col = "red",lty = "dashed")
lines(rList[2:length(rList)],faux_negatifsRateRho1Neg[(2:length(rList)),4,1],col = "red",lty = "dotted")
#############KL 10############################
lines(rList[2:length(rList)],faux_negatifsRateRho1Neg[(2:length(rList)),3,3],col = "orange")
lines(rList[2:length(rList)],faux_negatifsRateRho1Neg[(2:length(rList)),3,2],col = "orange",lty = "dashed")
lines(rList[2:length(rList)],faux_negatifsRateRho1Neg[(2:length(rList)),3,1],col = "orange",lty = "dotted")
#############KL 1############################
lines(rList[2:length(rList)],faux_negatifsRateRho1Neg[(2:length(rList)),2,3],col = "green")
lines(rList[2:length(rList)],faux_negatifsRateRho1Neg[(2:length(rList)),2,2],col = "green",lty = "dashed")
lines(rList[2:length(rList)],faux_negatifsRateRho1Neg[(2:length(rList)),2,1],col = "green",lty = "dotted")
#############KL 0############################
lines(rList[2:length(rList)],faux_negatifsRateRho1Neg[(2:length(rList)),1,3],col = "blue")
lines(rList[2:length(rList)],faux_negatifsRateRho1Neg[(2:length(rList)),1,2],col = "blue",lty = "dashed")
lines(rList[2:length(rList)],faux_negatifsRateRho1Neg[(2:length(rList)),1,1],col = "blue",lty = "dotted")


# log_ticks <- 10^seq(floor(log10(ylim[1])), ceiling(log10(ylim[2])), by = 1)
# axis(2, at = log_ticks,
#      labels = parse(text = paste0("10^", seq(floor(log10(ylim[1])),
#                                              ceiling(log10(ylim[2])), 1))))


# log_ticks <- 10^seq(-2, 2, by = 1)  # ticks de 10^-2 à 10^5
# axis(2, at = log_ticks,
#      labels = parse(text = paste0("10^", -2:2)))
# # Axe X tous les 5
axis(1, at = seq(5, max(rList), by = 5))


axis(2, at = seq(0, 100, by = 5),las = 1)


#######################Plot faux positifs rho1 <= 0.3 #############################
plot(rList, faux_positifsRateRho1Neg[,4,3],
     type = "l", lwd = 2, col = "red",
     #log = "y",              # axe Y en log10
     xlab = "", ylab = "",   # pas de noms d'axes
     yaxt = "n", xaxt = "n"
     #ylim = c(0, 1e2)) # on redessine les axes manuellement
     # ylim = c(0, 1e2)
)
# Axe Y en log10 avec labels 10^k
# ylim <- range(faux_negatifsRatek[,4,3], na.rm = TRUE)
# ylim[1] <- max(ylim[1], min(faux_negatifsRatek[faux_negatifsRatek[,4,3] > 0,4,3], na.rm = TRUE)) # > 0
#############KL 100############################
lines(rList,faux_positifsRateRho1Neg[,4,2],col = "red",lty = "dashed")
lines(rList,faux_positifsRateRho1Neg[,4,1],col = "red",lty = "dotted")
#############KL 10############################
lines(rList,faux_positifsRateRho1Neg[,3,3],col = "orange")
lines(rList,faux_positifsRateRho1Neg[,3,2],col = "orange",lty = "dashed")
lines(rList,faux_positifsRateRho1Neg[,3,1],col = "orange",lty = "dotted")
#############KL 1############################
lines(rList,faux_positifsRateRho1Neg[,2,3],col = "green")
lines(rList,faux_positifsRateRho1Neg[,2,2],col = "green",lty = "dashed")
lines(rList,faux_positifsRateRho1Neg[,2,1],col = "green",lty = "dotted")
#############KL 0############################
lines(rList,faux_positifsRateRho1Neg[,1,3],col = "blue")
lines(rList,faux_positifsRateRho1Neg[,1,2],col = "blue",lty = "dashed")
lines(rList,faux_positifsRateRho1Neg[,1,1],col = "blue",lty = "dotted")


# log_ticks <- 10^seq(floor(log10(ylim[1])), ceiling(log10(ylim[2])), by = 1)
# axis(2, at = log_ticks,
#      labels = parse(text = paste0("10^", seq(floor(log10(ylim[1])),
#                                              ceiling(log10(ylim[2])), 1))))


# log_ticks <- 10^seq(-2, 2, by = 1)  # ticks de 10^-2 à 10^5
# axis(2, at = log_ticks,
#      labels = parse(text = paste0("10^", -2:2)))
# # Axe X tous les 50
axis(1, at = seq(5, max(rList), by = 5))

axis(2, at = seq(0, 100, by = 5),las = 1)


###################################################
#Affiche contamination scenarios
#################################################
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

k =k1val[2];l=l1val[5];rho1 = rho1val[2];

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

file = paste0("contaminScenNearr",rate,".pdf")

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


k =k1val[3];l=l1val[6];rho1 = rho1val[3];

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

file = paste0("contaminScenMedr",rate,".pdf")

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


k =k1val[length(k1val)];l=l1val[length(l1val)];rho1 = rho1val[length(rho1val)];

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

file = paste0("contaminScenFarr",rate,".pdf")

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






######################################################################
#######Affichage des erreurs d'estimation de Sigma
######################################################################

# 
plot(1:n, erreursSigmaFar[,2,3,1],
     type = "l", lwd = 2, 
     xlab = "", ylab = "",   # Pas de label
     yaxt = "n", xaxt = "n", # On masque les axes par défaut
     log = "y",              # Échelle logarithmique Y
     ylim = c(1e-1, 1e10)     # Plage Y adaptée à tes ticks log
)
# 
# # Autres courbes
lines(1:n, erreursSigmaFar[,2,2,1], lty = "dashed")
lines(1:n, erreursSigmaFar[,2,1,1], lty = "dotted")
# 
# # Axe Y logarithmique
log_ticks <- 10^seq(-1, 10, by = 1)
axis(2, at = log_ticks,
     labels = parse(text = paste0("10^", -1:10)),
     las = 1)
x_ticks <- seq(0, n, by = 1000)
axis(1, at = x_ticks,las = 1)




plot(1:n, erreursSigmaFar[,5,3,1],
     type = "l", lwd = 2, 
     xlab = "", ylab = "",   # Pas de label
     yaxt = "n", xaxt = "n", # On masque les axes par défaut
     log = "y",              # Échelle logarithmique Y
     ylim = c(1e-1, 1e10)     # Plage Y adaptée à tes ticks log
)
# 
# # Autres courbes
lines(1:n, erreursSigmaFar[,5,2,1], lty = "dashed")
lines(1:n, erreursSigmaFar[,5,1,1], lty = "dotted")
# 
# # Axe Y logarithmique
log_ticks <- 10^seq(-1, 10, by = 1)
axis(2, at = log_ticks,
     labels = parse(text = paste0("10^", -1:10)),
     las = 1)
x_ticks <- seq(0, n, by = 1000)
axis(1, at = x_ticks,las = 1)



plot(1:n, erreursSigmaFar[,7,3,1],
     type = "l", lwd = 2, 
     xlab = "", ylab = "",   # Pas de label
     yaxt = "n", xaxt = "n", # On masque les axes par défaut
     log = "y",              # Échelle logarithmique Y
     ylim = c(1e-1, 1e10)     # Plage Y adaptée à tes ticks log
)
# 
# # Autres courbes
lines(1:n, erreursSigmaFar[,7,2,1], lty = "dashed")
lines(1:n, erreursSigmaFar[,7,1,1], lty = "dotted")
# 
# # Axe Y logarithmique
log_ticks <- 10^seq(-1, 10, by = 1)
axis(2, at = log_ticks,
     labels = parse(text = paste0("10^", -1:10)),
     las = 1)
x_ticks <- seq(0, n, by = 1000)
axis(1, at = x_ticks,las = 1) 



plot(1:n, erreursSigmaFar[,11,3,1],
     type = "l", lwd = 2, 
     xlab = "", ylab = "",   # Pas de label
     yaxt = "n", xaxt = "n", # On masque les axes par défaut
     log = "y",              # Échelle logarithmique Y
     ylim = c(1e-1, 1e10)     # Plage Y adaptée à tes ticks log
)
# 
# # Autres courbes
lines(1:n, erreursSigmaFar[,11,2,1], lty = "dashed")
lines(1:n, erreursSigmaFar[,11,1,1], lty = "dotted")
# 
# # Axe Y logarithmique
log_ticks <- 10^seq(-1,10, by = 1)
axis(2, at = log_ticks,
     labels = parse(text = paste0("10^", -1:10)),
     las = 1)
x_ticks <- seq(0, n, by = 1000)
axis(1, at = x_ticks,las = 1) 

#######################################################False positives and false negatives########################################



plot(rList[2:length(rList)], faux_negatifsId[2:length(rList),3,1]/((rList[2:length(rList)]/100)*n)*100,
     type = "l", lwd = 2, 
     xlab = "", ylab = "",   # Pas de label
     yaxt = "n", xaxt = "n", # On masque les axes par défaut
     #log = "y",              # Échelle logarithmique Y
     ylim = c(0, 100)     # Plage Y adaptée à tes ticks log
)
# 
# # Autres courbes
lines(rList[2:length(rList)], faux_negatifsId[2:length(rList),2,1]/((rList[2:length(rList)]/100)*n)*100, lty = "dashed")
lines(rList[2:length(rList)], faux_negatifsId[2:length(rList),1,1]/((rList[2:length(rList)]/100)*n)*100, lty = "dotted")
#  
# # Axe Y logarithmique
# log_ticks <- 10^seq(-1, 10, by = 1)
# axis(2, at = log_ticks,
#      labels = parse(text = paste0("10^", -1:10)),
#      las = 1)
x_ticks <- rList[2:length(rList)]
axis(2, at = seq(0, 100, by = 5), las = 1)

axis(1, at = x_ticks,las = 1)

# 
plot(1:n, erreursSigmaFar[,2,3,1],
     type = "l", lwd = 2, 
     xlab = "", ylab = "",   # Pas de label
     yaxt = "n", xaxt = "n", # On masque les axes par défaut
     log = "y",              # Échelle logarithmique Y
     ylim = c(1e-1, 1e10)     # Plage Y adaptée à tes ticks log
)
# 
# # Autres courbes
lines(1:n, erreursSigmaFar[,2,2,1], lty = "dashed")
lines(1:n, erreursSigmaFar[,2,1,1], lty = "dotted")
# 
# # Axe Y logarithmique
log_ticks <- 10^seq(-1, 10, by = 1)
axis(2, at = log_ticks,
     labels = parse(text = paste0("10^", -1:10)),
     las = 1)
x_ticks <- seq(0, n, by = 1000)
axis(1, at = x_ticks,las = 1)


setwd("~")

png("falsepositives_id.png", width = 1800, height = 1200, res = 200)


plot(rList[1:9], moyenne_faux_positifsId[1:9,3]/((1 - rList[1:9]/100)*n)*100,
     type = "l", lwd = 4, 
     col = "red",
     xlab = "", ylab = "",   # Pas de label
     yaxt = "n", xaxt = "n", # On masque les axes par défaut
     #log = "y",              # Échelle logarithmique Y
     ylim = c(0, 15)     # Plage Y adaptée à tes ticks log
)
# 
# # Autres courbes
lines(rList[1:9], moyenne_faux_positifsId[1:9,2]/((1 - rList[1:9]/100)*n)*100, ,lwd = 4,lty = "dashed",col = "blue")
lines(rList[1:9], moyenne_faux_positifsId[1:9,1]/((1 - rList[1:9]/100)*n)*100, lwd = 4,lty = "dotted",col = "darkgreen")
#  
# # Axe Y logarithmique
# log_ticks <- 10^seq(-1, 10, by = 1)
# axis(2, at = log_ticks,
#      labels = parse(text = paste0("10^", -1:10)),
#      las = 1)
x_ticks <- rList[2:length(rList)]
axis(2, at = seq(0, 15, by = 5), las = 1,cex.axis = 2.4)

axis(1, at = x_ticks,las = 1,cex.axis = 2.7)

dev.off()



setwd("~")

png("falsepositives_near.png", width = 1800, height = 1200, res = 200)


plot(rList[1:9], moyenne_faux_positifsNear[1:9,3]/((1 - rList[1:9]/100)*n)*100,
     type = "l", lwd = 4, 
     col = "red",
     xlab = "", ylab = "",   # Pas de label
     yaxt = "n", xaxt = "n", # On masque les axes par défaut
     #log = "y",              # Échelle logarithmique Y
     ylim = c(0, 15)     # Plage Y adaptée à tes ticks log
)
# 
# # Autres courbes
lines(rList[1:9], moyenne_faux_positifsNear[1:9,2]/((1 - rList[1:9]/100)*n)*100, ,lwd = 4,lty = "dashed",col = "blue")
lines(rList[1:9], moyenne_faux_positifsNear[1:9,1]/((1 - rList[1:9]/100)*n)*100, lwd = 4,lty = "dotted",col = "darkgreen")
#  
# # Axe Y logarithmique
# log_ticks <- 10^seq(-1, 10, by = 1)
# axis(2, at = log_ticks,
#      labels = parse(text = paste0("10^", -1:10)),
#      las = 1)
x_ticks <- rList[2:length(rList)]
axis(2, at = seq(0, 15, by = 5), las = 1,cex.axis = 2.4)

axis(1, at = x_ticks,las = 1,cex.axis = 2.7)

dev.off()



setwd("~")

png("falsepositives_med.png", width = 1800, height = 1200, res = 200)


plot(rList[1:9], moyenne_faux_positifsMed[1:9,3]/((1 - rList[1:9]/100)*n)*100,
     type = "l", lwd = 4, 
     col = "red",
     xlab = "", ylab = "",   # Pas de label
     yaxt = "n", xaxt = "n", # On masque les axes par défaut
     #log = "y",              # Échelle logarithmique Y
     ylim = c(0, 15)     # Plage Y adaptée à tes ticks log
)
# 
# # Autres courbes
lines(rList[1:9], moyenne_faux_positifsMed[1:9,2]/((1 - rList[1:9]/100)*n)*100, ,lwd = 4,lty = "dashed",col = "blue")
lines(rList[1:9], moyenne_faux_positifsMed[1:9,1]/((1 - rList[1:9]/100)*n)*100, lwd = 4,lty = "dotted",col = "darkgreen")
#  
# # Axe Y logarithmique
# log_ticks <- 10^seq(-1, 10, by = 1)
# axis(2, at = log_ticks,
#      labels = parse(text = paste0("10^", -1:10)),
#      las = 1)
x_ticks <- rList[2:length(rList)]
axis(2, at = seq(0, 15, by = 5), las = 1,cex.axis = 2.4)

axis(1, at = x_ticks,las = 1,cex.axis = 2.7)

dev.off()


setwd("~")

png("falsepositives_far.png", width = 1800, height = 1200, res = 200)


plot(rList[1:9], moyenne_faux_positifsFar[1:9,3]/((1 - rList[1:9]/100)*n)*100,
     type = "l", lwd = 4, 
     col = "red",
     xlab = "", ylab = "",   # Pas de label
     yaxt = "n", xaxt = "n", # On masque les axes par défaut
     #log = "y",              # Échelle logarithmique Y
     ylim = c(0, 15)     # Plage Y adaptée à tes ticks log
)
# 
# # Autres courbes
lines(rList[1:9], moyenne_faux_positifsFar[1:9,2]/((1 - rList[1:9]/100)*n)*100, ,lwd = 4,lty = "dashed",col = "blue")
lines(rList[1:9], moyenne_faux_positifsFar[1:9,1]/((1 - rList[1:9]/100)*n)*100, lwd = 4,lty = "dotted",col = "darkgreen")
#  
# # Axe Y logarithmique
# log_ticks <- 10^seq(-1, 10, by = 1)
# axis(2, at = log_ticks,
#      labels = parse(text = paste0("10^", -1:10)),
#      las = 1)
x_ticks <- rList[2:length(rList)]
axis(2, at = seq(0, 15, by = 5), las = 1,cex.axis = 2.4)

axis(1, at = x_ticks,las = 1,cex.axis = 2.7)

dev.off()






##################################Cumulative outlier detection###################################

########################Cumulative outlier detection###################


resFar5Naive = cumulativeOutlierDetection(labelsVraisFar[,2],majority_vote_Far[,2,1],5)
resFar5Online = cumulativeOutlierDetection(labelsVraisFar[,2],majority_vote_Far[,2,2],5)
resFar5Streaming = cumulativeOutlierDetection(labelsVraisFar[,2],majority_vote_Far[,2,3],5)

# Données supposées : n correspond à la taille de resFar5Naive$taux_outliers_detectes_vrais
x_vals <- 1:n

setwd("~")

file = paste0("cumulativeOutlierDetFarr5",".png")
png(file, width = 1800, height = 1200, res = 200)
# --- Plot principal (échelle log sur X)
plot(x_vals, resFar5Naive$taux_outliers_detectes_vrais,
     type = "l", lwd = 4, lty = "dotted",
     xlab = "", ylab = "",
     yaxt = "n", xaxt = "n",
     ylim = c(0, 100),
     col = "darkgreen",
     log = "x"        # <=== Échelle logarithmique sur l’axe X
)

# --- Autres lignes
lines(x_vals, resFar5Online$taux_outliers_detectes_vrais,
      lwd = 4, col = adjustcolor("blue",alpha.f = 0.5), lty = "dashed")
lines(x_vals, resFar5Streaming$taux_outliers_detectes_vrais,
      lwd = 4, col = adjustcolor("red",alpha.f = 0.8), lty = "solid")

# Ligne invisible (orange)
lines(x_vals, resFar5Naive$taux_outliers_vrais,
      lwd = 4, col = adjustcolor("orange", alpha.f = 0), type = "l")

# --- Axes manuels

# Axe X (logarithmique)
log_ticks_x <- 10^seq(0, floor(log10(n)), by = 1)   # 1, 10, 100, 1000, ...
axis(1, at = log_ticks_x,
     labels = parse(text = paste0("10^", seq(0, floor(log10(n))))),
     las = 1, cex.axis = 2)

# Axe Y (linéaire)
axis(2, at = seq(0, 100, by = 10), las = 1, cex.axis = 2.5)
dev.off() 

setwd("~")

file = paste0("cumulativeOutlierDetFarr20",".png")
png(file, width = 1800, height = 1200, res = 200)
resFar20Naive = cumulativeOutlierDetection(labelsVraisFar[,5],majority_vote_Far[,5,1],20)
resFar20Online = cumulativeOutlierDetection(labelsVraisFar[,5],majority_vote_Far[,5,2],20)
resFar20Streaming = cumulativeOutlierDetection(labelsVraisFar[,5],majority_vote_Far[,5,3],20)
x_vals <- 1:n

# --- Plot principal (échelle log sur X)
plot(x_vals, resFar20Naive$taux_outliers_detectes_vrais,
     type = "l", lwd = 4, lty = "dotted",
     xlab = "", ylab = "",
     yaxt = "n", xaxt = "n",
     ylim = c(0, 100),
     col = "darkgreen",
     log = "x"        # <=== Échelle logarithmique sur l’axe X
)

# --- Autres lignes
lines(x_vals, resFar20Online$taux_outliers_detectes_vrais,
      lwd = 4, col = adjustcolor("blue",alpha.f = 0.5), lty = "dashed")
lines(x_vals, resFar20Streaming$taux_outliers_detectes_vrais,
      lwd = 4, col = adjustcolor("red",alpha.f = 0.8), lty = "solid")

# Ligne invisible (orange)
lines(x_vals, resFar20Naive$taux_outliers_vrais,
      lwd = 4, col = adjustcolor("orange", alpha.f = 0), type = "l")

# --- Axes manuels

# Axe X (logarithmique)
log_ticks_x <- 10^seq(0, floor(log10(n)), by = 1)   # 1, 10, 100, 1000, ...
axis(1, at = log_ticks_x,
     labels = parse(text = paste0("10^", seq(0, floor(log10(n))))),
     las = 1, cex.axis = 2)

# Axe Y (linéaire)
axis(2, at = seq(0, 100, by = 10), las = 1, cex.axis = 2.5)

dev.off()

setwd("~")

file = paste0("cumulativeOutlierDetFarr30",".png")
png(file, width = 1800, height = 1200, res = 200)

resFar30Naive = cumulativeOutlierDetection(labelsVraisFar[,7],majority_vote_Far[,7,1],30)
resFar30Online = cumulativeOutlierDetection(labelsVraisFar[,7],majority_vote_Far[,7,2],30)
resFar30Streaming = cumulativeOutlierDetection(labelsVraisFar[,7],majority_vote_Far[,7,3],30)

# --- Plot principal (échelle log sur X)
plot(x_vals, resFar30Naive$taux_outliers_detectes_vrais,
     type = "l", lwd = 4, lty = "dotted",
     xlab = "", ylab = "",
     yaxt = "n", xaxt = "n",
     ylim = c(0, 100),
     col = "darkgreen",
     log = "x"        # <=== Échelle logarithmique sur l’axe X
)

# --- Autres lignes
lines(x_vals, resFar30Online$taux_outliers_detectes_vrais,
      lwd = 4, col = adjustcolor("blue",alpha.f = 0.5), lty = "dashed")
lines(x_vals, resFar30Streaming$taux_outliers_detectes_vrais,
      lwd = 4, col = adjustcolor("red",alpha.f = 0.8), lty = "solid")

# Ligne invisible (orange)
lines(x_vals, resFar30Naive$taux_outliers_vrais,
      lwd = 4, col = adjustcolor("orange", alpha.f = 0.7), type = "l")

# --- Axes manuels

# Axe X (logarithmique)
log_ticks_x <- 10^seq(0, floor(log10(n)), by = 1)   # 1, 10, 100, 1000, ...
axis(1, at = log_ticks_x,
     labels = parse(text = paste0("10^", seq(0, floor(log10(n))))),
     las = 1, cex.axis = 2)

# Axe Y (linéaire)
axis(2, at = seq(0, 100, by = 10), las = 1, cex.axis = 2.5)

dev.off()

resFar40Naive = cumulativeOutlierDetection(labelsVraisFar[,11],majority_vote_Far[,11,1],50)
resFar40Online = cumulativeOutlierDetection(labelsVraisFar[,11],majority_vote_Far[,11,2],50)
resFar40Streaming = cumulativeOutlierDetection(labelsVraisFar[,11],majority_vote_Far[,11,3],50)

setwd("~")

file = paste0("cumulativeOutlierDetFarr40",".png")
png(file, width = 1800, height = 1200, res = 200)


# --- Plot principal (échelle log sur X)



plot(x_vals, resFar40Naive$taux_outliers_detectes_vrais,
     type = "l", lwd = 4, lty = "dotted",
     xlab = "", ylab = "",
     yaxt = "n", xaxt = "n",
     ylim = c(0, 100),
     col = "darkgreen",
     log = "x"        # <=== Échelle logarithmique sur l’axe X
)

# --- Autres lignes
lines(x_vals, resFar40Online$taux_outliers_detectes_vrais,
      lwd = 4, col = adjustcolor("blue",alpha.f = 0.5), lty = "dashed")
lines(x_vals, resFar40Streaming$taux_outliers_detectes_vrais,
      lwd = 4, col = adjustcolor("red",alpha.f = 0.8), lty = "solid")

# Ligne invisible (orange)
lines(x_vals, resFar40Naive$taux_outliers_vrais,
      lwd = 4, col = adjustcolor("orange", alpha.f = 0.7), type = "l")

# --- Axes manuels

# Axe X (logarithmique)
log_ticks_x <- 10^seq(0, floor(log10(n)), by = 1)   # 1, 10, 100, 1000, ...
axis(1, at = log_ticks_x,
     labels = parse(text = paste0("10^", seq(0, floor(log10(n))))),
     las = 1, cex.axis = 2)

# Axe Y (linéaire)
axis(2, at = seq(0, 100, by = 10), las = 1, cex.axis = 2.5)

dev.off()


# ##############################################
# Temps calculs
# ################# #############################

temps_calculTout = res$temps_calcul
  

# Méthodes et indices souhaités
methodes <- c("sampleCovOnline", "samplecovTrimmed", "sampleCovOffline", "comedianeOffline",
              "comedianeOfflineShrinkage", "OGK", "FASTMCD", "offline", "online", "streaming")

# Indices à garder : 1, 2, 6 à 10
indices_gardes <- c(1, 6, 7, 8, 9, 10)

# Supposons que taux_index est défini
taux_index <- 3  # par exemple

# Extraction des données pour ce taux et méthodes sélectionnées
temps_sel <- temps_calculTout[taux_index, indices_gardes, ]  # dims : méthodes sélectionnées x runs

# Transformation en data frame long
df_temps <- melt(temps_sel, varnames = c("MethodeIndex", "Run"), value.name = "Temps")

# Remplacement par les noms des méthodes sélectionnées
df_temps$Methode <- factor(df_temps$MethodeIndex, 
                           levels = 1:length(indices_gardes), 
                           labels = methodes[indices_gardes])

# Plot boxplot
ggplot(df_temps, aes(x = Methode, y = Temps)) +
  geom_boxplot(fill = "lightblue", outlier.color = "red", outlier.shape = 1) +
  labs(
    title = "Boxplot of computation times",
    x = "Method",
    y = "Time (seconds)"
  ) +
  scale_y_log10() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

######Outliers labels######################

outliers_labelsTout = resMoyenne$outliersLabelsRec
outliers_labelsTout = res1runFarScenario$outliersLabelsRec
#outliers_labelsTout = res1runNearScenario$outliersLabelsRec
outliers_labelsTout = res100runNearesScenario$outliersLabelsRec
#outliers_labelsTout = res1run$outliersLabelsRec[,,,1]

# Fonction pour vote majoritaire
majority_vote <- function(x) {
  ux <- unique(x)
  counts <- tabulate(match(x, ux))
  max_count <- max(counts)
  major_values <- ux[counts == max_count]
  
  if (length(major_values) == 1) {
    return(major_values)
  } else {
    return(1)
  }
}
# Appliquer vote majoritaire sur la 4e dimension
outliers_majority <- apply(outliers_labelsTout, c(1, 2, 3), majority_vote)

#outliers_majority = res$outliersLabelsRec[,,,1]

cumulativeOutlierDetection <- function(labelsVrais,outlier_labels , pourcentage) {
  total_points <- length(labelsVrais)
  
  
  total_outliers_theoriques <- pourcentage / 100 * total_points
  
  nb_outliers_detectes <- 0
  
  nb_outliers_detectes_vrais <- 0
  
    nb_outliers_vrais <- 0
  
  taux_outliers_detectes <- rep(100,total_points)
  
  taux_outliers_vrais <- rep(100,total_points)
  
  taux_outliers_detectes_vrais = rep(100,total_points)
  
  for (i in 1:total_points) {
    
    if (labelsVrais[i] == 1 & as.numeric(outlier_labels[i]) == 1) {
      nb_outliers_detectes <- nb_outliers_detectes + 1
      nb_outliers_vrais <- nb_outliers_vrais + 1
      nb_outliers_detectes_vrais = nb_outliers_detectes_vrais + 1
      taux_outliers_detectes[i] = nb_outliers_detectes/nb_outliers_vrais*100
      taux_outliers_detectes_vrais[i] = nb_outliers_detectes_vrais/nb_outliers_vrais*100
      
          }
    if (labelsVrais[i] == 0 & as.numeric(outlier_labels[i]) == 1) {
      nb_outliers_detectes <- nb_outliers_detectes + 1
      if(nb_outliers_vrais != 0){
        taux_outliers_detectes[i] = nb_outliers_detectes/nb_outliers_vrais*100
        taux_outliers_detectes_vrais[i] = nb_outliers_detectes_vrais/nb_outliers_vrais*100}
      
    }
    
    
    if (labelsVrais[i] == 0 & as.numeric(outlier_labels[i]) == 0) {
      #nb_outliers_detectes <- nb_outliers_detectes + 1
      if(nb_outliers_vrais != 0){
        taux_outliers_detectes[i] = nb_outliers_detectes/nb_outliers_vrais*100
        taux_outliers_detectes_vrais[i] = nb_outliers_detectes_vrais/nb_outliers_vrais*100}
      
    }
    
    
    if (labelsVrais[i] == 1 & as.numeric(outlier_labels[i]) == 0) {
      #nb_outliers_detectesCov <- nb_outliers_detectesCov + 1
      nb_outliers_vrais <- nb_outliers_vrais + 1
      #nb_outliers_detectes_vraisCov = nb_outliers_detectes_vraisCov + 1
      taux_outliers_detectes[i] = nb_outliers_detectes/nb_outliers_vrais*100
      taux_outliers_detectes_vrais[i] = nb_outliers_detectes_vrais/nb_outliers_vrais*100
    }
    
    # 
    # if (labelsVrais[i] == 1 & as.numeric(outlier_labelsOnline[i]) == 1) {
    #   nb_outliers_detectesOnl <- nb_outliers_detectesOnl + 1
    #   nb_outliers_vrais <- nb_outliers_vrais + 1
    #   nb_outliers_detectes_vraisOnl = nb_outliers_detectes_vraisOnl + 1
    #   taux_outliers_detectesOnl[i] = nb_outliers_detectesOnl/nb_outliers_vrais*100
    #   taux_outliers_detectes_vraisOnl[i] = nb_outliers_detectes_vraisOnl/nb_outliers_vrais*100
    # }
    # 
    # if (labelsVrais[i] == 1 & as.numeric(outliers_labelsStrm[i]) == 1) {
    #   nb_outliers_detectesStrm <- nb_outliers_detectesStrm + 1
    #   nb_outliers_vrais <- nb_outliers_vrais + 1
    #   nb_outliers_detectes_vraisStrm = nb_outliers_detectes_vraisStrm + 1
    #   taux_outliers_detectesStrm[i] = nb_outliers_detectesStrm/nb_outliers_vrais*100
    #   taux_outliers_detectes_vraisStrm[i] = nb_outliers_detectes_vraisStrm/nb_outliers_vrais*100
    # }
    #   #taux_outliers_vrais[i] <- nb_outliers_detectes_vrais
    #     if (labelsVrais[i] == 1 && as.numeric(outlier_labelsCovOnline[i]) == 0) {
    #   #nb_outliers_detectes <- nb_outliers_detectes + 1
    #   nb_outliers_vrais = nb_outliers_vrais + 1
    #   
    #     taux_outliers_detectesCov[i] = nb_outliers_detectesCov/nb_outliers_vrais*100
    #     taux_outliers_detectes_vraisCov[i] = nb_outliers_detectes_vraisCov/nb_outliers_vrais*100
    #   
    # }
    # 
    #print(paste0("taux outliers vrais detectes Cov ",taux_outliers_detectes_vraisCov))
    #else {taux_outliers_detectes[i] <-  100
    #taux_outliers_detectes[i] <- nb_outliers_detectes
  }
    
    
    
    
    #taux_outliers_vrais[i] <- 100
    
    # print(paste("nb_outliers_detectes_vrais ",nb_outliers_detectes_vrais ))
    # print(paste("nb_outliers_vrais ",nb_outliers_vrais ))
    # print(paste("taux_outliers_vrais[i] ",taux_outliers_vrais[i] ))
    # print(paste("taux_outliers_vrais[i] ",taux_outliers_vrais[i] ))
    # print(taux_outliers_vrais[1:10])
  return(list(taux_outliers_vrais = taux_outliers_vrais,taux_outliers_detectes = taux_outliers_detectes,taux_outliers_detectes_vrais = taux_outliers_detectes_vrais,nb_outliers_detectes_vrais = nb_outliers_detectes_vrais,nb_outliers_vrais = nb_outliers_vrais))
  } 
  
  
#outliers_majority = outliers_labelsTout[,,1]
# res1run$faux_negatifsRec[9900:1e4]
# table(res1run$labelsVraisRec[,5],res1run$outliersLabelsRec[,5,9,1])
#dim(res10run$outliersLabelsRec)

########Graphs

# pCumOutDetRateFarcOnl5 = cumulativeOutlierDetection(labelsVraisFar[,2],outliersLabelsFar[,2,3],5)
# 
# taux_outliers_detectesCov5 = pCumOutDetRateFarcOnl5$taux_outliers_detectes
# 
# taux_outliers_detectesVraisCov5 = pCumOutDetRateFarcOnl5$taux_outliers_detectes_vrais
# 
# taux_outliers_detectesOnl5 = pCumOutDetRateFarcOnl5$taux_outliers_detectes
# 
# taux_outliers_detectesVraisOnl5 = pCumOutDetRateFarcOnl5$taux_outliers_detectes_vrais
# 
# taux_outliers_detectesStrm5 = pCumOutDetRateFarcOnl5$taux_outliers_detectes
# 
# taux_outliers_detectesVraisStrm5 = pCumOutDetRateFarcOnl5$taux_outliers_detectes_vrais
# 
# 
# taux_outliers_detectesVraisCov5 = pCumOutDetRateFarcOnl5$taux_outliers_detectes_vrais
# 
# 
# pCumOutDetRateNearScOnl10 = cumulativeOutlierDetection(labelsVraisFar[,3],outliersLabelsFar[,3,],10,"")
# 
# pCumOutDetRateNearScOnl20 = cumulativeOutlierDetection(labelsVraisFar[,5],outliersLabelsFar[,5,],20,"")
# 
# pCumOutDetRateNearScOnl35 = cumulativeOutlierDetection(labelsVraisFar[,7],outliersLabelsFar[,7,],35,"")

plot_cumulative_outlier = function(labelsVrais,outliers_labels,rate){

  cumoutDetCov = cumulativeOutlierDetection(labelsVrais,outliers_labels[,1],rate)
  
  taux_outliers_detectesCov = cumoutDetCov$taux_outliers_detectes
  taux_outliers_detectesVraisCov = cumoutDetCov$taux_outliers_detectes_vrais
  
  
  cumoutDetOnl = cumulativeOutlierDetection(labelsVrais,outliers_labels[,2],rate)
  
  taux_outliers_detectesOnl = cumoutDetOnl$taux_outliers_detectes
  taux_outliers_detectesVraisOnl = cumoutDetOnl$taux_outliers_detectes_vrais
  
  
  cumoutDetStrm = cumulativeOutlierDetection(labelsVrais,outliers_labels[,3],rate)
  
  taux_outliers_detectesStrm = cumoutDetStrm$taux_outliers_detectes
  taux_outliers_detectesVraisStrm = cumoutDetStrm$taux_outliers_detectes_vrais
# Préparer les données
df <- data.frame(
  index = 1:total_points,
  True_outliers = taux_outliers_vrais,
  Detected_Cov = taux_outliers_detectesCov,
  TP_Cov = taux_outliers_detectesVraisCov,
  Detected_Onl = taux_outliers_detectesOnl,
  TP_Onl = taux_outliers_detectesVraisOnl,
  Detected_Strm = taux_outliers_detectesStrm,
  TP_Strm = taux_outliers_detectesVraisStrm
)

# Construction du graphique
p <- ggplot(df, aes(x = index)) +
  # Courbes Cov (lignes pleines)
  geom_line(aes(y = Detected_Cov, color = "Detected rate", linetype = "Cov"), size = 0.5) +
  geom_line(aes(y = TP_Cov, color = "True positives rate", linetype = "Cov"), size = 0.5) +
  # Courbes Onl (lignes pointillées)
  geom_line(aes(y = Detected_Onl, color = "Detected rate", linetype = "Onl"), size = 0.5) +
  geom_line(aes(y = TP_Onl, color = "True positives rate", linetype = "Onl"), size = 0.5) +
  # Courbes Strm (lignes en pointillés courts)
  geom_line(aes(y = Detected_Strm, color = "Detected rate", linetype = "Strm"), size = 0.5) +
  geom_line(aes(y = TP_Strm, color = "True positives rate", linetype = "Strm"), size = 0.5) +
  # True outliers (une seule courbe, donc pas besoin de linetype)
  geom_line(aes(y = True_outliers, color = "True outliers rate"), size = 0.5) +
  # Couleurs partagées
  scale_color_manual(values = c(
    "True outliers rate" = "red",
    "Detected rate" = "orange",
    "True positives rate" = "purple"
  )) +
  # Légende des types de ligne
  scale_linetype_manual(values = c(
    "Cov" = "dotted",
    "Onl" = "dashed",
    "Strm" = "solid"
  )) +
  ylim(c(0, 250)) +
  labs(
    title = paste(rate, "% of outliers"),
    x = "",
    y = "",
    color = "Metric",
    linetype = "Method"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Affichage
return(p)
}

pfar5 = plot_cumulative_outlier(labelsVraisFar[,2],outliersLabelsFar[,2,],5 )

pfar20 = plot_cumulative_outlier(labelsVraisFar[,5],outliersLabelsFar[,5,],20 )

pfar30 = plot_cumulative_outlier(labelsVraisFar[,7],outliersLabelsFar[,7,],30 )


pfar40 = plot_cumulative_outlier(labelsVraisFar[,9],outliersLabelsFar[,9,],40 )
library(patchwork)

combined_plot <- (
  (pfar5 | pfar20) /
    (pfar30 | pfar40)
) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

# Ajouter un titre global
final_plot <- combined_plot + 
  plot_annotation(
    title = "Cumulative Outlier Detection F_1 : (k,l,rho1) = (8.59,32,0.975)",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )
  )

# Afficher
print(final_plot)

pNear5 = plot_cumulative_outlier(labelsVraisNear[,2],outliersLabelsNear[,2,],5 )

pNear20 = plot_cumulative_outlier(labelsVraisNear[,5],outliersLabelsNear[,5,],20 )

pNear30 = plot_cumulative_outlier(labelsVraisNear[,7],outliersLabelsNear[,7,],30 )


pNear40 = plot_cumulative_outlier(labelsVraisNear[,9],outliersLabelsNear[,9,],40 )

combined_plot <- (
  (pNear5 | pNear20) /
    (pNear30 | pNear40)
) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

# Ajouter un titre global
final_plot <- combined_plot + 
  plot_annotation(
    title = "Cumulative Outlier Detection F_1 : (k,l,rho1) = (0.86,0.56,0.6)",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )
  )

# Afficher
print(final_plot)


# library(ggplot2)
# library(cowplot)
# 
# # Couleurs des lignes
# rate_colors <- c(
#   "True and false outliers rate" = "orange",
#   "True outliers rate" = "red",
#   
#   "True positive rates" = "purple"
# )
# 
# # Fonction de nettoyage : enlever la légende et le titre de l'axe Y
# clean_plot <- function(p) {
#   p +
#     theme(
#       legend.position = "none",
#       axis.title.y = element_blank()  # Enlever uniquement le titre
#     )
# }
# 
# # Appliquer à chaque graphe
# p1 <- clean_plot(pCumOutDetRateFarcOnl5[[1]])
# p2 <- clean_plot(pCumOutDetRateNearScOnl10[[1]])
# p3 <- clean_plot(pCumOutDetRateNearScOnl20[[1]])
# p4 <- clean_plot(pCumOutDetRateNearScOnl35[[1]])
# 
# # Grille 2x2 de graphes
# main_grid <- plot_grid(
#   p1, p2,
#   p3, p4,
#   ncol = 2,
#   align = "hv"
# )
# 
# # Légende personnalisée
# legend_plot <- ggplot() +
#   annotate("point", x = 1, y = 1, color = rate_colors[1], size = 3) +
#   annotate("point", x = 2, y = 1, color = rate_colors[2], size = 3) +
#   annotate("point", x = 3, y = 1, color = rate_colors[3], size = 3) +
#   annotate("text", x = 1.2, y = 1, label = names(rate_colors)[1], hjust = 0, size = 4) +
#   annotate("text", x = 2.2, y = 1, label = names(rate_colors)[2], hjust = 0, size = 4) +
#   annotate("text", x = 3.2, y = 1, label = names(rate_colors)[3], hjust = 0, size = 4) +
#   theme_void() +
#   xlim(0.5, 4)
# 
# # Combiner graphiques + légende
# body_with_legend <- plot_grid(
#   main_grid,
#   legend_plot,
#   ncol = 1,
#   rel_heights = c(10, 1)
# )
# 
# # Titre principal
# final_plot <- plot_grid(
#   ggdraw() + draw_label("(k,l,rho1) = (8.59,32,0.975) online outlier detection", fontface = "bold", size = 14, hjust = 0.5),
#   body_with_legend,
#   ncol = 1,
#   rel_heights = c(0.08, 0.92)
# )
# 
# # Affichage
# print(final_plot)

###################
#Plot distances and outliers
######################
# plot_outliers_distances <- function(dist_sample, dist_online, dist_streaming, labels_vrais, rate, method = "", ylab = "", dim = 100) {
#   n <- length(dist_sample)
#   
#   # Check longueur cohérente
#   if (!(length(dist_online) == n && length(dist_streaming) == n && length(labels_vrais) == n)) {
#     stop("Tous les vecteurs doivent être de même longueur.")
#   }
#   
#   # Préparer les données
#   df <- data.frame(
#     index = rep(1:n, times = 3),
#     distance = c(dist_sample, dist_online, dist_streaming),
#     Method = factor(rep(c("Sample", "Online", "Streaming"), each = n)),
#     Label = factor(rep(labels_vrais, times = 3))
#   )
#   
#   # Tracer
#   p <- ggplot(df, aes(x = index, y = distance, color = Label, shape = Method)) +
#     geom_point(size = 2, alpha = 0.8) +
#     geom_hline(
#       yintercept = qchisq(0.95, df = dim),
#       linetype = "dashed",
#       color = "black"
#     ) +
#     scale_color_manual(
#       values = c("0" = "blue", "1" = "red"),
#       labels = c("Inliers", "True outliers"),
#       name = "Point type"
#     ) +
#     scale_shape_manual(
#       values = c("Sample" = 15, "Online" = 17, "Streaming" = 8),
#       name = "Method"
#     ) +
#     labs(
#       title = "",
#       subtitle = paste0(rate, "% outliers"),
#       x = "",
#       y = ylab
#     ) +
#     ylim(0, 4000) +
#     theme_minimal() +
#     theme(legend.position = "bottom")
#   
#   return(p)
# }
# # 1. Fonction modifiée pour les plots (sans légende)
plot_outliers_distances_one_method <- function(distances, labels_vrais, rate, ylab = "", cutoff) {
 p =  ggplot(data.frame(index = seq_along(distances),
                    distance = distances,
                    label = factor(labels_vrais)),
         aes(x = index, y = distance)) +
    geom_point(aes(color = label), size = 1.5, show.legend = FALSE) +
    geom_hline(yintercept = cutoff, linetype = "dashed", color = "black", show.legend = FALSE) +
    labs(subtitle = paste0(rate, "% of outliers"), x = "", y = ylab) +
    scale_color_manual(values = c("0" = "blue", "1" = "red")) +
    ylim(c(0, 400)) +
    theme_minimal()
return(p)
 }

#k = 0.86;l=0.56;rho1 = 0.6
#k = 6.57;l = 19.02;rho1 = 0.845
#k = 0.66;l=0.82;rho1 = 0.415
#k = 8.59;l=32;rho1 = 0.975
#method = "streaming"
contamin_rate = c(0,10,20,30,40,50)
sim = 1
plot_list <- list()  # Liste pour stocker les plots

########plot outlier distance One method#################################

for (r in contamin_rate){
  print(paste0("r = ",r))
  dataFile <- paste0('SimData-d', d, '-n', n, '-k', k, '-l', l, '-rho', rho1,'-r',r , '-sim', sim,".RData")
  setwd(simDir)
  contParam = ParmsF1(m1, k, l, rho1)
  data = genererEchantillon(n,n,mu1 = mu0,mu2 = contParam$mu1,Sigma1 = Sigma0,Sigma2 = contParam$Sigma1,r)
  save(data,file  = dataFile)
  #load(dataFile)
  print(paste0("nb_outliers = " ,sum(data$labelsVrais == 1)))
  setwd(resDir)
  resUsStreaming= StreamingOutlierDetection(data$Z,batch = sqrt(ncol(data$Z)),cutoff = 1.38*qchisq(0.95,df = d))
  fitFile = paste0('FitParms-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
  #load(fitFile)
  #resUsStreaming = StreamingOutlierDetection(data$Z,batch = ncol(data$Z))
  if(r == contamin_rate[1] | r == contamin_rate[3]| r == contamin_rate[5]){
    ylab = "Distances"
  }
  
  
  
  
  p = plot_outliers_distances_one_method(distances = resUsStreaming$distances,labels_vrais = data$labelsVrais,ylab = ylab,rate = r,cutoff = 1.38*qchisq(0.95,df = d))
  ylab = ""
  plot_list[[length(plot_list) + 1]] <- p
  print(p)
}
create_standalone_legend <- function() {
  # Création d'un mini dataframe pour la légende
  legend_data <- data.frame(
    x = c(1, 2),
    y = c(1, 1),
    type = factor(c("Inliers", "Outliers"), levels = c("Inliers", "Outliers")),
    cutoff = "Cutoff"
  )
  
  # Construction du plot de légende
  ggplot(legend_data) +
    geom_point(aes(x, y, color = type), size = 4) +
    geom_hline(aes(yintercept = 0.5, linetype = cutoff), color = "black", linewidth = 0.8) +
    scale_color_manual(
      name = "",
      values = c("blue", "red"),
      labels = c("Inliers", "Outliers")
    ) +
    scale_linetype_manual(
      name = "",
      values = "dashed"
    ) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11),
      legend.key.width = unit(1.5, "cm"),
      legend.spacing.x = unit(0.5, "cm")
    )
  
}
  
  # 4. Création de la légende
  legend_plot <- create_standalone_legend()
  legend_grob <- get_legend(legend_plot)

  # 5. Assemblage final
  final_plot <- wrap_plots(plot_list, ncol = 2) /
    legend_grob +
    plot_layout(heights = c(10, 1))

  # Ajout du titre
  final_plot <- final_plot +
    plot_annotation(
      title = paste0("k =", k,    " l = " ,l," rho1 = ",    rho1, " d = ",    d ,    " Streaming method"),
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )

  # Affichage
  print(final_plot)

#   
#   main_title <- "k = 6.57    l = 19.02    rho1 = 0.845    d = 100    Streaming method"
#   # 2. Légende graphique (points + ligne)
#   legend_grob <- grobTree(
#     pointsGrob(x = unit(0.1, "npc"), y = unit(0.6, "npc"), pch = 16, size = unit(5, "pt"), gp = gpar(col = "blue")),
#     textGrob("Inliers", x = unit(0.15, "npc"), y = unit(0.6, "npc"), just = "left", gp = gpar(fontsize = 11)),
#     
#     pointsGrob(x = unit(0.4, "npc"), y = unit(0.6, "npc"), pch = 16, size = unit(5, "pt"), gp = gpar(col = "red")),
#     textGrob("True Outliers", x = unit(0.45, "npc"), y = unit(0.6, "npc"), just = "left", gp = gpar(fontsize = 11)),
#     
#     linesGrob(x = unit(c(0.75, 0.8), "npc"), y = unit(c(0.6, 0.6), "npc"),
#               gp = gpar(lty = "dashed", lwd = 1.5, col = "black")),
#     textGrob("Cutoff", x = unit(0.82, "npc"), y = unit(0.6, "npc"), just = "left", gp = gpar(fontsize = 11))
#   )
#   
#   # 3. Assembler les plots avec patchwork
#   combined_plots <- wrap_plots(plot_list, ncol = 2)
#   
#   # 4. Créer une page vide et placer chaque élément manuellement
#   grid.newpage()
#   pushViewport(viewport(layout = grid.layout(nrow = 3, heights = unit(c(1, 18, 2), "null"))))
#   
#   # Titre
#   print(main_title, vp = viewport(layout.pos.row = 1))
#   
#   # Plots
#   print(combined_plots, vp = viewport(layout.pos.row = 2))
#   
#   # Légende
#   grid.draw(editGrob(legend_grob, vp = viewport(layout.pos.row = 3)))
# # 
# # Make one plot show the legend
# plot_list[[1]] <- plot_list[[1]] + 
#   theme(legend.position = "bottom") +
#   labs(color = "Observation Type") +
#   scale_color_manual(
#     values = c("0" = "blue", "1" = "red"),
#     labels = c("Inliers", "True outliers")
#   )
# 
# # Création d'un plot spécial juste pour la légende
# legend_plot <- ggplot(data.frame(
#   x = c(1, 1),
#   y = c(1, 1),
#   label = factor(c(0, 1)),
#   linetype = factor("Cutoff")
# )) +
# #  geom_point(aes(x, y, color = label), size = 3) +
#  # geom_hline(aes(yintercept = 0.5, linetype = linetype), 
#   #           color = "black", size = 0.8) +
#   scale_color_manual(
#     values = c("0" = "blue", "1" = "red"),
#     labels = c("Inliers", "True outliers"),
#     name = "Points") +
#   scale_linetype_manual(
#     values = "dashed",
#     name = "Threshold",
#     labels = "Cutoff") +
#   theme_minimal() +
#   theme(
#     legend.position = "bottom",
#     legend.box = "horizontal",
#     legend.spacing.x = unit(0.5, 'cm'),
#     legend.title = element_text(size = 10),
#     legend.text = element_text(size = 9)
#   )
# 
# # Extraction de la légende
# shared_legend <- cowplot::get_legend(legend_plot)
# 
# # Combinaison finale
# full_title <- "k = 6.57    l = 19.02    rho1 = 0.845    d = 100    Streaming method"
# 
# final_plot <- wrap_plots(plot_list, ncol = 2) / 
#   legend_plot +
#   plot_layout(heights = c(10, 1)) +  # 90% pour les plots, 10% pour la légende
#   plot_annotation(
#     title = full_title,
#     theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
#   )
# 
# print(final_plot)
# # Combine with shared legend
# final_plot <- wrap_plots(plot_list, ncol = 2) + 
#   plot_layout(guides = "collect") &
#   theme(legend.position = "bottom")
# 
# final_plot <- final_plot + 
#   plot_annotation(title = full_title)

# 
# # 
# # for (r in contamin_rate){
# # print(paste0("r = ",r))
# # dataFile <- paste0('SimData-d', d, '-n', n, '-k', k, '-l', l, '-rho', rho1,'-r',r , '-sim', sim,".RData")
# # setwd(simDir)
# # contParam = ParmsF1(m1, k, l, rho1)
# # data = genererEchantillon(n,n,mu1 = mu0,mu2 = contParam$mu1,Sigma1 = Sigma0,Sigma2 = contParam$Sigma1,r)
# # save(data,file  = dataFile)
# # #load(dataFile)
# # print(paste0("nb_outliers = " ,sum(data$labelsVrais == 1)))
# # setwd(resDir)
# # fitFile = paste0('FitParms-d', d,  '-n', n, '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
# # load(fitFile)
# # 
# # if(r == contamin_rate[1] | r == contamin_rate[3]| r == contamin_rate[5]){
# # ylab = "Distances"
# # }
# # 
# 
# 
# 
# p = plot_outliers_distances(dist_sample = fitNaif$distances,dist_online = fitUsOnline$distances,dist_streaming = fitUSStreaming$distances,labels_vrais = data$labelsVrais,ylab = ylab,dim = 100,rate = r)
# ylab = ""
# plot_list[[length(plot_list) + 1]] <- p
# print(p)
# }
# Combine les plots
combined_plot <- wrap_plots(plot_list, ncol = 2)

legend_data <- data.frame(
  x = c(1, 2),
  y = c(1, 1),
  label = factor(c("Inliers", "True outliers"), levels = c("Inliers", "True outliers"))
)

full_title <- "k = 6.57    l = 19.02    rho1 = 0.845    d = 100    Streaming method"

# Texte de légende en dur (sous forme de plot vide contenant du texte)
manual_legend <- grid::textGrob(
  "Blue = Inliers    Red = True Outliers    Black Dashed Line = Cutoff Threshold",
  gp = grid::gpar(fontsize = 10)
)

# Combiner les plots + légende texte
final_plot <- (wrap_plots(plot_list, ncol = 2)) +
  plot_annotation(
    title = full_title,
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )
print(final_plot)


# pvraies = plot_outliers_distances(labels_vrais = data$labelsVrais,distances = vraiesDistances,ylab = "Distance",method = "vraies distances",dim = 100,rate = 40)
# 
# pnaif = plot_outliers_distances(labels_vrais = data$labelsVrais,distances = resNaif$distances,ylab = "Distance",method = "naive",dim = 100,rate = 40)
# 
# pStream = plot_outliers_distances(labels_vrais = data$labelsVrais,distances = resUsStreaming$distances,ylab = "Distance",method = "streaming",dim = 100,rate = 40)
# 
# 
# p1 = plot_outliers_distances(labels_vrais = data$labelsVrais,distances = resUsStreaming$distances,ylab = "Distance",method = method,rate = 10)
# 
# p2 = plot_outliers_distances(labels_vrais = data$labelsVrais,distances = resUsStreaming$distances,method = method,rate = 20)
# 
# p3 = plot_outliers_distances(labels_vrais = data$labelsVrais,distances = resUsStreaming$distances,ylab = "Distance",method = method,rate = 30)
# 
# p4 = plot_outliers_distances(labels_vrais = data$labelsVrais,distances = resUsStreaming$distances,method = method,rate = 30)
# 
# combined_plot <- (
#   (p1 | p2) /
#     (p3 | p4)
# ) + 
#   plot_layout(guides = "collect") & 
#   theme(legend.position = "bottom")


#################################################################################
#######Affichage erreurs estimation val propres Sigma0 et SigmaEst###############
#################################################################################
library(ggplot2)
library(patchwork)

# Paramètres
p <- nrow(Sigma0)                # dimension de la matrice
n_iter <- dim(resUsStreaming$Sigma)[1]  # nombre d'itérations

# 1. Calcul des vraies valeurs propres triées (lambda_0)
true_eigenvalues <- sort(eigen(Sigma0, symmetric = TRUE)$values, decreasing = TRUE)

# 2. Initialiser la matrice d'erreurs (p × n_iter)
df_errors <- matrix(NA, nrow = p, ncol = n_iter)

# 3. Calculer l'erreur signée pour chaque itération et chaque composante
for (iter in 1:n_iter) {
  # Eigenvalues estimées triées à l'itération iter
  est_eigvals <- sort(eigen(resUsStreaming$Sigma[iter, , ], symmetric = TRUE)$values, decreasing = TRUE)
  # Différence signée (hat(lambda) - lambda0)
  df_errors[, iter] <- est_eigvals - true_eigenvalues
}

# 4. Création des graphiques
plot_list <- lapply(1:p, function(i) {
  df_i <- data.frame(
    Iteration = 1:n_iter,
    Error = df_errors[i, ]
  )
  
  ggplot(df_i, aes(x = Iteration, y = Error)) +
    geom_line(color = "steelblue", size = 0.7) +  # Estimated eigenvalue error
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.6) +  # True eigenvalue reference
    ylim(c(-0.9,0))+
    labs(
      title = paste("Eigenvalue", i),
      x = ifelse(i > 5, "Iteration", ""),
      y = ifelse(i %% 5 == 1, expression(hat(lambda) - lambda[0]), "")
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.title.x = element_text(size = 9),
      axis.title.y = element_text(size = 9)
    )
})

# 5. Légende manuelle
legend_df <- data.frame(
  x = c(1, 2),
  y = c(1, 1),
  type = c("Estimated eigenvalue", "True eigenvalue")
)

legend_plot <- ggplot(legend_df, aes(x = x, y = y, linetype = type, color = type)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("Estimated eigenvalue" = "steelblue", "True eigenvalue" = "black")) +
  scale_linetype_manual(values = c("Estimated eigenvalue" = "solid", "True eigenvalue" = "dashed")) +
  guides(color = guide_legend(title = NULL), linetype = guide_legend(title = NULL)) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 11)
  )

# 6. Assemblage final
combined <- wrap_plots(plot_list, ncol = 5, nrow = 2)

final_plot <- combined / legend_plot +
  plot_layout(heights = c(10, 1)) +
  plot_annotation(
    title = "(k,l,rho1) = (0.86,0.56;0.6) 40 % of outliers streaming method",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  )

print(final_plot)


