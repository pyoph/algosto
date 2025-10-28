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
for(kl in 2:klNb){
  parms1 <- ParmsF1(m1=m1, k1=k1val[kl], l1=l1valup1[kl], rho1=rho1val[kl])
  X1 <- rmvnorm(B, mean=parms1$mu, sigma=parms1$Sigma)
  d1 <- rowSums((X1%*%Sigma0invHalf)^2)
  f1 <- density(log10(d1))
  lines(f1, col= klcol[compt+1]) 
  q1_95 <- quantile(d1, probs = .95)
  abline(v=log10(quantile(d1, probs=.95)), lty=2,col = klcol[(compt+1)])
  text(
    x = log10(quantile(d1, probs = .95)) - 0.05,
    y = 2.0 - 0.2 * compt,
    labels = paste0(round(q1_95, 2)),  # texte numérique seulement
    col = klcol[compt + 1],           # couleur appliquée ici
    pos = 4,
    cex = 0.8
  )
  text(0.98*f1$x[which.max(f1$y)], 1.02*max(f1$y), labels= kl-1, col=klcol[compt+1])
compt = compt + 1
  }
dev.off()

contamin_rate = c(0,5,20,40)
k = k1val[3]; l = l1valup1[3]; rho1 = rho1val[3]
# Distribution (log scale)
U <- matrix(rnorm(d*B), B, d)
X0 <- rmvnorm(B, mean=mu0, sigma=Sigma0)
d0 <- rowSums((X0%*%Sigma0invHalf)^2)
f0 <- density(log10(d0))

cutoffcorr = rep(0,n)
distancescorr = rep(0,n)
xMin <- 1e-1; xMax <- 1e11
distances = rep(0,n)
outliers =matrix(0,n,length(contamin_rate))
labelsvrais = matrix(0,n,length(contamin_rate))


for (m in seq_along(contamin_rate) ){
r = contamin_rate[m]
  contParam = ParmsF1(m1, k, l, rho1)
data = genererEchantillon(n,n,mu1 = mu0,mu2 = contParam$mu1,Sigma1 = Sigma0,Sigma2 = contParam$Sigma1,r)
resUsStreaming = StreamingOutlierDetection(data$Z,batch = ncol(data$Z))

labelsvrais[,m] = data$labelsVrais

file = paste0("distancedensity-r",r,".pdf")
setwd("~")
pdf(file,width = 8,height = 6)

plot(f0, xlim=c(log10(xMin), log10(xMax)), ylim=c(0, 2.2),col = "black", main='', lwd = "4")  

d1 <- resUsStreaming$distances

f1 = density(log10(d1))
lines(f1$x,f1$y, col= adjustcolor("darkgreen", alpha.f = 0.6),lwd = "4") 

for(j in (1:n)){
  
  cutoffcorr[j] = qchisq(p = .5,df = d)/median(d1[1:j])
  
  if (d1[j] > qchisq(p = .95,df = d)){outliers[j,m] = 1}
  
  distancescorr[j] = cutoffcorr[j]*d1[j]
}


d1corr = distancescorr

f1corr = density(log10(d1corr))
lines(f1corr$x,f1corr$y, col= adjustcolor("darkred", alpha.f = 0.6),lwd = "4")
dev.off()

}#

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

