# KL divergence btw F0 and F1

rm(list=ls())
library(nlme); library(flexmix); library(fields)

# Functions
KL <- function(parms1, parms2){
  invSigma2 <- solve(parms2$Sigma)
  0.5*(log(det(parms2$Sigma)/det(parms1$Sigma)) - d + sum(diag(invSigma2%*%parms1$Sigma)) +
         t(parms2$mu-parms1$mu)%*%invSigma2%*%(parms2$mu-parms1$mu))[1, 1]
}

# Contamination parms: F1
ParmsF1 <- function(m1, k1, l1, rho1){
  d <- length(m1)
  mu1 <- k1*m1
  sigmaSq1 <- l1*sigmaSq0
  Sigma1 <- diag(sqrt(sigmaSq1)) %*% toeplitz(rho1^(0:(d-1))) %*% diag(sqrt(sigmaSq1))
  return(list(mu1=mu1, Sigma1=Sigma1))
}


# Dims
d <- 10

# Null parms: F0
mu0 <- rep(0, d)
sigmaSq0 <- (1:d); sigmaSq0 <- sigmaSq0 / mean(sigmaSq0)
rho0 <- 0.3
Sigma0 <- diag(sqrt(sigmaSq0)) %*% toeplitz(rho0^(0:(d-1))) %*% diag(sqrt(sigmaSq0))
parms0 <- list(mu=mu0, Sigma=Sigma0)

# # Null (null) parms: F0
# mu0 <- rep(0, d)
# sigmaSq0 <- rep(1, d)
# rho0 <- 0
# Sigma0 <- diag(sqrt(sigmaSq0)) %*% toeplitz(rho0^(0:(d-1))) %*% diag(sqrt(sigmaSq0))
# parms0 <- list(mu=mu0, Sigma=Sigma0)

# Contamination
m1 <- rep(1/sqrt(d), d)
KL(parms1=parms0, parms2=ParmsF1(m1, 1, 1, 0.5))

par(mfrow=c(3, 2))
k1grid <- 2^(-4:5); 
plot(k1grid, sapply(k1grid, function(k1){KL(parms0, ParmsF1(m1, k1, 1, rho0))}), log='xy')
l1grid <- 2^(-5:5); 
plot(l1grid, sapply(l1grid, function(l1){KL(parms0, ParmsF1(m1, 0, l1, rho0))}), log='xy')
rho1grid <- seq(0.005, 0.995, by=0.005); 
plot(rho1grid, sapply(rho1grid, function(rho1){KL(parms0, ParmsF1(m1, 0, 1, rho1))}), log='y')

KLk1l1 <- sapply(k1grid, function(k1){sapply(l1grid, function(l1){KL(parms0, ParmsF1(m1, k1, l1, rho0))})})
image.plot(log10(l1grid), log10(k1grid), log10(KLk1l1))
KLk1rho1 <- sapply(k1grid, function(k1){sapply(rho1grid, function(rho1){KL(parms0, ParmsF1(m1, k1, 1, rho1))})})
image.plot(rho1grid, log10(k1grid), log10(KLk1rho1))
KLrho1l1 <- sapply(rho1grid, function(rho1){sapply(l1grid, function(l1){KL(parms0, ParmsF1(m1, 0, l1, rho1))})})
KLrho1l1[KLrho1l1==0] <- NA
image.plot(log10(l1grid), rho1grid, log10(KLrho1l1))


##########################
#Tests effets sur norme de Frobenius variation k
#############################

contamin_rate = seq(0,40,5)

methodes = c("streaming","online")

erreurNormeFrobenius = array(0, dim = c(length(contamin_rate), length(methodes),length(k1grid)))
faux_positifs = array(0, dim = c(length(contamin_rate), length(methodes),length(k1grid)))
faux_negatifs = array(0, dim = c(length(contamin_rate), length(methodes),length(k1grid)))

for (i in seq_along(k1grid))
{
for (j in seq_along(contamin_rate))
{
  
  r = contamin_rate[j]
  print(paste("r ", r))
  print(paste("k ",k))
  k = k1grid[i]

    data <- genererEchantillon(n,d,mu1,mu2 = k*mu1,p1 = 1- r/100,r/100,Sigma1,Sigma2 = Sigma1,contamin = "moyenne",cluster = FALSE)
  Z = data$Z
  compt = 1
for (m in methodes)  
{
  
  if(m == "online"){
    print (m)
    resultats = StreamingOutlierDetection(Z,batch = 1)
  Sigma = resultats$Sigma[nrow(Z),,]
  outliers = resultats$outlier_labels
  }
  if(m == "streaming"){
    print (m)
    resultats = StreamingOutlierDetection(Z,batch = ncol(Z))
    Sigma = resultats$Sigma[nrow(Z),,]
    outliers = resultats$outlier_labels
  }
  
  
  print(paste("i = ",i))
  print(paste("compt = ",compt))
  print(paste("j = ",j))
  #print("m = ",m)
  erreurNormeFrobenius[j,compt,i] = norm(Sigma - Sigma1,"F")  
  
  
  
  tc <- table(data$labelsVrais[1:(nrow(Z))], as.numeric(outliers)[1:(nrow(Z))])
  if ("0" %in% rownames(tc)){
  if((tc["0","0"] + tc["0","1"]) != 0)
  {faux_positifs[j,compt,i]   <-(tc["0", "1"]/(tc["0", "1"] + tc["0", "0"]))*100}
  }
  if ("1" %in% rownames(tc)){
  if((tc["1","0"] + tc["1","1"]) != 0)
  {faux_negatifs[j,compt,i]   <-(tc["1", "1"]/(tc["1", "1"] + tc["1", "0"]))*100}}
  compt = compt + 1
} 
  
}

}


plot(k1grid,faux_positifs[,1,])

fp_seuil =array(0,dim = c(length(contamin_rate), length(methodes)))

for (k in seq_along(contamin_rate))
{
 r = 0 
 compt = 1
  data <- genererEchantillon(n,d,mu1,mu2 = k*mu1,p1 = 1- r/100,r/100,Sigma1,Sigma2 = Sigma1,contamin = "moyenne",cluster = FALSE)
  Z = data$Z
  for (m in methodes){
  if(m == "online"){
    print (m)
    resultats = StreamingOutlierDetection(Z,batch = 1)
    Sigma = resultats$Sigma[nrow(Z),,]
    outliers = resultats$outlier_labels
    cutoffCorr = qchisq(.95,df = d)*median(resultats$distances)/qchisq(.5,df = d)
    resultats = StreamingOutlierDetection(Z,batch = 1,cutoff = cutoffCorr)
    outliers = resultats$outlier_labels
    
  }
  if(m == "streaming"){
    print (m)
    resultats = StreamingOutlierDetection(Z,batch = ncol(Z))
    Sigma = resultats$Sigma[nrow(Z),,]
    outliers = resultats$outlier_labels
    cutoffCorr = qchisq(.95,df = d)*median(resultats$distances)/qchisq(.5,df = d)
    resultats = StreamingOutlierDetection(Z,batch = 1,cutoff = cutoffCorr)
    outliers = resultats$outlier_labels
    
    
  }
  
  
  tc <- table(data$labelsVrais[1:(nrow(Z))], as.numeric(outliers)[1:(nrow(Z))])
  if ("0" %in% rownames(tc)){
    if((tc["0","0"] + tc["0","1"]) != 0)
    {fp_seuil[k,compt]   <-(tc["0", "1"]/(tc["0", "1"] + tc["0", "0"]))*100
    print(fp_seuil[k,compt])}
  }
  if ("1" %in% rownames(tc)){
    if((tc["1","0"] + tc["1","1"]) != 0)
    {fp_seuil[k,compt]   <-(tc["1", "1"]/(tc["1", "1"] + tc["1", "0"]))*100
    print(fp_seuil[k,compt])}}
  compt = compt + 1}
}  
  



# Paramètres
alpha <- 0.05
conf_level <- 1 - alpha
p_seq <- seq(0.6, 1, length.out = 200)
N <- nrow(Z)  # N doit être défini au préalable

# Calcul des n0 et des x théoriques
n0 <- floor(p_seq * N)
x <- round(n0 * alpha)

# Calcul de l’intervalle de confiance exact (Clopper-Pearson)
conf <- binom.confint(x = x, n = n0, conf.level = conf_level, methods = "exact")


# Tracé des bornes
plot(p_seq, conf$lower, type = "l", col = "blue", lwd = 2,
     ylim = c(0, 0.15), xlab = "p (proportion pour n0 = p * N)",
     ylab = "Intervalle de confiance de la proportion",
     main = "IC à 95% pour Bin(n0, alpha = 0.05)")
lines(p_seq, conf$upper, col = "red", lwd = 2)

fp_positions = c(1.00, 0.95, 0.9,  0.85, 0.80, 0.75, 0.70, 0.65,0.60)
length(fp_positions)
fp_streaming_values = fp_seuil[,1]/100
length(fp_seuil[,2])
points(fp_positions, fp_streaming_values, col = "orange", pch = 19)

