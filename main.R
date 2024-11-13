#install.packages("microbenchmark")
#install.packages("matlib")
#install.packages("MASS")
#install.packages("reshape2")
#install.packages("corrplot")
#install.packages("dlpyr")
#install.packages("devtools")
#install.packages("usethis")

library(reshape2)
library(RobRegression)
library(Gmedian)
library(ggplot2)
library(far)
library(gridExtra)
library(microbenchmark)
library("matlib")
library("MASS")
library("corrplot")
library("dplyr")
source("~/work/algosto/parametres.R")
source("~/work/algosto/simulations.R")
source("~/work/algosto/algorithmes.R")
source("~/work/algosto/resultats.R")
#usethis::create_package("C:/Users/Paul/Documents/codeRTheseLPSM")
#devtools::load_all("C:/Users/Paul/Documents/codeRTheseLPSM")


n <- 1e4

d <- 10

mu1 <- rep(0,d)

mu2 <- 5*rep(1,d)

Sigma1 <- diag(sqrt(1:10))


Sigma2 <- 10*diag(1:10)

set.seed(123)

#Création d'une matrice de covariance non nécessairement diagonale
A <- matrix(runif(100, min=-1, max=1), nrow=10)
A <- (A + t(A)) / 2  
Sigma1 <- A %*% t(A)
Sigma1=(Sigma1+diag(10))/10
eigen(Sigma1)$values

#Pourcentage de non outliers
p1 <- 0.98

#Pourcentage d'outliers
p2 <- 1-p1

resultsSimul <- genererEchantillon(n,d,mu1,mu2,p1,p2,Sigma1 =Sigma1 ,Sigma2)


Z <- resultsSimul$Z

params <- initialiser_parametres(
  Y = Z,
  c = sqrt(10),
  Sigma = Sigma1,
  #Sigma = Sigma1,
  r = 1.5,
  k = 1
)
#cov(Z)
params$Sigma

#eigen(resultsSimul$Vvrai)$vectors

#resultsSimul$labelsVrais
#Z
results <- estimMVOutliers(Z,params$c ,params$n,params$d,params$d,params$r,aa = 0.75,niter = 1e4)

#Détection des outliers après 10^3 itérations
results2 <- estimMVOutliers(Z,params$c ,n = 1e4,params$d,params$d,params$r,aa = 0.75,niter = 1e4,depart = 1e2 + 1,minit = results$m, Vinit = results$V,
                            U = results$U,stat = results$stat,vpMCM = results$vpMCM, lambda = results$lambda,lambdatilde = results$lambdatilde)

resultsStr <- streaming(Z,params$n/params$k,params$k,params$c,params$n,params$d,params$q,params$r)

#lambdaResultat <- RobbinsMC2(1e4,results$vpMCM[99,],10^(-8),alpha = 0.75,c= 2,w=2,samp=1e4, initbarre = results$vpMCM[99,],cbarre=0)

#lambdaResultat$lambdalist

#eigen(params$Sigma1)$values

#Affichage de l'estimation des erreurs pour l'estimation de mbarre

affiche_erreurs_mediane_geometrique(calculErreursM(results$miter,rep(0,d),params$n/params$k),params$n/params$k)

#Affichage de l'estimation des erreurs pour l'estimation de Vbarre

affiche_erreurs_mcm(calculErreursV(results$VIter,resultsSimul$Vvrai,params$n/params$k),params$n/params$k,params$c)

#Affichage de l'erreur d'estimation des valeurs propres de la MCM

affiche_erreurs_vp_mcm(calculErreursVpMCM(results$vpMCM,resultsSimul$VpvraiesV,params$n/params$k),params$n/params$k)

#results$U[99,,]

#Affichage de l'erreur d'estimation des vecteurs propres de la MCM

affiche_erreurs_vectp_mcm(calculErreursVectpMCM(results$U,eigen(params$Sigma)$vectors,params$n/params$k,params$d),params$n/params$k,params$c)
#Phi = results$phijn
#Phi[9999,]

#affiche_erreurs_vectp_mcm(calculErreursVectpMCM(results$Uphi,eigen(params$Sigma)$vectors,params$n/params$k,params$d),params$n/params$k,params$c)
#results$Uphi[9000:9999,,]
#resultsSimul$VcovVrai
#Affichage de l'erreur d'estimation des valeurs propres de la matrice de covariance
affiche_erreurs_cov(calculErreursValPropreCov(results$lambdaIter,eigen(params$Sigma)$values,params$n/params$k),params$n/params$k)
#Tracer avec les vraies valeurs 

#Comparaison stat khi-deux avec densité théorique Khi deux

densiteHistKhi2(results$stat,params$d)
#results$outlier_labels
table_contingence(resultsSimul$labelsVrais[1:9999],results$outlier_labels[1:9999])
plot(1:(n-1),cumsum(results$outlier_labels)/(1:(n-1)),ylim = c(0,1),'l');abline(h = ((p1)*0.05 +1- p1))
#length(results2$outlier_labels)

#affiche_erreurs_Sigma(calculErreursSigma(params$Sigma,results$U,results$lambdaIter,params$n/params$k,params$d),params$n/params$k,params$c)
VP <- results$U[999999,,] %*% diag(1/sqrt(colSums(results$U[999999,,]^2)))
VP <- results$U[999999,,] %*% diag(1/sqrt(colSums(results$U[999999,,]^2)))
SigmaEstim <- VP %*%diag(results$lambdaIter[999999,])%*%t(VP)

plot(params$Sigma,SigmaEstim) 
abline(0,1)
erreursSigma <- calculErreursSigma(params$Sigma,results$U,results$lambdaIter,params$n/params$k,params$d)
erreursSigmaMoyenne <- erreursSigma$erreursSigmaMoyenne
affiche_erreurs_Sigma(erreursSigmaMoyenne,params$n/params$k,params$c)

erreursSigmaMoyenne[9999]
erreursSigmaMoyenne[999999]
erreursSigmaMoyenne[9999]
erreursSigmaMoyenneRep <- matrix(0,n,10)
resultsRep <- rep(0,10)
for (i in 1:10){
  resultsRep[i] <- estimMVOutliers(Z,params$c,params$n,params$d,params$d,params$r)
  erreursSigmaMoyenneRep[,i] <- calculErreursSigma(params$Sigma,resultsRep[i]$U,resultsRep[i]$lambdaIter,params$n/params$k,params$d)
  affiche_erreurs_Sigma(erreursSigmaMoyenneRep[,i],params$n/params$k,params$c)
  
}
#densiteHistKhi2(statTheorique,params$d)
outliervrais <- numeric(n)  
# Boucle pour calculer les valeurs des outliers
for (i in 1:n) {
  # Vérifiez que Z, resultsSimul, et les éléments passés à Outlier sont correctement définis
  outliervrais[i] <- Outlier(
    Z[i, ],                     # L'élément de Z à traiter
    0.05,                       # Niveau de confiance
    resultsSimul$VcovVrai$vectors,  # Vecteurs propres de la matrice
    rep(0, d),                  
    eigen(params$Sigma)$values   # Valeurs propres de la matrice
  )
}
# Conversion en vecteurs pour éviter les erreurs de type 'list'
resultsSimul$labelsVrais <- unlist(resultsSimul$labelsVrais)  # Conversion en vecteur
outliervrais <- unlist(outliervrais)                          # Conversion en vecteur

# Vérification et appel de table_contingence
table_contingence(resultsSimul$labelsVrais, outliervrais)

table_contingence(as.vector(resultsSimul$labelsVrais),as.vector(outliervrais))
length(outliervrais)

invSigma <- solve(params$Sigma)

# Boucle pour calculer les valeurs des outliers
for (i in 1:n) {
  # Vérifiez que Z, resultsSimul, et les éléments passés à Outlier sont correctement définis
  outliervrais[i] <- Outlier(
    Z[i, ],                     # L'élément de Z à traiter
    0.05,                       # Niveau de confiance
    resultsSimul$VcovVrai$vectors,  # Vecteurs propres de la matrice
    rep(0, d),                  
    eigen(params$Sigma)$values   # Valeurs propres de la matrice
  )
}
# Conversion en vecteurs pour éviter les erreurs de type 'list'
resultsSimul$labelsVrais <- unlist(resultsSimul$labelsVrais)  # Conversion en vecteur
outliervrais <- unlist(outliervrais)                          # Conversion en vecteur

#params$Sigma
# Vérification et appel de table_contingence
table_contingence(resultsSimul$labelsVrais, outliervrais)

table_contingence(as.vector(resultsSimul$labelsVrais),as.vector(outliervrais))
length(outliervrais)
dim(Z)
outliervrais <- numeric(n)  
# Boucle pour calculer les valeurs des outliers

var(Z[1,])
var(Z)


invSigma <- solve(params$Sigma)



params$Sigma


S <- sum(Z1^2)

pchisq(S, df = params$d,lower.tail =FALSE)
#print(pha

VPSigma <- results$U[i,,] %*% diag(1/sqrt(colSums(VPSigma^2)))
#VPSigma <- VPSigma %*% diag(1/sqrt(colSums(VPSigma^2)))
lambda <- results$lambdaIter[9999,]

for (i in 1:n) {

  #VPSigma <- eigen(params$Sigma)$vectors
  #VPSigma <- orthonormalization(VPSigma,basis = TRUE, norm = TRUE)
  
  #lambda <- eigen(params$Sigma)$values
  
  
  S <- 0
  

      for (j in 1:length(lambda)) {
     S <- S + 1/lambda[j] * (t(VPSigma[,j]) %*% (Z[i,]))^2 
    
      }
  
  
  #S <- sum(S^2)
  
  #Zi <- 1/sqrt(var(Z[i,])) %*% Z[i,]
  
  
  phat <- pchisq(S, df = params$d,lower.tail =FALSE)
  #print(phat)
  # Détection de l'outlier
  #print(phat)
  
  
  
  if (phat < 0.05) {
    outliervrais[i] <- 1  # Indiquer qu'il s'agit d'un outlier
  } else {
    outliervrais[i] <- 0  # Indiquer qu'il ne s'agit pas d'un outlier
  }
  
  
}
# Conversion en vecteurs pour éviter les erreurs de type 'list'
resultsSimul$labelsVrais <- unlist(resultsSimul$labelsVrais)  # Conversion en vecteur
outliervrais <- unlist(outliervrais)                          # Conversion en vecteur

# Vérification et appel de table_contingence
table_contingence(resultsSimul$labelsVrais, outliervrais)


#Calcul de la distance de Mahalanobis
mahaDist <- mahalanobis(x = Z, center = rep(0,d), cov = params$Sigma)

# Cutoff values from chi-square distribution to identify outliers
cutoff <- qchisq(p = 0.95, df = d)

# Créer un data frame combinant Z et mahaDist
df <- data.frame(Z, Mahalanobis_Distance = mahaDist, labels_vrais = resultsSimul$labelsVrais)

# Afficher le data frame
print(df)

#Calcul des p-values
df$p <- pchisq(df$Mahalanobis_Distance,df = d,lower.tail =   FALSE )

#Identification des outliers
df <- df %>%
  mutate(Outlier = ifelse(Mahalanobis_Distance > cutoff, 'Yes','No'))
# Calculer la proportion d'outliers
proportion_outliers <- df %>%
  summarise(Proportion = mean(Outlier == 'Yes')) %>%
  pull(Proportion)

# Afficher la proportion d'outliers
print(proportion_outliers)

# Création de la table de contingence
table_contingence <- table(df$Outlier, df$labels)

# Afficher la table de contingence
print(table_contingence)
