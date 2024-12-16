#install.packages("microbenchmark")
#install.packages("matlib")
#install.packages("MASS")
#install.packages("reshape2")
#install.packages("corrplot")
#install.packages("dlpyr")
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
#Fonction de détection des outliers prend en paramètre une nouvelle donnée
Outlier <- function(donnee, seuil_p_value, VP, m, lambda) {
  
  
  vectPV <- VP %*% diag(1/sqrt(colSums(VP^2)))
  
  cutoff <- qchisq(p = 0.975, df = 10)
  
  S <- 0
  
  # Calcul de la statistique S
    
    for (j in 1:length(lambda)) {
      #S <- S + 1/lambda[j] * ((donnee - m)%*% vectPV[,j])^2
      S <- S + 1/lambda[j] * (sum((donnee - m)* vectPV[,j]))^2
    }
  
  
  # Calcul de la p-value basée sur la statistique du Chi2
  phat <- pchisq(S, df = length(lambda),lower.tail =FALSE)
  #print(phat)
  # Détection de l'outlier
  #if (phat < seuil_p_value) {
   # outlier_label <- 1  # Indiquer qu'il s'agit d'un outlier
   #} else {
    #outlier_label <- 0  # Indiquer qu'il ne s'agit pas d'un outlier
  #}
  
  if (S > cutoff) {outlier_label <- 1}
  else {outlier_label <- 0}
  
  # Retourner le label, la p-value et la statistique S
  return(list(outlier_label = outlier_label, phat = phat, S = S))
}
#Estimation des valeurs propres de Sigma

RobbinsMC=function(mc_sample_size=10000,vp,epsilon=10^(-8),alpha=0.75,c=2,w=2,samp=mc_sample_size,init=detoile)
{
  p=length(vp)
  vp2=init
  lambda=init
  lambdalist=c()
  vplist=c()
  Y=matrix(rnorm(mc_sample_size*p),ncol=p)
  # X2=matrix(rnorm(mc_sample_size*p),ncol=p)
  slog=1
  for (i in 1:mc_sample_size)
  {
    Z=Y[i,]
    # Z2=X2[i,]
    E1=    Z^2*(sum(( (vp)-lambda*(Z^2))^2) + sum((lambda * Z^2)%*%t((lambda * Z^2))) - sum(((lambda * Z^2)^2))  )^(-0.5)
    vp0=vp2
    E2=    (sum(( (vp)-lambda*(Z^2))^2) + sum((lambda * Z^2)%*%t((lambda * Z^2))) - sum(((lambda * Z^2)^2))  )^(-0.5)
    lambda =lambda  - c*i^(-alpha)*lambda*E1 + c*i^(-alpha)* (vp)*E2
    slog=slog+log(i+1)^w
    vp2=vp2+log(i+1)^w *((slog)^(-1)) *(lambda - vp2)
    eps=sqrt(sum((vp2-vp0)^2))
    if (length(which(samp==i) > 0))
    {
      lambdalist=rbind(lambdalist,lambda)
      vplist=rbind(vplist,vp2)
    }
    #   if ( eps < epsilon ) break;
  }
  return(list(vp=vp2,niter=i, lambdalist=lambdalist, vplist=vplist))
}



#Détection offline des outliers

detectionOffline <- function(Z,SigmaEstim,m,seuil_p_value, cutoff =qchisq(p = 0.975, df = ncol(Z)))

{
  
  outliers_labels <- rep(0,n)
  
   #cutoff <- qchisq(p = 0.975, df = ncol(Z))
  
    
  vectPV <- eigen(SigmaEstim)$vectors
  lambda <- eigen(SigmaEstim)$values
  vectPV <- vectPV %*% diag(1/sqrt(colSums(vectPV^2)))
  #cutoff <- qchisq(p = 0.975, df = ncol(Z))

  
  
  #m = rep(0,d)
  for (i in 1:nrow(Z)) 
  {
    
    
    S <- 0
    
    # Calcul de la statistique S
    
    for (j in 1:length(lambda)) {
     S <- S + 1/lambda[j] * (sum((Z[i,] - m)* vectPV[,j]))^2 
    
    }
    
    
    #S <- t(Z[i,] - m)%*%solve(SigmaEstim)%*%(Z[i,] - m) 
    
    # Calcul de la p-value basée sur la statistique du Chi2
    phat <- pchisq(S, df = ncol(Z),lower.tail =FALSE)
    
    # Calcul de la p-value basée sur la statistique du Chi2
    #phat <- pchisq(S, df = ncol(Z),lower.tail =FALSE)
    #print(phat)
    # Détection de l'outlier
    #if (phat < seuil_p_value) {
     # outliers_labels[i] <- 1  # Indiquer qu'il s'agit d'un outlier
      #print("OK")
     #} else {
      #outliers_labels[i] <- 0  # Indiquer qu'il ne s'agit pas d'un outlier
    #}
    
    if (S > cutoff) {outliers_labels[i] <- 1}
    else {outliers_labels[i] <- 0}
    #}
  }
    
  #print(which(outliers_labels == 1))
  return (outliers_labels)
  
}


#Détection des outliers à partir d'un vecteur de distances et d'un cutoff


detectionOutliers <- function(distances, cutoff =qchisq(p = 0.975))
  
{
  
  outliers_labels <- rep(0,n)
  
  #cutoff <- qchisq(p = 0.975, df = ncol(Z))
  
  
  
  
  
  #m = rep(0,d)
  for (i in 1:nrow(Z)) 
  {
    
    
    
    if (distances[i] > cutoff) {outliers_labels[i] <- 1}
    else {outliers_labels[i] <- 0}
    #}
  }
  
  #print(which(outliers_labels == 1))
  return (outliers_labels)
  
}


#Estimation des valeurs propres de Sigma

RobbinsMC2=function(mc_sample_size=10000,vp,epsilon=10^(-8),alpha=0.75,c=2,w=2,samp=mc_sample_size,init=vp,initbarre = vp,cbarre = 0,ctilde = 0)
{
  p=length(vp)
  vp2=initbarre
  lambda=init
  lambdalist=c()
  vplist=c()
  Y=matrix(rnorm(mc_sample_size*p),ncol=p)
  # X2=matrix(rnorm(mc_sample_size*p),ncol=p)
  slog=1
  for (i in 1:mc_sample_size)
  {
    Z=Y[i,]
    # Z2=X2[i,]
    E1=    Z^2*(sum(( (vp)-lambda*(Z^2))^2) )^(-0.5)
    vp0=vp2
    E2=    (sum(( (vp)-lambda*(Z^2))^2)  )^(-0.5)
    lambda =lambda  - c*(i+ctilde)^(-alpha)*lambda*E1 + c*(i+ctilde)^(-alpha)* (vp)*E2
    slog=slog+log(i+1+ cbarre)^w
    vp2=vp2+log(i+1+ cbarre)^w *((slog)^(-1)) *(lambda - vp2)
    eps=sqrt(sum((vp2-vp0)^2))
    if (length(which(samp==i) > 0))
    {
      lambdalist=rbind(lambdalist,lambda)
      vplist=rbind(vplist,vp2)
    }
    #   if ( eps < epsilon ) break;
  }
  return(list(vp=vp2,niter=i, lambdalist=lambdalist, vplist=vplist,lambda = lambda))
}



#Fonction qui estime m, et V et détecte la présence d'outliers online
estimMVOutliers <- function(Y,c,n,d,q,r,aa,depart = 0,niter = n,minit = r*rnorm(d), Vinit = diag(d),
                            U = array(1, dim = c(n, q, q)),stat = rep(0,n-1),vpMCM = matrix(0,n,d), lambda = rep(1,d),lambdatilde = rep(1,d),methode = "eigen",seuil_p_value = 0.05,cutoff =qchisq(p = 0.975, df = d),niterRMon = d^2 )
{

  sampsize = d
  
  #Initialisation du vecteur avec les vrais paramètres

  #Sigma <- diag(sqrt(1:d))

  #Stockage des estimations des valeurs propres de la matrice de covariance
  lambdaIter = matrix(0,n,d)

  #Yv <- mvtnorm::rmvnorm(1000000,sigma=Sigma)

  ##Estimation de mn

  #Initialisation de m

  m <- minit


  #Sigma <- diag(d)

  #Initialisation de V

  V <- Vinit

  #Calcul de la vraie matrice de covariance médiane

  #Vvrai <- WeiszfeldCov(Yv,nitermax = 1000)$covmedian

  #VcovVrai <- GmedianCov(Yv, scores = q)


  #Calcul de la vraie médiane géométrique

  mvrai <- rep(0,d)

  #Stockage des outliers

  stat <- stat

  #Calcul de la vraie MCM par l'algorithme de Weiszfeld

  moyennem = m

  #VvectproprCovV <- VcovVrai$vectors

  moyenneV <- V

  #Stockage des itérations de matrices dans un tableau
  VIter <- array(0, dim = c(d, d, n))

  #Valeurs propres de V
  #vpMCM <- matrix(0,n,d)

  #Initialisation de lambda

  #lambda <- eigen(Vvrai)$values
  #lambda <- rep(1,d)
  #lambdatilde <- rep(1,d)

  #vp2 <- eigen(Vvrai)$values
  slog <- 1
  w <- 2
  #Compteur pour les non outliers

  #j <- 0


  #Stockage des vecteurs propres dans un tableau

#Si départ = 0 initialisation de U sur la sphère unité
if(depart == 0){
  
  #Vecteurs propres pour l'algorithme SGA
  Uphi <- array(1, dim = c(n, q, q))
  
  #Stockage des phijn
  phijn <- matrix(0,n,d)

  #Initialisation de U avec des vecteurs propres sur la sphère unité

  #U[1,,] <-  1.5*diag(q)
  
  matrix_random <- matrix(0, d, d)
  
  # Générer chaque colonne aléatoire sur la sphère unité
  for (i in 1:d) {
    v <- rnorm(d)  # Tirer un vecteur d composantes normales
    matrix_random[, i] <- v / sqrt(sum(v^2))  # Normaliser
  }
  
  U[1,,] <- matrix_random

  Uphi[1,,] <- matrix_random
}
  
  #Stockage des itérations
  miter = matrix(0,n,d)

  # Vecteur pour stocker les labels des outliers
  outlier_labels <- rep(0, niter-1)

  phat <- rep(0,n-1)

  for (i  in 1:(niter-1))
  {

    gamma = c/(i+depart)^(0.75)
    #gamma = c/i



    #Test pour gérer la non négativité de Vn

    if (gamma >= norm((Y[i+1,] - moyennem) %*% t(Y[i+1,] - moyennem) - V,type = "F"))
    {
    gamma = norm((Y[i+1,] - moyennem) %*% t(Y[i+1,] - moyennem) - V,type = "F")
    }




    #Mise à jour de mn

    m = m + gamma*((Y[i+1,] - m))/sqrt(sum((Y[i+1,]  - m)^2))





    #Mise à jour de V


    V = V + gamma*((Y[i+1,] - moyennem) %*% t(Y[i+1,] - moyennem) - V)/norm((Y[i+1,] - moyennem) %*% t(Y[i+1,] - moyennem) - V,type = "F")

    moyennem = moyennem*i/(i+1) + 1/(i+1)*m

    miter[i,] = moyennem

    #Mise à jour de moyenneV
  
    moyenneV <- (i+depart)/(i + depart +1)*moyenneV + 1/(i + 1)*V


    VIter[,,i] = moyenneV
   #U[i,,] <- orthonormalization(U[i,,],basis = TRUE, norm = TRUE)
    #Calcul des phijn
    #for (l in (1:q))
    #{
     # phijn[i,l] <- sum(t(Y[i+1,] - m)*Uphi[i,,l])
      #print(phijn[i,l])
    #}
    
    for (l in (1:q))
    {
      Un =  U[i,,l]/sqrt(sum(U[i,,l]^2))
      U[i+1,,l] <-  (1 - 1/(i + 1)^(aa))*U[i,,l] + 1/(i + 1)^(aa)*(moyenneV %*% Un)
      #U[i+1,,l] <- U[i,,l] + gamma*((Y[i+1,] - m)%*%t((Y[i+1,] - m)) %*% U[i,,l])
      
      #if (l  >= 2){
        
       # S <- rep(0,d)
      #for (iter in (1:(l - 1)))
      #{
       #S <- S + phijn[i,iter]*Uphi[i,,iter] 
      #}
      #}
      #Uphi[i+1,,l]= Uphi[i,,l] + gamma * phijn[i,l]*((Y[i+1,] - m) - phijn[i,l]*Uphi[i,,l] )
      
    }



    #matrice <- U[i+1,,]

    #Tri des vecteurs en ordre décroissant de norme
    #normes_colonnes <- apply(matrice, 2, function(col) sqrt(sum(col^2)))
    
    # order() pour obtenir les indices de tri par ordre décroissant des normes
    #indices_tri <- order(normes_colonnes, decreasing = TRUE)

    # Trier la matrice selon les normes décroissantes des colonnes
    #matrice <-matrice[,indices_tri]

    #U[i+1,,] <- matrice
    #Orthogonalisation des vecteurs propres
    U[i+1,,] <- orthonormalization(U[i+1,,],basis = TRUE, norm = FALSE)

    #vp = apply(U[i+1,,],2,norm)

    #Estimation des valeurs propres de la MCM
    #vpMCM[i,] = apply(U[i,,], 2, function(col) sqrt(sum(col^2)))
    
    #Rédiger partie simulation
    #vp  = lambdaEstim (précédente valeur de lambda) sample 100

    #Récupération des résultats de la fonction Outlier
    #resultatOutlier <- Outlier(donnee = Y[i, ],seuil_p_value = 0.05,VP = U[i,,],m = moyennem,lambdatilde)
    if(methode == "eigen"){
    VPropresV <- eigen(V)$vectors
    valPV <- eigen(V)$values
    lambdaResultat <- RobbinsMC2(1e2,vp=valPV,samp=1:sampsize,init = lambdatilde,initbarre = lambda,ctilde = sampsize*(i-1),cbarre = sampsize*(i-1))
    lambda <- lambdaResultat$vp
    lambdatilde <- lambdaResultat$lambda
    #print(lambda)
    lambdaIter[i,] <- lambdatilde
    
    resultatOutlier <- Outlier(donnee = Y[i, ],seuil_p_value = 0.05,VP = eigen(V)$vectors,m = moyennem,lambdatilde)
    
    }
    else {
      vpMCM[i,] = sqrt(colSums(U[i,,]^2))
      #Prendre norme de U[]
      #Z=Y[i,]
      # Z2=X2[i,]
      #E1=    Z^2*(sum(( (vp)-lambda*(Z^2))^2) + sum((lambda * Z^2)%*%t((lambda * Z^2))) - sum(((lambda * Z^2)^2))  )^(-0.5)
      #slog=slog+log(i+1)^w
      #vp2=vp2+log(i+1)^w *((slog)^(-1)) *(lambda - vp2)
      #vp0 = vp2
      #E2=    (sum(( (vp)-lambda*(Z^2))^2) + sum((lambda * Z^2)%*%t((lambda * Z^2))) - sum(((lambda * Z^2)^2))  )^(-0.5)
      #lambda = lambda  - c*i^(-0.75)*lambda*E1 + c*i^(-0.75)* (vp)*E2
      #Convergence si initialisation à vpMCM[i,] (= delta) au lieu de lambda
      lambdaResultat <- RobbinsMC2(1e2,vp=vpMCM[i,],samp=1:sampsize,init = lambdatilde,initbarre = lambda,ctilde = sampsize*(i-1),cbarre = sampsize*(i-1))
      lambda <- lambdaResultat$vp
      lambdatilde <- lambdaResultat$lambda
      #print(lambda)
      lambdaIter[i,] <- lambdatilde
      
      
      resultatOutlier <- Outlier(donnee = Y[i, ],seuil_p_value = 0.05,VP = U[i,,],m = moyennem,lambdatilde)}
    # Vérifier si l'entrée est un outlier
    if (resultatOutlier$outlier_label) {

      #Estimation des valeurs propres selon une procédure de Robbins-Monro

      outlier_labels[i] <- 1


    }
    stat[i] <- resultatOutlier$S

    phat[i] <- resultatOutlier$phat
    #print(phat[i])

  }


  miter[n,] = miter[n-1,]
  VIter[,,n] = VIter[,,n-1]




  return (list(m=m,V=V,lambdatilde = lambdatilde,lambdaIter = lambdaIter,moyennem=moyennem,moyenneV=moyenneV,miter = miter,VIter = VIter,U = U,vpMCM = vpMCM,outlier_labels = outlier_labels, stat = stat,phat = phat))
}


streaming <- function(Y,t,k,c,n,d,q,r)
{
  
  sampsize = d
  
  
  
  #Initialisation de m
  m = r*rnorm(d)
  
  V = diag(d)
  
  
  #Stockage des estimations des valeurs propres de la matrice de covariance
  lambdaIter = matrix(0,n,d)
  #lambda <- eigen(Vvrai)$values
  lambda <- rep(1,d)
  lambdatilde <- rep(1,d)
  
  #vp2 <- eigen(Vvrai)$values
  slog <- 1
  w <- 2
  
  beta = 1/2
  
  #Yv <- mvtnorm::rmvnorm(1000000,sigma=Sigma)
  
  ##Estimation de mn
  
  #Initialisation de m
  
  m <- r*rnorm(d)
  
  
  #Sigma <- diag(d)
  
  
  
  
  moyennem = m
  
  
  
  
  moyenneV <- V
  
  #Stockage des itérations de médiane géométrique 
  
  miter = matrix(0,t,d)
  
  stat <- rep(0,t-1)
  phat <- rep(0,t-1)
  #Stockage des itérations de matrices dans un tableau
  VIter <- array(0, dim = c(d, d, t))
  
  #Stockage des vecteurs propres
  
  U <- array(1, dim = c(t, q, q))
  
  
  #Valeurs propres de V
  vpMCM <- matrix(0,t,d)
  
  
  #Initialisation de U avec des vecteurs propres pas trop éloignés de ceux de cov(X)
  
  U[1,,] <-  1.5*diag(q)
  
  # Vecteur pour stocker les labels des outliers
  outlier_labels <- rep(0, t-1)
  
  phat <- rep(0,t-1)
  
  
  finboucle = t - 1
  
  for (i in (1:finboucle)) 
  {
    
    #print(i)  
    
    #Somme des médianes géométriques
    S <- rep(0,d)
    
    #Somme matrices 
    
    Smatr <- matrix(0,d,d)
    
    gamma = c/i^(0.75) 
    
    #Test pour gérer la non négativité de Vn
    
    if (gamma >= norm((Y[i+1,] - moyennem) %*% t(Y[i+1,] - moyennem) - V,type = "F"))
    {
      gamma = norm((Y[i+1,] - moyennem) %*% t(Y[i+1,] - moyennem) - V,type = "F")
    }
    
    
    for (j in 1:k) 
    {
      
      S <- S + (Y[(i - 1) * k + j, ] - m)/ sqrt(sum((Y[(i - 1) * k + j, ] - m)^2))
    }
    
    #S <- rowSums(sapply(1:k, function(j) (Y[(i - 1) * k + j, ] - m) 
    
    #Mise à jour de m
    
    m <- m + gamma*k^(beta)*S/k
    
    
    moyennem = moyennem*i/(i+1) + 1/(i+1)*m
    
    miter[i,] = moyennem
    
    
    
    
    
    for (j in 1 : k) 
    {
      Smatr <- Smatr +  ((Y[(i - 1) * k + j, ] - moyennem) %*% t(Y[(i - 1) * k + j, ] - m)-V)/norm((Y[(i - 1) * k + j, ] - moyennem) %*% t(Y[(i - 1) * k + j, ] - m)-V,type = "F") 
      
    }
    
    
    V <- V + gamma/k*k^(beta)*Smatr
    
    
    
    
    
    
    
    #Mise à jour de moyenneV
    
    moyenneV <- i/(i+1)*moyenneV + 1/(i + 1)*V
    
    
    #Stockage des itérations de matrices dans un tableau
    VIter[,,i] <- moyenneV
    
    #Estimation des vecteurs propres de Vt et orthonormalisation
    
    for (l in (1:d)) 
    {
      Un =  U[i,,l]/sqrt(sum(U[i,,l]^2))
      
      U[i+1,,l] <-  U[i,,l] + 1/(i + 1)*(moyenneV %*% Un - U[i,,l])
    }
    
    
    #matrice <- U[i+1,,]
    
    #Tri des vecteurs en ordre décroissant de norme
    #normes_colonnes <- apply(matrice, 2, function(col) sqrt(sum(col^2)))
    
    ## order() pour obtenir les indices de tri par ordre décroissant des normes
    #indices_tri <- order(normes_colonnes, decreasing = TRUE)
    
    # Trier la matrice selon les normes décroissantes des colonnes
    #matrice <-matrice[,indices_tri]
    
    #U[i+1,,] <- matrice
    
    
    
    #Orthogonalisation des vecteurs propres
    U[i+1,,] <- orthonormalization(U[i+1,,], basis=TRUE, norm=FALSE)
    
    #Estimation des valeurs propres de la MCM
    #vpMCM[i,] = apply(U[i,,], 2, function(col) sqrt(sum(col^2)))
    vpMCM[i,] = sqrt(colSums(U[i,,]^2))
    #Prendre norme de U[]
    #Z=Y[i,]
    # Z2=X2[i,]
    #E1=    Z^2*(sum(( (vp)-lambda*(Z^2))^2) + sum((lambda * Z^2)%*%t((lambda * Z^2))) - sum(((lambda * Z^2)^2))  )^(-0.5)
    #slog=slog+log(i+1)^w
    #vp2=vp2+log(i+1)^w *((slog)^(-1)) *(lambda - vp2)
    #vp0 = vp2
    #E2=    (sum(( (vp)-lambda*(Z^2))^2) + sum((lambda * Z^2)%*%t((lambda * Z^2))) - sum(((lambda * Z^2)^2))  )^(-0.5)
    #lambda = lambda  - c*i^(-0.75)*lambda*E1 + c*i^(-0.75)* (vp)*E2
    #Convergence si initialisation à vpMCM[i,] (= delta) au lieu de lambda
    lambdaResultat <- RobbinsMC2(sampsize,vp=vpMCM[i,],samp=1:sampsize,init = lambdatilde,initbarre = lambda,ctilde = sampsize*(i-1),cbarre = sampsize*(i-1))
    lambda <- lambdaResultat$vp
    lambdatilde <- lambdaResultat$lambda
    #print(lambda)
    lambdaIter[i,] <- lambdatilde
    
    #Rédiger partie simulation
    #vp  = lambdaEstim (précédente valeur de lambda) sample 100
    
    #Récupération des résultats de la fonction Outlier
    #resultatOutlier <- Outlier(donnee = Y[i, ],seuil_p_value = 0.05,VP = U[i,,],m = moyennem,lambdatilde)
    resultatOutlier <- Outlier(donnee = Y[i, ],seuil_p_value = 0.05,VP = eigen(V)$vectors,m = moyennem,lambdatilde)
    
    # Vérifier si l'entrée est un outlier
    if (resultatOutlier$outlier_label) {
      
      #Estimation des valeurs propres selon une procédure de Robbins-Monro
      
      outlier_labels[i] <- 1
      
      
    }
    stat[i] <- resultatOutlier$S
    
    phat[i] <- resultatOutlier$phat
    #print(phat[i])
    

  }  
  
  
  return (list(m=m,V=V,lambdatilde = lambdatilde,lambdaIter = lambdaIter,moyennem=moyennem,moyenneV=moyenneV,miter = miter,VIter = VIter,U = U,vpMCM = vpMCM,outlier_labels = outlier_labels, stat = stat,phat = phat))
}



