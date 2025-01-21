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
estimMV <- function(Y,c, exposantPas = 0.75,aa = 1,r = 1.5, 
                    minit = r*rnorm(ncol(Y)),Vinit = diag(ncol(Y))
                    ,U = array(1, dim = c(nrow(Y), ncol(Y),ncol(Y))),
                    vpMCM = matrix(0,n,ncol(Y))
                    ,methode = "eigen",depart = 0,cutoff =qchisq(p = 0.95, df = d),niterRMon = ncol(Y))
{

  #Initialisations 
  Sigma = array(0, dim = c(nrow(Y), ncol(Y),ncol(Y)))
  sampsize = ncol(Y)
  lambda = rep(1,ncol(Y))
  lambdatilde = rep(1,ncol(Y))
  
  #Stockage des estimations des valeurs propres de la matrice de covariance
  lambdaIter = matrix(0,nrow(Y),ncol(Y))

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

  mvrai <- rep(0,ncol(Y))

  #Stockage des distances
  
  distances <- rep(0,nrow(Y))
  

  #Calcul de la vraie MCM par l'algorithme de Weiszfeld

  moyennem = m

  #VvectproprCovV <- VcovVrai$vectors

  moyenneV <- V

  #Stockage des itérations de matrices dans un tableau
  VIter <- array(0, dim = c(ncol(Y), ncol(Y), nrow(Y)))

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

  
  matrix_random <- matrix(0, ncol(Y), ncol(Y))

  # Générer chaque colonne aléatoire sur la sphère unité
  for (i in (1:ncol(Y))) {
    v <- rnorm(ncol(Y))  # Tirer un vecteur d composantes normales
    matrix_random[, i] <- v / sqrt(sum(v^2))  # Normaliser
  }

  U[1,,] <- matrix_random

}

  #Stockage des itérations
  miter = matrix(0,nrow(Y),ncol(Y))

  # Vecteur pour stocker les labels des outliers
  outlier_labels <- rep(0, nrow(Y)-1)

  phat <- rep(0,nrow(Y)-1)

  for (i  in 1:(nrow(Y)-1))
  {

    gamma = c/(i+depart)^(exposantPas)
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

    for (l in (1:ncol(Y)))
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
    VPropresV <- VPropresV %*% diag(1/sqrt(colSums(VPropresV^2)))

    valPV <- eigen(V)$values
    lambdaResultat <- RobbinsMC2(niterRMon,vp=valPV,samp=1:sampsize,init = lambdatilde,initbarre = lambda,ctilde = sampsize*(i-1),cbarre = sampsize*(i-1))
    lambda <- lambdaResultat$vp
    lambdatilde <- lambdaResultat$lambda
    #print(lambda)
    lambdaIter[i,] <- lambdatilde
    #Calcul de Sigma
    #Sigma <- VPropresV %*% diag(valPV) %*% t(VPropresV)
    #distances <- rep(0,nrow(Z))
    #distances <- calcule_vecteur_distances(Z,m,Sigma)
    #cutoff <- calcule_cutoff(distances,d)
    #resultatOutlier <- Outlier(donnee = Y[i, ],seuil_p_value = 0.05,VP = eigen(V)$vectors,m = moyennem,lambdatilde)

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
      lambdaResultat <- RobbinsMC2(niterRMon,vp=vpMCM[i,],samp=1:sampsize,init = lambdatilde,initbarre = lambda,ctilde = sampsize*(i-1),cbarre = sampsize*(i-1))
      lambda <- lambdaResultat$vp
      lambdatilde <- lambdaResultat$lambda
      #print(lambda)
      lambdaIter[i,] <- lambdatilde


      #VPropresV <- U[i,,] %*% diag(1/sqrt(colSums(U[i,,]^2)))

      # for (j in (1:d))
      #  {
      #   SigmaEstim <- 1/lambda[i,j]*U[i,,j]%*%t(U[i,,j])
      #}
      #Sigma <- (VPropresV) %*% diag(lambdatilde)%*% t(VPropresV)

      #distances <- rep(0,nrow(Z))
      #distances <- calcule_vecteur_distances(Z,m, Sigma)
      #cutoff <- calcule_cutoff(distances,d)


      #resultatOutlier <- Outlier(donnee = Y[i, ],seuil_p_value = 0.05,VP = U[i,,],m = moyennem,lambdatilde)
      }
    # Vérifier si l'entrée est un outlier
    #if (resultatOutlier$outlier_label) {

      #Estimation des valeurs propres selon une procédure de Robbins-Monro

     # outlier_labels[i] <- 1


    #}
    #stat[i] <- resultatOutlier$S

    #phat[i] <- resultatOutlier$phat
    #print(phat[i])
#Estimation de Sigma 
    VP <- U[i,,] %*% diag(1/sqrt(colSums(U[i,,]^2)))
    Sigma[i,,] <-  VP%*% diag(lambdaIter[i,]) %*% t(VP) 
  }


  miter[nrow(Y),] = miter[nrow(Y)-1,]
  VIter[,,nrow(Y)] = VIter[,,nrow(Y)-1]




  return (list(m=m,V=V,lambdatilde = lambdatilde,lambdaIter = lambdaIter,moyennem=moyennem,moyenneV=moyenneV,miter = miter,VIter = VIter,U = U,vpMCM = vpMCM,outlier_labels = outlier_labels,distances = distances, Sigma = Sigma))
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
    #if (resultatOutlier$outlier_label) {

      #Estimation des valeurs propres selon une procédure de Robbins-Monro

      #outlier_labels[i] <- 1


#    }
 #   stat[i] <- resultatOutlier$S

  #  phat[i] <- resultatOutlier$phat
    #print(phat[i])


  }


  return (list(m=m,V=V,lambdatilde = lambdatilde,lambdaIter = lambdaIter,moyennem=moyennem,moyenneV=moyenneV,miter = miter,VIter = VIter,U = U,vpMCM = vpMCM))
}


#Estimation de la matrice de covariance par shrinkage

covCor <- function(Y, k = -1) {
  dim.Y <- dim(Y)
  N <- dim.Y[1]
  p <- dim.Y[2]
  if (k < 0) {    # demean the data and set k = 1
    Y <- scale(Y, scale = F)
    k <- 1
  }
  n <- N - k    # effective sample size
  c <- p / n    # concentration ratio
  #sample <- (t(Y) %*% Y) / n   
  sample <- covComed(Y)$cov
  # compute shrinkage target
  samplevar <- diag(sample)
  sqrtvar <- sqrt(samplevar)
  rBar <- (sum(sample / outer(sqrtvar, sqrtvar)) - p) / (p * (p - 1))
  target <- rBar * outer(sqrtvar, sqrtvar)
  diag(target) <- samplevar
  
  # estimate the parameter that we call pi in Ledoit and Wolf (2003, JEF)
  Y2 <- Y^2
  sample2 <- (t(Y2) %*% Y2) / n   
  #sample2 <- covComed(Y2)
  piMat <- sample2 - sample^2
  pihat <- sum(piMat)
  
  # estimate the parameter that we call gamma in Ledoit and Wolf (2003, JEF)
  gammahat <- norm(c(sample - target), type = "2")^2
  
  # diagonal part of the parameter that we call rho 
  rho_diag <- sum(diag(piMat))
  
  # off-diagonal part of the parameter that we call rho 
  term1 <- (t(Y^3) %*% Y) / n;
  term2 <- rep.row(samplevar, p) * sample;
  term2 <- t(term2)
  thetaMat <- term1 - term2
  diag(thetaMat) <- 0
  rho_off <- rBar * sum(outer(1/sqrtvar, sqrtvar) * thetaMat)
  
  # compute shrinkage intensity
  rhohat <- rho_diag + rho_off
  kappahat <- (pihat - rhohat) / gammahat
  shrinkage <- max(0, min(1, kappahat / n))
  
  # compute shrinkage estimator
  sigmahat <- shrinkage * target + (1 - shrinkage) * sample
  
  
  return (sigmahat)
  
}

rep.row <- function(x, n){
  matrix(rep(x, each = n), nrow = n)
}

