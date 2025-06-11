#install.packages("microbenchmark")
#install.packages("matlib")
#install.packages("MASS")
#install.packages("reshape2")
#install.packages("corrplot")
#install.packages("dlpyr")


RobbinsMC2=function(mc_sample_size=10000,vp,epsilon=10^(-8),alpha=0.75,c=length(vp),w=2,samp=mc_sample_size,init=vp,initbarre = vp,cbarre = 0,ctilde = 0,slog=1)
{
  p=length(vp)
  vp2=initbarre
  lambda=init
  lambdalist=c()
  vplist=c()
  Y=matrix(rnorm(mc_sample_size*p),ncol=p)
  # X2=matrix(rnorm(mc_sample_size*p),ncol=p)
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
    #print(i)
    #print(" ")
    #print(vp2)
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

#Estimation online
estimMVOnline <- function(Y,c = sqrt(ncol(Y)), exposantPas = 0.75,aa = 1,r = 1.5, w=2,cMC=ncol(Y),
                    minit = r*rnorm(ncol(Y)),Vinit = diag(d)
                    ,U = array(1, dim = c(nrow(Y), ncol(Y),ncol(Y))),
                    vpMCM = matrix(0,ncol(Y),ncol(Y)),lambdaInit =  rep(1,ncol(Y))
                    ,SigmaInit = diag(d),methode = "eigen",depart = 100,cutoff =qchisq(p = 0.95, df = ncol(Y)),niterRMon = nrow(Y)*ncol(Y))
{
  
  lambdatilde = rep(1,ncol(Y))
  lambdaIter = matrix(0,nrow(Y),ncol(Y))
  U = array(0, dim = c(nrow(Y), ncol(Y),ncol(Y)))
  Sigma = array(0, dim = c(nrow(Y), ncol(Y),ncol(Y)))
  distances <- rep(0,nrow(Y))
  #Initialisations
  m <- minit
  V <- Vinit
  outlier_labels = rep(0,nrow(Y))
  
  #Si départ = 0, intialisation de U sur la sphère unité  
  if (depart == 0)
  {
    
    matrix_random <- matrix(0, ncol(Y), ncol(Y))
    
    #Générer chaque colonne aléatoire sur la sphère unité
    for (i in (1:ncol(Y))) {
      v <- rnorm(ncol(Y))  # Tirer un vecteur d composantes normales
      matrix_random[, i] <- v / sqrt(sum(v^2))  # Normaliser
    }
    
    U[1,,] <- matrix_random
    m <- minit
    V <- Vinit
    moyennem <- m
    
  }
  if (depart > 0)
  {
    #resoff=RobVar(Y[1:depart,],mc_sample_size = nrow(Y),c=ncol(Y),w=2)
    # resoff = WeiszfeldCov_init(Y[1:depart,],init = r*rnorm(ncol(Y[1:depart,])),init_cov = covComed(Y[1:depart,])$cov)
    # med <- resoff$median
    # V <- WeiszfeldCov_init(Y[1:depart,],init = med,init_cov = covComed(Y[1:depart,])$cov)$covmedian
    # eig_init = eigen(V)
    # #eig_init=eigen(resoff$variance)
    # valPV <- eig_init$values 
    # #print(valPV)
    # valPV = apply(cbind(valPV,rep(10^(-4),length(valPV))),MARGIN=1, FUN=max)
    # #print(valPV)
    # #lambdaInit <- RobbinsMC2(c=cMC,mc_sample_size = niterRMon,w=w,vp=valPV,samp=1:niterRMon,init = valPV,initbarre = valPV,ctilde = niterRMon*(niterr-1),cbarre =niterRMon*(niterr-1),slog=sum((log(1:((niterRMon*(niterr-1))+1))^w)))
    # #lambdaResultat <- RobbinsMC2(c=cMC,mc_sample_size = niterRMon,w=w,vp=valPV,samp=1:niterRMon,init = valPV,initbarre = valPV,ctilde = niterRMon,cbarre =niterRMon)
    # #lambdaInit <- lambdaResultat$vp
    # #lambdaInit=eig_init$values
    # #lambdatilde=lambdaInit
    # lambdaIter[1:depart,]=matrix(rep(lambdaInit,depart),byrow=T,nrow=depart)
    # 
    
    resoff = WeiszfeldCov_init(Y[1:depart,],r*rnorm(ncol(Y)),init_cov = covComed(Y[1:100,])$cov,nitermax = 100)
    minit <- resoff$median
    Vinit <- WeiszfeldCov_init(Y[1:depart,],minit,init_cov = covComed(Y[1:100,])$cov,nitermax = 100)$covmedian
    #V <-  GmedianCov(Y, init = med,scores = ncol(Y))$covmedian
    #eig_init = eigen(Vinit)
    eig_init = eigs_sym(Vinit,k = ncol(Y))
    #eig_init = eigen( WeiszfeldCov(Y, nitermax = 1000)$covmedian)
    #eig_init=eigen(resoff$variance)
    valPV <- eig_init$values 
    #print(valPV)
    valPV = apply(cbind(valPV,rep(10^(-4),length(valPV))),MARGIN=1, FUN=max)
    #print(valPV)
    #lambdaInit <- RobbinsMC2(c=cMC,mc_sample_size = niterRMon,w=w,vp=valPV,samp=1:niterRMon,init = valPV,initbarre = valPV,ctilde = niterRMon*(niterr-1),cbarre =niterRMon*(niterr-1),slog=sum((log(1:((niterRMon*(niterr-1))+1))^w)))
    lambdaInit = valPV
    lambdatilde = valPV
    
    for (i in 1:nrow(Y[1:100,]))
      
    {
      sampsize = ncol(Y)
      
      lambda = lambdaInit
      lambdatilde = lambdatilde
      
      lambdaResultat <- RobbinsMC2(sampsize,c = cMC, vp=valPV,w=w,samp=1:sampsize,init = lambdatilde,initbarre = lambda,ctilde = sampsize*(i-1),cbarre =sampsize*(i-1),slog=sum((log(1:((sampsize*(i-1))+1))^w)))
      #ctilde = sampsize*(i-1),cbarre =sampsize*(i-1)
      lambda <- lambdaResultat$vp
      lambdatilde <- lambdaResultat$vp
      
    }
    lambdaInit <- lambda
    
    #lambdaResultat <- RobbinsMC2(c=cMC,mc_sample_size = niterRMon,w=w,vp=valPV,samp=1:niterRMon,init = valPV,initbarre = valPV,ctilde = 0,cbarre =0)
    #lambdaInit <- lambdaResultat$vp
    #lambdaInit=eig_init$values
    #lambdatilde=lambdaInit
    #lambdaIter[1:depart,]=matrix(rep(lambdaInit,depart),byrow=T,nrow=depart)
    
    #m <- minit
    
    #V <- Vinit
    
    #moyennem <- m
    #moyenneV <- V
    VPropresV <- eig_init$vectors
    
    VP <- VPropresV %*% diag(1/sqrt(colSums(VPropresV^2)))
    #print(lambdaIter[i,])
    #print(lambdaInit)
    varianc= VP %*% diag(lambdaInit) %*% t(VP) 
    
    
    #SigmaOffline <- varianc
    #V <- resoff$covmedian
    #U <- eigen(V)$vectors
    #m <- minit
    
    #V <- Vinit
    
    moyennem <- minit
    moyenneV <- Vinit
    VPropresV <- eig_init$vectors
    VP <- VPropresV %*% diag(1/sqrt(colSums(VPropresV^2)))
    #print(lambdaIter[i,])
    #print(lambdaInit)
    varianc= VP %*% diag(lambdaInit) %*% t(VP) 
    
    for (l in 1 :depart)
    {
      U[l,,] <- eig_init$vectors
      
      #Faire calcul t(P) %*% diag(D) %*% P
      #Sigma[l,,] <- resoff$variance
      Sigma[l,,] <- varianc
      S <- 0
      for(j in (1:ncol(Y)))
      {
        S <- S + 1/(lambdaInit[j])*sum((Y[l,] - moyennem)*(eig_init$vectors)[,j])^2
      }
      
      #print(solve(Sigma[i,,]))
      #distances[l] <- as.numeric(Y[l,] - m) %*% solve(Sigma[l,,]) %*% (as.numeric(t(Y[l,] - m)))
      distances[l] <- S
      
    }
  if (distances[l] > cutoff) {outlier_labels[l] = 1}
    }
  
  
  sampsize = niterRMon
  lambda = lambdaInit
  lambdatilde = lambdatilde
  
  
  #m <- minit
  
  #V <- Vinit
  
  
  #Calcul de la vraie MCM par l'algorithme de Weiszfeld
  
  #moyennem = m
  
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
  #Compteur pour les non outliers
  
  #j <- 0
  
  
  #Stockage des vecteurs propres dans un tableau
  
  
  #Stockage des itérations
  miter = matrix(0,nrow(Y),ncol(Y))
  
  # Vecteur pour stocker les labels des outliers
  #outlier_labels <- rep(0, nrow(Y)-1)
  
  phat <- rep(0,nrow(Y)-1)
  
  for (i  in (depart+1):(nrow(Y)-1))
  {
    
    gamma = c/(i)^(exposantPas)
    #gamma = c/i
    
    ##Test pour gérer la non négativité de Vn
    
    
    if (gamma >= norm(c(Y[i+1,] - moyennem) %*% t(c((Y[i+1,] - moyennem))) - V,type = "F"))
    {
      gamma = norm(c(Y[i+1,] - moyennem) %*% t(c((Y[i+1,] - moyennem))) - V,type = "F")
    }
    
    
    
    
    #Mise à jour de mn
    
    m = m + gamma*((Y[i+1,] - m))/sqrt(sum((Y[i+1,]  - m)^2))
    
    
    
    
    
    #Mise à jour de V
    
    V = V + gamma*(c(Y[i+1,] - moyennem) %*% t(c((Y[i+1,] - moyennem))) - V)/norm(c(Y[i+1,] - moyennem) %*% t(c((Y[i+1,] - moyennem)))- V,type = "F")
    
    
    moyennem = moyennem*(i - depart)/(i - depart + 1) + 1/(i - depart +1)*m
    
    miter[i,] = moyennem
    
    #Mise à jour de moyenneV
    
    moyenneV <- (i - depart)/(i - depart  +1)*moyenneV + 1/(i- depart  + 1)*V
    
    
    VIter[,,i] = moyenneV
    #U[i,,] <- orthonormalization(U[i,,],basis = TRUE, norm = TRUE)
    #Calcul des phijn
    #for (l in (1:q))
    #{
    # phijn[i,l] <- sum(t(Y[i+1,] - m)*Uphi[i,,l])
    #print(phijn[i,l])
    #}
    
    #vp = apply(U[i+1,,],2,norm)
    
    #Estimation des valeurs propres de la MCM
    #vpMCM[i,] = apply(U[i,,], 2, function(col) sqrt(sum(col^2)))
    
    #Rédiger partie simulation
    #vp  = lambdaEstim (précédente valeur de lambda) sample 100
    
    #Récupération des résultats de la fonction Outlier
    #resultatOutlier <- Outlier(donnee = Y[i, ],seuil_p_value = 0.05,VP = U[i,,],m = moyennem,lambdatilde)
    if(methode == "eigen"){
      elpropresV = eigs_sym(moyenneV,k = ncol(Y))
      VPropresV <- elpropresV$vectors
      #VPropresV <- VPropresV %*% diag(1/sqrt(colSums(VPropresV^2)))
      
      valPV <- elpropresV$values  
      lambdaResultat <- RobbinsMC2(niterRMon,c = cMC, vp=valPV,w=w,samp=1:sampsize,init = lambdatilde,initbarre = lambda,ctilde = sampsize*(i-1),cbarre =sampsize*(i-1),slog=sum((log(1:((sampsize*(i-1))+1))^w)))
      #ctilde = sampsize*(i-1),cbarre =sampsize*(i-1)
      lambda <- lambdaResultat$vp
      lambdatilde <- lambdaResultat$lambda
      #print(lambda)
      #lambdaIter[i,] <- lambdatilde
      lambdaIter[i,] <- lambda
      #Calcul de Sigma
      #Sigma <- VPropresV %*% diag(valPV) %*% t(VPropresV)
      #distances <- rep(0,nrow(Z))
      #distances <- calcule_vecteur_distances(Z,m,Sigma)
      #cutoff <- calcule_cutoff(distances,d)
      #resultatOutlier <- Outlier(donnee = Y[i, ],seuil_p_value = 0.05,VP = eigen(V)$vectors,m = moyennem,lambdatilde)
      VP <- VPropresV %*% diag(1/sqrt(colSums(VPropresV^2)))
      #print(lambdaIter[i,])
      Sigma[i,,] <-  VP %*% diag(lambdaIter[i,]) %*% t(VP) 
      #print(Sigma[i,,])
      #Calcul de la distance i+1 sur la donnée i+1
      S <- 0
      for(j in (1:ncol(Y)))
      {
        
        S<- S + 1/(lambdaIter[i,j])*sum(t(Y[i+1,] - moyennem)*(VP[,j]))^2
      #print(solve(Sigma[i,,]))
      #distances[i] <- as.numeric(Y[i,] - m) %*% solve(Sigma[i,,]) %*% (as.numeric(t(Y[i,] - m)))
      
      }      
      distances[i+1] <- S
      if (distances[i+1] > cutoff) {outlier_labels[i+1] = 1}
      
      }
    else if(methode == "CPP"){
      elPropres <- calculValeursEtVecteursPropres(moyenneV)
      VPropresV <- elPropres$vecteurs_propres
      #VPropresV <- VPropresV %*% diag(1/sqrt(colSums(VPropresV^2)))
      
      valPV <- elPropres$valeurs_propres  
      lambdaResultat <- RobbinsMC2(niterRMon,c = cMC, vp=valPV,w=w,samp=1:sampsize,init = lambdatilde,initbarre = lambda,ctilde = sampsize*(i-1),cbarre =sampsize*(i-1),slog=sum((log(1:((sampsize*(i-1))+1))^w)))
      #ctilde = sampsize*(i-1),cbarre =sampsize*(i-1)
      lambda <- lambdaResultat$vp
      lambdatilde <- lambdaResultat$lambda
      #print(lambda)
      #lambdaIter[i,] <- lambdatilde
      lambdaIter[i,] <- lambda
      #Calcul de Sigma
      #Sigma <- VPropresV %*% diag(valPV) %*% t(VPropresV)
      #distances <- rep(0,nrow(Z))
      #distances <- calcule_vecteur_distances(Z,m,Sigma)
      #cutoff <- calcule_cutoff(distances,d)
      #resultatOutlier <- Outlier(donnee = Y[i, ],seuil_p_value = 0.05,VP = eigen(V)$vectors,m = moyennem,lambdatilde)
      VP <- VPropresV %*% diag(1/sqrt(colSums(VPropresV^2)))
      #print(lambdaIter[i,])
      Sigma[i,,] <-  VP %*% diag(lambdaIter[i,]) %*% t(VP) 
      #print(Sigma[i,,])
      #Calcul de la distance i
      S <- 0
      for(j in (1:ncol(Y)))
      {
        S<- S + 1/(lambdaIter[i,j])*sum(t(Y[i,] - moyennem)*(VP[,j]))^2
        
        #print(solve(Sigma[i,,]))
        #distances[i] <- as.numeric(Y[i,] - m) %*% solve(Sigma[i,,]) %*% (as.numeric(t(Y[i,] - m)))
       
        
      }
      distances[i] <- S
      if(distances[i] > cutoff) {outlier_labels[i] = 1}
      }
    else {
      
      print("méthode sans eigen")
      
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
        
        #U[i+1,,] <- matrice
        #Orthogonalisation des vecteurs propres
        U[i+1,,] <- orthonormalization(U[i+1,,],basis = TRUE, norm = FALSE)
        
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
        lambdaResultat <- RobbinsMC2(c = cMC,niterRMon,vp=vpMCM[i,],samp=1:sampsize,init = lambdatilde,initbarre = lambda,ctilde = sampsize*(i-1),cbarre = sampsize*(i-1),slog=sum((log(1:((sampsize*(i-1))+1))^w)))
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
      #distances[i] <- as.numeric(Y[i,] - m) %*% solve(Sigma[i,,]) %*% (as.numeric(t(Y[i,] - m)))
      #Calcul de la distance i
      S <- 0
      for(j in (1:ncol(Y)))
      {
        S <- S + 1/(lambdaIter[i,j])*(t(Y[i,] - moyennem)%*%VP[,j])^2
      }
      
      #print(solve(Sigma[i,,]))
      #distances[i] <- as.numeric(Y[i,] - m) %*% solve(Sigma[i,,]) %*% (as.numeric(t(Y[i,] - m)))
      distances[i] <- S
      
    }
    
    
    
    #matrice <- U[i+1,,]
    
    #Tri des vecteurs en ordre décroissant de norme
    #normes_colonnes <- apply(matrice, 2, function(col) sqrt(sum(col^2)))
    
    # order() pour obtenir les indices de tri par ordre décroissant des normes
    #indices_tri <- order(normes_colonnes, decreasing = TRUE)
    
    # Trier la matrice selon les normes décroissantes des colonnes
    #matrice <-matrice[,indices_tri]
    
  }
  
  
  miter[nrow(Y),] = miter[nrow(Y)-1,]
  VIter[,,nrow(Y)] = VIter[,,nrow(Y)-1]
  Sigma[nrow(Y),,] = Sigma[(nrow(Y) - 1),,]
  
  
  
  return (list(m=m,V=V,lambdatilde = lambdatilde,lambdaIter = lambdaIter,moyennem=moyennem,moyenneV=moyenneV,miter = miter,VIter = VIter,U = U,vpMCM = vpMCM,outlier_labels = outlier_labels,distances = distances, Sigma = Sigma))


}

#Fonction qui prend en paramètre une matrice de données, et une méthode d'estimation et renvoie les paramètres estimés, les distances et les outliers
estimation <- function(Y,c = ncol(Y), exposantPas = 0.75,aa = 1,r = 1.5,sampsize = ncol(Y),  cMC=ncol(Y),w= 2,
                      minit = r*rnorm(ncol(Y)),Vinit = diag(ncol(Y))
                      ,U = array(1, dim = c(nrow(Y), ncol(Y),ncol(Y))),
                      vpMCM = matrix(0,n,ncol(Y)),methodeOnline = "eigen", methodeEstimation = "offline",depart_online = 100,niterRMon = nrow(Y)*ncol(Y))
{
  distances <- rep(0, nrow(Y))
  if (methodeEstimation == "offline")
  {
    #print(Y)
    #results <- RobVar(Y,mc_sample_size = nrow(Y)*ncol(Y),c=ncol(Y),w=2)
    
    #Récupération de la médiane de Sigma des vecteurs propres (variable U), et des valeurs propres
    #resoff=RobVar(Y,mc_sample_size = nrow(Y)*ncol(Y),c=ncol(Y),w=2)
    
    resoff = WeiszfeldCov_init(Y,r*rnorm(ncol(Y)),init_cov = covComed(Y)$cov,nitermax = 100)
    med <- resoff$median
    V <- WeiszfeldCov_init(Y,med,init_cov = covComed(Y)$cov,nitermax = 100)$covmedian
    #V <-  GmedianCov(Y, init = med,scores = ncol(Y))$covmedian
    eig_init = eigs_sym(V,k = ncol(Y))
    #eig_init = eigen( WeiszfeldCov(Y, nitermax = 1000)$covmedian)
    #eig_init=eigen(resoff$variance)
    valPV <- eig_init$values 
    #print(valPV)
    valPV = apply(cbind(valPV,rep(10^(-4),length(valPV))),MARGIN=1, FUN=max)
    #print(valPV)lam
    #lambdaInit <- RobbinsMC2(c=cMC,mc_sample_size = niterRMon,w=w,vp=valPV,samp=1:niterRMon,init = valPV,initbarre = valPV,ctilde = niterRMon*(niterr-1),cbarre =niterRMon*(niterr-1),slog=sum((log(1:((niterRMon*(niterr-1))+1))^w)))
    lambdaInit = valPV
    lambdatilde = valPV
    distances <- rep(0,nrow(Y))
    lambda = lambdaInit
    lambdatilde = lambdatilde
    VPropresV <- eig_init$vectors
    
    VP <- VPropresV %*% diag(1/sqrt(colSums(VPropresV^2)))
    #print(lambdaIter[i,])
    #varianc= VP %*% diag(lambdaInit) %*% t(VP) 
    
    
    #SigmaOffline <- varianc
    V <- resoff$covmedian
    U <- eigen(V)$vectors
    #lambda <- RobbinsMC()
    #distances <- calcule_vecteur_distances(Y, med, SigmaOffline) 
    for (i in 1:nrow(Y))
      
    {
      mc_sample_size = ncol(Y)
      UU=matrix(rnorm(mc_sample_size*ncol(Y)),ncol=ncol(Y))
      lambdaResultat <-   robbinsMC(U=UU, c=cMC, w=w,delta=valPV,init = lambdatilde,
                                    init_bar = lambda,c_tilde = mc_sample_size*(i-1),
                                    c_bar =mc_sample_size*(i-1), sumlog=sum((log(1:((mc_sample_size*(i-1))+1))^w)))
      #lambdaResultat = robbinsMC(U = ,delta = valPV)
      #lambdaResultat <- RobbinsMC2(sampsize,c = cMC, vp=valPV,w=w,samp=1:sampsize,init = lambdatilde,initbarre = lambda,ctilde = sampsize*(i-1),cbarre =sampsize*(i-1),slog=sum((log(1:((sampsize*(i-1))+1))^w)))
      #ctilde = sampsize*(i-1),cbarre =sampsize*(i-1)
      lambda <- lambdaResultat$vp
      lambdatilde <- lambdaResultat$vp
      
    }
    
    #Calcul des distances 
    for (i in 1:nrow(Y))
      
    {
      S <- 0
      for(j in (1:ncol(Y)))
      {
        S <- S + 1/(lambda[j])*sum(t(Y[i,] - med)*eig_init$vectors[,j])^2
      }
      
      #print(solve(Sigma[i,,]))
      #distances[i] <- as.numeric(Y[i,] - m) %*% solve(Sigma[i,,]) %*% (as.numeric(t(Y[i,] - m)))
      distances[i] <- S
      
    }
    
    varianc= VP %*% diag(as.vector(lambda)) %*% t(VP) 
    
    
    SigmaOffline <- varianc
    
    #lambdaInit <- lambda
  
    #lambdaResultat <- RobbinsMC2(c=cMC,mc_sample_size = niterRMon,w=w,vp=valPV,samp=1:niterRMon,init = valPV,initbarre = valPV,ctilde = 0,cbarre =0)
    #lambdaInit <- lambdaResultat$vp
    #lambdaInit=eig_init$values
    #lambdatilde=lambdaInit
    #lambdaIter[1:depart,]=matrix(rep(lambdaInit,depart),byrow=T,nrow=depart)
    
    #m <- minit
    
    #V <- Vinit
    
    #moyennem <- m
    #moyenneV <- V
   
    
  }
  else if (methodeEstimation == "online")
  { 
    
    if (methodeOnline == "eigen"){
    #results <- estimMVOnline(Y, depart = depart_online,niterRMon = ncol(Y))
    results = StreamingRobustVariance(Y,batch = 1,mc_sample_size = floor(sqrt(ncol(Y))))
      }
    else {results <- estimMVOnline(Y, depart = depart_online,methode = "CPP",niterRMon = ncol(Y))}
    #Retour des résultats
    med <- results$moyennem
    miter <- results$miter
    SigmaOnline <- results$Sigma
    U <- results$U
    lambda <- results$lambda
    V <- results$moyenneV
    distances <- results$distances
    outliers = results$outlier_labels
    # results$Sigma[101,,]
  }
  else if (methodeEstimation == "streaming")
  { 
    #results <- StreamingMV(Y,batch = ncol(Y),depart = depart_online,niterRMon = ncol(Y))
    results = StreamingRobustVariance(Y,batch =ncol(Y) ,mc_sample_size = ncol(Y)*floor(sqrt(ncol(Y))))
    
    #Retour des résultats
    med <- results$moyennem
    SigmaStreaming <- results$Sigma
    U <- results$U
    lambda <- results$lambda
    V <- results$moyenneV
    distances <- results$distances
    miter <- results$miter
    outliers = results$outlier_labels
    
    # results$Sigma[101,,]
  }
  
  if (methodeEstimation  == "online") {
    return(list(med = med, SigmaOnline = SigmaOnline[nrow(Y)- 1,,], SigmaOnlineIter = SigmaOnline,U = U ,lambda = lambda,V = V, distances = distances,outliers = outliers))
}
  else if (methodeEstimation  == "offline"){return(list(med = med, SigmaOffline = SigmaOffline ,V = V, distances = distances))}
else {return(list(med = med, SigmaStreamingIter = SigmaStreaming,SigmaStreaming = SigmaStreaming[nrow(Y)-1,,] ,V = V, distances = distances,outliers = outliers))}
  }

#Estimation streaming
StreamingMV <- function(Y,c = sqrt(ncol(Y)), exposantPas = 0.75,aa = 1,r = 1.5,w=2, cMC=ncol(Y),
                        minit = r*rnorm(ncol(Y)),Vinit = diag(d)
                        ,U = array(1, dim = c(nrow(Y), ncol(Y),ncol(Y))), batch=ncol(Y),
                        vpMCM = matrix(0,ncol(Y),ncol(Y)),lambdaInit =  rep(1,ncol(Y))
                        ,SigmaInit = diag(d),methode = "eigen",depart = 0,cutoff =qchisq(p = 0.95, df = ncol(Y)),niterRMon = batch)
{
  compt=0
  
  niterr=0
  #print(nrow(Y))
  lambdatilde = rep(1,ncol(Y))
  lambdaIter = matrix(0,nrow(Y),ncol(Y))
  U = array(0, dim = c(nrow(Y), ncol(Y),ncol(Y)))
  Sigma = array(0, dim = c(nrow(Y), ncol(Y),ncol(Y)))
  #V <- Vinit
  
  distances <- rep(0,nrow(Y))
  
  outlier_labels = rep(0,nrow(Y))
  #Initialisations 
  #Si départ = 0, intialisation de U sur la sphère unité  
  if (depart == 0)
  {
    
    matrix_random <- matrix(0, ncol(Y), ncol(Y))
    
    #Générer chaque colonne aléatoire sur la sphère unité
    for (i in (1:ncol(Y))) {
      v <- rnorm(ncol(Y))  # Tirer un vecteur d composantes normales
      matrix_random[, i] <- v / sqrt(sum(v^2))  # Normaliser
    }
    
    U[1,,] <- matrix_random
    
  }
  if (depart > 0)
  {
    #resoff=RobVar(Y[1:depart,],mc_sample_size = nrow(Y),c=ncol(Y),w=2)
    # resoff = WeiszfeldCov_init(Y[1:depart,],init = r*rnorm(ncol(Y[1:depart,])),init_cov = covComed(Y[1:depart,])$cov)
    # med <- resoff$median
    # V <- WeiszfeldCov_init(Y[1:depart,],init = med,init_cov = covComed(Y[1:depart,])$cov)$covmedian
    # eig_init = eigen(V)
    # #eig_init=eigen(resoff$variance)
    # valPV <- eig_init$values 
    # #print(valPV)
    # valPV = apply(cbind(valPV,rep(10^(-4),length(valPV))),MARGIN=1, FUN=max)
    # #print(valPV)
    # #lambdaInit <- RobbinsMC2(c=cMC,mc_sample_size = niterRMon,w=w,vp=valPV,samp=1:niterRMon,init = valPV,initbarre = valPV,ctilde = niterRMon*(niterr-1),cbarre =niterRMon*(niterr-1),slog=sum((log(1:((niterRMon*(niterr-1))+1))^w)))
    # lambdaResultat <- RobbinsMC2(c=cMC,mc_sample_size = niterRMon,w=w,vp=valPV,samp=1:niterRMon,init = valPV,initbarre = valPV,ctilde = niterRMon,cbarre =niterRMon)
    # lambdaInit <- lambdaResultat$vp
    # #lambdaInit=eig_init$values
    # #lambdatilde=lambdaInit
    # lambdaIter[1:depart,]=matrix(rep(lambdaInit,depart),byrow=T,nrow=depart)
    # 
    # m <- minit
    # 
    # V <- Vinit
    # 
    # moyennem <- m
    # moyenneV <- V
    # VPropresV <- eig_init$vectors
    # VP <- VPropresV %*% diag(1/sqrt(colSums(VPropresV^2)))
    # #print(lambdaIter[i,])
    # #print(lambdaInit)
    # varianc= VP %*% diag(lambdaInit) %*% t(VP) 
# 
#     resoff = WeiszfeldCov_init(Y[1:100,],r*rnorm(ncol(Y)),init_cov = covComed(Y[1:100,])$cov,nitermax = 100)
#     minit <- resoff$median
#     Vinit <- WeiszfeldCov_init(Y[1:100,],minit,init_cov = covComed(Y[1:100,])$cov,nitermax = 100)$covmedian
#     #V <-  GmedianCov(Y, init = med,scores = ncol(Y))$covmedian
#     eig_init = eigen(Vinit)
#     #eig_init = eigen( WeiszfeldCov(Y, nitermax = 1000)$covmedian)
#     #eig_init=eigen(resoff$variance)
#     valPV <- eig_init$values 
#     #print(valPV)
#     valPV = apply(cbind(valPV,rep(10^(-4),length(valPV))),MARGIN=1, FUN=max)
#     #print(valPV)
#     #lambdaInit <- RobbinsMC2(c=cMC,mc_sample_size = niterRMon,w=w,vp=valPV,samp=1:niterRMon,init = valPV,initbarre = valPV,ctilde = niterRMon*(niterr-1),cbarre =niterRMon*(niterr-1),slog=sum((log(1:((niterRMon*(niterr-1))+1))^w)))
#     lambdaInit = valPV
#     lambdatilde = valPV
#     
#     for (i in 1:nrow(Y[1:100,]))
#       
#     {
#       sampsize = ncol(Y)
#       
#       lambda = lambdaInit
#       lambdatilde = lambdatilde
#       
#       lambdaResultat <- RobbinsMC2(sampsize,c = cMC, vp=valPV,w=w,samp=1:sampsize,init = lambdatilde,initbarre = lambda,ctilde = sampsize*(i-1),cbarre =sampsize*(i-1),slog=sum((log(1:((sampsize*(i-1))+1))^w)))
#       #ctilde = sampsize*(i-1),cbarre =sampsize*(i-1)
#       lambda <- lambdaResultat$vp
#       lambdatilde <- lambdaResultat$vp
#       
#     }
#     lambdaInit <- lambda
#     
#     #lambdaResultat <- RobbinsMC2(c=cMC,mc_sample_size = niterRMon,w=w,vp=valPV,samp=1:niterRMon,init = valPV,initbarre = valPV,ctilde = 0,cbarre =0)
#     #lambdaInit <- lambdaResultat$vp
#     #lambdaInit=eig_init$values
#     #lambdatilde=lambdaInit
#     #lambdaIter[1:depart,]=matrix(rep(lambdaInit,depart),byrow=T,nrow=depart)
#     
#     m <- minit
#     
#     V <- Vinit
#     
#     #moyennem <- m
#     #moyenneV <- V
#     VPropresV <- eig_init$vectors
#     
#     VP <- VPropresV %*% diag(1/sqrt(colSums(VPropresV^2)))
#     #print(lambdaIter[i,])
#     #print(lambdaInit)
#     varianc= VP %*% diag(lambdaInit) %*% t(VP) 
#     
#     
#     #SigmaOffline <- varianc
#     #V <- resoff$covmedian
#     U <- eigen(V)$vectors
#     
#     moyennem <- minit
#     moyenneV <- Vinit
#     VPropresV <- eig_init$vectors
#     VP <- VPropresV %*% diag(1/sqrt(colSums(VPropresV^2)))
#     #print(lambdaIter[i,])
#     #print(lambdaInit)
#     varianc= VP %*% diag(lambdaInit) %*% t(VP) 
#     
#     
#         
#      for (l in 1 :depart)
#     {
#       U[l,,] <- eig_init$vectors
#       
#       #Faire calcul t(P) %*% diag(D) %*% P
#       #Sigma[l,,] <- resoff$variance
#       Sigma[l,,] <- varianc
#       S <- 0
#       for(j in (1:ncol(Y)))
#       {
#         S <- S + 1/(lambdaInit[j])*sum((Y[l,] - moyennem)*(eig_init$vectors)[,j])^2
#       }
#       
#       #print(solve(Sigma[i,,]))
#       #distances[l] <- as.numeric(Y[l,] - m) %*% solve(Sigma[l,,]) %*% (as.numeric(t(Y[l,] - m)))
#       distances[l] <- S
#     }
    #minit=resoff$median
    #Vinit=resoff$covmedian

    resoff = WeiszfeldCov_init(Y[1:100,],r*rnorm(ncol(Y)),init_cov = covComed(Y[1:100,])$cov,nitermax = 100)
    minit <- resoff$median
    Vinit <- WeiszfeldCov_init(Y[1:100,],minit,init_cov = covComed(Y[1:100,])$cov,nitermax = 100)$covmedian
    #V <-  GmedianCov(Y, init = med,scores = ncol(Y))$covmedian
    #eig_init = eigen(Vinit)
    eig_init = eigs_sym(Vinit,k = ncol(Y))
    #eig_init = eigen( WeiszfeldCov(Y, nitermax = 1000)$covmedian)
    #eig_init=eigen(resoff$variance)
    valPV <- eig_init$values 
    #print(valPV)
    valPV = apply(cbind(valPV,rep(10^(-4),length(valPV))),MARGIN=1, FUN=max)
    #print(valPV)
    #lambdaInit <- RobbinsMC2(c=cMC,mc_sample_size = niterRMon,w=w,vp=valPV,samp=1:niterRMon,init = valPV,initbarre = valPV,ctilde = niterRMon*(niterr-1),cbarre =niterRMon*(niterr-1),slog=sum((log(1:((niterRMon*(niterr-1))+1))^w)))
    lambdaInit = valPV
    lambdatilde = valPV
    
    for (i in 1:nrow(Y[1:100,]))
      
    {
      sampsize = ncol(Y)
      
      lambda = lambdaInit
      lambdatilde = lambdatilde
      
      lambdaResultat <- RobbinsMC2(sampsize,c = cMC, vp=lambda,w=w,samp=1:sampsize,init = lambdatilde,initbarre = lambda,ctilde = sampsize*(i-1),cbarre =sampsize*(i-1),slog=sum((log(1:((sampsize*(i-1))+1))^w)))
      #ctilde = sampsize*(i-1),cbarre =sampsize*(i-1)
      lambda <- lambdaResultat$vp
      lambdatilde <- lambdaResultat$vp
      
    }
    lambdaInit <- lambda
    
    #lambdaResultat <- RobbinsMC2(c=cMC,mc_sample_size = niterRMon,w=w,vp=valPV,samp=1:niterRMon,init = valPV,initbarre = valPV,ctilde = 0,cbarre =0)
    #lambdaInit <- lambdaResultat$vp
    #lambdaInit=eig_init$values
    #lambdatilde=lambdaInit
    #lambdaIter[1:depart,]=matrix(rep(lambdaInit,depart),byrow=T,nrow=depart)
    
    #m <- minit
    
    #V <- Vinit
    
    #moyennem <- m
    #moyenneV <- V
    VPropresV <- eig_init$vectors
    
    VP <- VPropresV %*% diag(1/sqrt(colSums(VPropresV^2)))
    #print(lambdaIter[i,])
    #print(lambdaInit)
    varianc= VP %*% diag(lambdaInit) %*% t(VP) 
    
    
    #SigmaOffline <- varianc
    #V <- resoff$covmedian
    #U <- eigen(V)$vectors
    m <- minit
    
    V <- Vinit
    
    moyennem <- minit
    moyenneV <- Vinit
    VPropresV <- eig_init$vectors
    VP <- VPropresV %*% diag(1/sqrt(colSums(VPropresV^2)))
    #print(lambdaIter[i,])
    #print(lambdaInit)
    varianc= VP %*% diag(lambdaInit) %*% t(VP) 
    
    for (l in 1 :depart)
    {
      U[l,,] <- eig_init$vectors
      
      #Faire calcul t(P) %*% diag(D) %*% P
      #Sigma[l,,] <- resoff$variance
      Sigma[l,,] <- varianc
      S <- 0
      for(j in (1:ncol(Y)))
      {
        S <- S + 1/(lambdaInit[j])*sum((Y[l,] - moyennem)*eig_init$vectors[,j])^2
      }
      
      #print(solve(Sigma[i,,]))
      #distances[l] <- as.numeric(Y[l,] - m) %*% solve(Sigma[l,,]) %*% (as.numeric(t(Y[l,] - m)))
      distances[l] <- S
      
      if (S > cutoff) {outlier_labels[l] = 1}
    
    
    }
      }
  
  #print(V)
  sampsize = niterRMon
  lambda = lambdaInit
  lambdatilde = lambdatilde
  
  
  #m <- minit
  
  #V <- Vinit
  
  #Calcul de la vraie MCM par l'algorithme de Weiszfeld
  
  #moyennem = m
  
  #VvectproprCovV <- VcovVrai$vectors
  
#  moyenneV <- V
  
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
  #Compteur pour les non outliers
  
  #j <- 0
  sslog = 1
  
  #Stockage des vecteurs propres dans un tableau
  
  
  #Stockage des itérations
  miter = matrix(0,nrow(Y),ncol(Y))
  
  # Vecteur pour stocker les labels des outliers
  outlier_labels <- rep(0, nrow(Y)-1)
  
  phat <- rep(0,nrow(Y)-1)
  
  for (i  in (depart+1):(nrow(Y)-1))
  {
    niterr=niterr+1
    if (niterr*batch+depart > nrow(Y)){
      break
    }
    gamma = c*sqrt(batch)/(i)^(exposantPas)
    #gamma = c/i
    
    ##Test pour gérer la non négativité de Vn
    gradV=0
    normm=0
    for (l in 1:batch)
    {
      normm=normm+norm(c(Y[depart + (niterr-1)*batch + l,] - moyennem) %*% t(c((Y[depart + (niterr-1)*batch + l,] - moyennem))) - V,type = "F")^(-1)
      gradV=gradV+(c(Y[depart + (niterr-1)*batch + l,] - moyennem) %*% t(c((Y[depart + (niterr-1)*batch + l,] - moyennem))) - V )/ norm(c(Y[depart + (niterr-1)*batch + l,] - moyennem) %*% t(c((Y[depart + (niterr-1)*batch + l,] - moyennem))) - V,type = "F")
    }
    gradV=gradV/(batch)
    #    if (gamma >= normm^(-1))
    #    {
    #      gamma = normm^(-1)
    #    }
    
    
    
    
    gradm=0
    #Mise à jour de mn
    for (l in 1:batch)
    {
      gradm=gradm+(Y[depart + (niterr-1)*batch + l,] - m)  / sqrt(sum((Y[depart + (niterr-1)*batch + l,] - m)^2))
    }
    
    m = m + gradm/(batch)* gamma
    
    
    
    
    
    #Mise à jour de V
    
    V = V + gamma*gradV
    
    sslog=sslog+log(niterr+1)^w
    #moyennem = moyennem*(i - depart)/(i - depart + 1) + 1/(i - depart +1)*m
    moyennem = moyennem + log(niterr+1)^w/sslog *(m - moyennem)
    
    
    # moyenneV <- (i - depart)/(i - depart  +1)*moyenneV + 1/(i- depart  + 1)*V
    moyenneV= moyenneV + log(niterr+1)^w/sslog *(V - moyenneV)
    for (l in 1:batch)
    {
      miter[depart + (niterr-1)*batch + l,] = moyennem
      VIter[,,depart + (niterr-1)*batch + l] = moyenneV
    }
    #U[i,,] <- orthonormalization(U[i,,],basis = TRUE, norm = TRUE)
    #Calcul des phijn
    #for (l in (1:q))
    #{
    # phijn[i,l] <- sum(t(Y[i+1,] - m)*Uphi[i,,l])
    #print(phijn[i,l])
    #}
    
    #vp = apply(U[i+1,,],2,norm)
    
    #Estimation des valeurs propres de la MCM
    #vpMCM[i,] = apply(U[i,,], 2, function(col) sqrt(sum(col^2)))
    
    #Rédiger partie simulation
    #vp  = lambdaEstim (précédente valeur de lambda) sample 100
    
    #Récupération des résultats de la fonction Outlier
    #resultatOutlier <- Outlier(donnee = Y[i, ],seuil_p_value = 0.05,VP = U[i,,],m = moyennem,lambdatilde)
    if(methode == "eigen"){
      #VPropresV <- eigen(moyenneV)$vectors
      elPropresV <- eigs_sym(moyenneV,k = ncol(Y))
      #VPropresV <- VPropresV %*% diag(1/sqrt(colSums(VPropresV^2)))
      VPropresV <- elPropresV$vectors
      valPV <- elPropresV$values 
      valPV = apply(cbind(valPV,rep(10^(-4),length(valPV))),MARGIN=1, FUN=max)
      lambdaResultat <- RobbinsMC2(c=cMC,mc_sample_size = niterRMon,w=w,vp=valPV,samp=1:niterRMon,init = lambdatilde,initbarre = lambda,ctilde = niterRMon*(niterr-1),cbarre =niterRMon*(niterr-1),slog=sum((log(1:((niterRMon*(niterr-1))+1))^w)))
        #lambdaResultat <- RobbinsMC2(c=cMC,mc_sample_size = niterRMon*(floor(log(batch*niterr))+1),w=w,vp=valPV,init = lambdatilde,initbarre = lambda,ctilde = compt,cbarre =compt,slog=sum((log(1:(compt+1))^w)))
      compt=compt+niterRMon*(floor(log(batch*niterr))+1)
      #ctilde = sampsize*(i-1),cbarre =sampsize*(i-1)
      lambda <- lambdaResultat$vp
      lambda = apply(cbind(lambda,rep(10^(-4),length(lambda))),MARGIN=1, FUN=max)
      lambdatilde <- lambdaResultat$lambda
      lambdatilde = apply(cbind(lambdatilde,rep(10^(-4),length(lambdatilde))),MARGIN=1, FUN=max)
      #print(lambda)
      #lambdaIter[i,] <- lambdatilde
      for (l in 1:batch){
        lambdaIter[depart + (niterr-1)*batch + l,] <- lambda
      }
      #Calcul de Sigma
      #Sigma <- VPropresV %*% diag(valPV) %*% t(VPropresV)
      #distances <- rep(0,nrow(Z))
      #distances <- calcule_vecteur_distances(Z,m,Sigma)
      #cutoff <- calcule_cutoff(distances,d)
      #resultatOutlier <- Outlier(donnee = Y[i, ],seuil_p_value = 0.05,VP = eigen(V)$vectors,m = moyennem,lambdatilde)
      VP <- VPropresV %*% diag(1/sqrt(colSums(VPropresV^2)))
      #print(lambdaIter[i,])
      varian= VP %*% diag(lambda) %*% t(VP) 
      for (l in 1:batch)
{        Sigma[depart + (niterr-1)*batch + l,,] <- varian 
      S <- 0
      for(j in (1:ncol(Y)))
     {
        S<- S + 1/(lambdaIter[depart + (niterr-1)*batch + l,j])*sum(t(Y[depart + (niterr-1)*batch + l,] - moyennem)*(VP[,j]))^2
      }
      #iterations <- depart + (niterr-1)*batch + l
      #print(l)
      distances[depart + (niterr-1)*batch + l] <- S
      if (distances[depart + (niterr-1)*batch + l] > cutoff) {outlier_labels[depart + (niterr-1)*batch + l] = 1}
      
      }
      }
    else {
      
      print("méthode sans eigen")
      
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
        
        #U[i+1,,] <- matrice
        #Orthogonalisation des vecteurs propres
        U[i+1,,] <- orthonormalization(U[i+1,,],basis = TRUE, norm = FALSE)
        
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
        lambdaResultat <- RobbinsMC2(c=cMC,mc_sample_size =  niterRMon,vp=vpMCM[i,],samp=1:sampsize,init = lambdatilde,initbarre = lambda,ctilde = sampsize*(i-1),cbarre = sampsize*(i-1),slog=sum((log(1:((sampsize*(i-1))+1))^w)))
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
      #distances[i] <- as.numeric(Y[i,] - m) %*% solve(Sigma[i,,]) %*% (as.numeric(t(Y[i,] - m)))
      #Calcul de la distance i
      S <- 0
      for(j in (1:ncol(Y)))
      {
        S <- S + 1/(lambdaIter[i,j])*(t(Y[i+1,] - moyennem)%*%VP[,j])^2
      }
      
      #print(solve(Sigma[i,,]))
      #distances[i] <- as.numeric(Y[i,] - m) %*% solve(Sigma[i,,]) %*% (as.numeric(t(Y[i,] - m)))
      distances[i] <- S
      
    }
    
    
    
    #matrice <- U[i+1,,]
    
    #Tri des vecteurs en ordre décroissant de norme
    #normes_colonnes <- apply(matrice, 2, function(col) sqrt(sum(col^2)))
    
    # order() pour obtenir les indices de tri par ordre décroissant des normes
    #indices_tri <- order(normes_colonnes, decreasing = TRUE)
    
    # Trier la matrice selon les normes décroissantes des colonnes
    #matrice <-matrice[,indices_tri]
    
  }
  
  
  miter[nrow(Y),] = miter[nrow(Y)-1,]
  VIter[,,nrow(Y)] = VIter[,,nrow(Y)-1]
  Sigma[nrow(Y),,] = Sigma[(nrow(Y) - 1),,]
  
  return (list(m=m,V=V,lambdatilde = lambdatilde,lambdaIter = lambdaIter,moyennem=moyennem,moyenneV=moyenneV,miter = miter,VIter = VIter,U = U,vpMCM = vpMCM,outlier_labels = outlier_labels,distances = distances, Sigma = Sigma,niter=niterr,VP=VP))
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

SampleCovOnline = function(Z)
{
  nblignes = nrow(Z)
  
  mean = 1.5*rnorm(ncol(Z))
  Sigma = matrix(0,ncol(Z),ncol(Z))
  
  meanIter =   matrix(0,nrow(Z),ncol(Z))
  SigmaIter = array(0, dim = c(nrow(Z),ncol(Z),ncol(Z)))
  
  for (i in (1:(nblignes-1)))
  {
    mean   = mean   + (1.0 / (i + 1)) * (Z[i+1,] - mean);
    Sigma = (i + 1)/i*Sigma + 1/(i + 1)*((Z[i+1,] - meanIter[i,])%*%t(Z[i+1,] - meanIter[i,]) - Sigma)
    meanIter[i,] = mean
    SigmaIter[i,,] = Sigma
    
    
    
  }
  SigmaIter[nrow(Z),,] = Sigma
  return(list(mean = mean, Sigma = Sigma, meanIter = meanIter, SigmaIter = SigmaIter))
}


#Shrinkage pour la médiane géométrique 

shrinkage_med <- function(Y,muVrai = mu1)
{
  
  #Calcul de la médiane composante par composante
  muCCm <- colMedians(Y)
  
  #Calcul de target
  
  nu <- (muCCm %*% rep(1,ncol(Y)))/ncol(Y)
  
  target <- rep(nu,ncol(Y))
  
  VarTot <- 0
  
  for (j in 1:ncol(Y))
  {
    VarTot <- VarTot + var(Y[,j])  
  }
  
  #Calcul de eta 
  
  eta <- pi/(2*ncol(Y))*VarTot
  
  eta <- eta/(eta + sum(diag(t(muVrai - target)%*%(muVrai - target))))
  
  muShrink <- eta*target + (1 - eta)*muCCm
  
  return (list(muCCm = muCCm, nu = nu,eta = eta,muShrink = muShrink))
}


#Reprise fonction package Michael Wolf changement de matrice de covariance empirique par SCCM

shrinkage_SCCM <- function(Y, k = -1) {
  dim.Y <- dim(Y)
  N <- dim.Y[1]  
  p <- dim.Y[2]  
  
  if (k < 0) {    
    Y <- scale(Y, scale = FALSE)  # Centrage des données
    k <- 1
  }
  n <- N - k    # Taille d'échantillon effective
  c <- p / n    # Ratio de concentration
  
  #  Calcul de la Matrice de Comédiane (COM)
  COM <- matrix(0, p, p)
  for (j in 1:p) {
    for (t in 1:p) {
      COM[j, t] <- median((Y[, j] - median(Y[,j])) * (Y[, t]  - median(Y[,t]))) 
    }
  }
  
  # Ajustement SCCM
  SCCM <- 2.198 * COM
  
  #  Calcul de la cible de shrinkage
  SCCM_diag <- diag(SCCM)
  sqrtvar <- sqrt(SCCM_diag)
  rBar <- (sum(SCCM / outer(sqrtvar, sqrtvar)) - p) / (p * (p - 1))
  #target <- rBar * outer(sqrtvar, sqrtvar)
  #diag(target) <- SCCM_diag  # Conservation de la diagonale
  target <- diag(ncol(SCCM))
  nu <- sum(diag(SCCM))/ncol(SCCM) 
  
  target <- nu*target
  
  # Estimation de pi (variance des éléments hors-diagonale)
  
  # SCCM2 <- 2.198 * COM2
  piMat <- SCCM - SCCM^2
  pihat <- sum(piMat)
  
  #  Estimation de gamma (distance entre SCCM et la cible)
  gammahat <- norm((SCCM - target), type = "F")^2
  
  # ⃣ Partie diagonale de rho 
  rho_diag <- sum(diag(piMat))
  
  #  Partie hors-diagonale de rho 
  thetaMat <- matrix(0, p, p)
  for (j in 1:p) {
    for (t in 1:p) {
      if (j != t) {
        thetaMat[j, t] <- median(Y[, j] * Y[, t]) - SCCM_diag[j] * SCCM_diag[t]
      }
    }
  }
  rho_off <- rBar * sum(outer(1/sqrtvar, sqrtvar) * thetaMat)
  
  #  Calcul de l'intensité du shrinkage
  rhohat <- rho_diag + rho_off
  kappahat <- (pihat - rhohat) / gammahat
  shrinkage <- max(0, min(1, kappahat / n))
  
  #  Calcul de l'estimateur shrinké
  sigmahat <- shrinkage * target + (1 - shrinkage) * SCCM
  
  return(list(SCCM_shrinked = sigmahat,SCCM = SCCM ,shrinkage_intensity = shrinkage))
}
