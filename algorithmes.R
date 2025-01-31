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




estimMV <- function(Y,c = sqrt(ncol(Y)), exposantPas = 0.75,aa = 1,r = 1.5, 
                    minit = r*rnorm(ncol(Y)),Vinit = diag(ncol(Y))
                    ,U = array(1, dim = c(nrow(Y), ncol(Y),ncol(Y))),
                    vpMCM = matrix(0,ncol(Y),ncol(Y)),lambdaInit =  rep(1,ncol(Y))
                    ,SigmaInit = diag(d),methode = "eigen",depart = 0,cutoff =qchisq(p = 0.95, df = ncol(Y)),niterRMon = ncol(Y))
{
  print(depart)
  
  lambdatilde = rep(1,ncol(Y))
  lambdaIter = matrix(0,nrow(Y),ncol(Y))
  U = array(0, dim = c(nrow(Y), ncol(Y),ncol(Y)))
  Sigma = array(0, dim = c(nrow(Y), ncol(Y),ncol(Y)))
  #Stockage des distances
  
  distances <- rep(0,nrow(Y)) 
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
    print("depart > 0")
    resoff=RobVar(Y[1:100,],mc_sample_size = nrow(Y),c=ncol(Y),w=2)
    eig_init=eigen(resoff$variance)
    lambdaInit=eig_init$values
    lambdatilde=lambdaInit
    lambdaIter[1:depart,]=matrix(rep(lambdaInit,depart),byrow=T,nrow=depart)
    moyennem <- resoff$median
    minit=resoff$median
    Vinit=resoff$covmedian
    
    for (l in 1 :depart)
    {
      U[l,,] <- eig_init$vectors
      Sigma[l,,] <- resoff$variance
      #S <- 0
      #for(j in (1:ncol(Y)))
      #{
       # S <- S + 1/(lambdaInit[j])*(sum((Y[l,] - minit)*(eig_init$vectors)[,j]))^2
        #}
      
      
      #print(solve(Sigma[i,,]))
      distances[l] <- as.numeric(Y[l,] - minit) %*% solve(Sigma[l,,]) %*% (as.numeric(t(Y[l,] - minit)))
      #distances[l] <- S
      #distances[l] <- calcule_vecteur_distances(Z,minit,Sigma[l,,])
      print(distances[l])
      
      
      }
  
  }
  
  
  sampsize = niterRMon
  lambda = lambdaInit
  lambdatilde = lambdatilde
  
  
  #Stockage des estimations des valeurs propres de la matrice de covariance
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
  
  
  #Stockage des itérations
  miter = matrix(0,nrow(Y),ncol(Y))
  
  # Vecteur pour stocker les labels des outliers
  outlier_labels <- rep(0, nrow(Y)-1)
  
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
    
    moyenneV <- (i - depart)/(i - depart  +1)*moyenneV + 1/(i -depart + 1)*V
    
    
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
      VPropresV <- eigen(moyenneV)$vectors
      #VPropresV <- VPropresV %*% diag(1/sqrt(colSums(VPropresV^2)))
      
      valPV <- eigen(V)$values  
      lambdaResultat <- RobbinsMC2(niterRMon,vp=valPV,samp=1:sampsize,init = lambdatilde,initbarre = lambda,ctilde = sampsize*(i-1),cbarre =sampsize*(i-1),slog=sum((log(1:((sampsize*(i-1))+1))^w)))
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
      #S <- 0
      
       # for(j in (1:ncol(Y)))
          
        #{
          
         # S<- S + 1/(lambdaIter[i,j])*sum(t(Y[i+1,] - moyennem)*(VP[,j]))^2
          
        #}
      #S <- 0
      #for(j in (1:ncol(Y)))
      #{
       # S <- S + 1/(lambdaIter[i,j])*(sum((Y[i,] - moyennem)*(eig_init$vectors)[,j]))^2
      #}
        
      
      
      #print(solve(Sigma[i,,]))
      distances[i] <- as.numeric(Y[i,] - m) %*% solve(Sigma[i,,]) %*% (as.numeric(t(Y[i,] - m)))
      #distances[i] <- S
      
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
        lambdaResultat <- RobbinsMC2(niterRMon,vp=vpMCM[i,],samp=1:sampsize,init = lambdatilde,initbarre = lambda,ctilde = sampsize*(i-1),cbarre = sampsize*(i-1),slog=sum((log(1:((sampsize*(i-1))+1))^w)))
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
  
  
  
  return (list(m=m,V=V,lambdatilde = lambdatilde,lambdaIter = lambdaIter,moyennem=moyennem,moyenneV=moyenneV,miter = miter,VIter = VIter,U = U,vpMCM = vpMCM,outlier_labels = outlier_labels,distances = distances, Sigma = Sigma))
}


#Fonction qui prend en paramètre une matrice de données, et une méthode d'estimation et renvoie les paramètres estimés, les distances et les outliers
detection <- function(Y,c = ncol(Y), exposantPas = 0.75,aa = 1,r = 1.5,sampsize = ncol(Y), 
                      minit = r*rnorm(ncol(Y)),Vinit = diag(ncol(Y))
                      ,U = array(1, dim = c(nrow(Y), ncol(Y),ncol(Y))),
                      vpMCM = matrix(0,n,ncol(Y)),methodeOnline = "eigen", methodeEstimation = "offline",depart_online = 100,niterRMon = 100)
{
  distances <- rep(0, nrow(Y))
  if (methodeEstimation == "offline")
  {
   
    results <- RobVar(Y,mc_sample_size = nrow(Z),c=ncol(Z),w=2)
    
    #Récupération de la médiane de Sigma des vecteurs propres (variable U), et des valeurs propres
    med <- results$median
    SigmaOffline <- results$variance
    V <- results$covmedian
    U <- eigen(V)$vectors
    #lambda <- RobbinsMC()
    distances <- calcule_vecteur_distances(Y, med, SigmaOffline)
  }
  else if (methodeEstimation == "online")
  { 
    results <- estimMV(Y, depart = depart_online)
    #Retour des résultats
    med <- results$moyennem
    SigmaOnline <- results$Sigma
    U <- results$U
    lambda <- results$lambda
    V <- results$moyenneV
    distances <- results$distances
   # results$Sigma[101,,]
  }
  if (methodeEstimation  == "online") {
    return(list(med = med, SigmaOnline = SigmaOnline,U = U ,lambda = lambda,V = V, distances = distances))
}
  else {return(list(med = med, SigmaOffline = SigmaOffline ,V = V, distances = distances))}
}


#Fonction streaming en paramètre matrice de données renvoie les paramètres estimés et des distances
streaming <- function(Y,k = 1,c = sqrt(ncol(Y)),r = 1.5,exposantPas = 0.75,aa = 1, 
                      minit = r*rnorm(ncol(Y)),Vinit = diag(ncol(Y))
                      ,U = array(1, dim = c(nrow(Y), ncol(Y),ncol(Y))),
                      vpMCM = matrix(0,ncol(Y),ncol(Y)),lambdaInit =  rep(1,ncol(Y))
                      ,SigmaInit = diag(d),methode = "eigen",depart = 100/k,cutoff =qchisq(p = 0.95, df = ncol(Y)),niterRMon = ncol(Y))
{

  sampsize = ncol(Y)
  
  t <- nrow(Y)/k
 
  print(t)  
  
  Sigma = array(0, dim = c(t, ncol(Y),ncol(Y)))
  

  #Stockage des estimations des valeurs propres de la matrice de covariance
  lambdaIter = matrix(0,t,ncol(Y))
  #lambda <- eigen(Vvrai)$values
  lambda <- rep(1,ncol(Y))
  lambdatilde <- rep(1,ncol(Y))

  #vp2 <- eigen(Vvrai)$values
  slog <- 1
  w <- 2

  beta = 1/2



  #Stockage des itérations de médiane géométrique

  miter = matrix(0,t,ncol(Y))

  stat <- rep(0,(t-1))
  phat <- rep(0,(t-1))
  #Stockage des itérations de matrices dans un tableau
  VIter <- array(0, dim = c(ncol(Y), ncol(Y), t))

  #Stockage des vecteurs propres

  U <- array(1, dim = c(t, ncol(Y), ncol(Y)))


  #Valeurs propres de V
  vpMCM <- matrix(0,t,ncol(Y))


  #Initialisation de U avec des vecteurs propres pas trop éloignés de ceux de cov(X)

  #U[1,,] <-  1.5*diag(ncol(Y))
  #Si départ = 0, 
  #intialisation de U et de m sur la sphère unité  
  if (depart == 0)
  {
    
    matrix_random <- matrix(0, ncol(Y), ncol(Y))
    
    #Générer chaque colonne aléatoire sur la sphère unité
    for (i in (1:ncol(Y))) {
      v <- rnorm(ncol(Y))  # Tirer un vecteur d composantes normales
      matrix_random[, i] <- v / sqrt(sum(v^2))  # Normaliser
    }
    
    
    U[1,,] <- matrix_random
    #Initialisation de m
    
    m <- r*rnorm(ncol(Y))
    
    V = diag(ncol(Y))
    
    moyennem = m
    
    
    
    
    moyenneV <- V
    
  }
  if (depart > 0)
  {
    print("depart > 0")
    resoff=RobVar(Y[1:depart/k,],mc_sample_size = nrow(Y),c=ncol(Y),w=2)
    eig_init=eigen(resoff$variance)
    lambdaInit=eig_init$values
    lambdatilde=lambdaInit
    print("test avec lambdaIter")
    lambdaIter[1:(depart/k),]=matrix(rep(lambdaInit,(depart/k)),byrow=T,nrow=(depart/k))
    minit=resoff$median
    moyennem <- t(minit)
    m <- t(minit)
    Vinit=resoff$covmedian
    V <- Vinit
  
    moyennem = m
  
    moyenneV <- V
    
    
    for (l in 1 :(depart/k))
    {
      U[l,,] <- eig_init$vectors
      Sigma[l,,] <- resoff$variance
      #S <- 0
      #for(j in (1:ncol(Y)))
      #{
      # S <- S + 1/(lambdaInit[j])*(sum((Y[l,] - minit)*(eig_init$vectors)[,j]))^2
      #}
      
      
      #print(solve(Sigma[i,,]))
      distances[l] <- as.numeric(Y[l,] - minit) %*% solve(Sigma[l,,]) %*% (as.numeric(t(Y[l,] - minit)))
      #distances[l] <- S
      #distances[l] <- calcule_vecteur_distances(Z,minit,Sigma[l,,])
      print(distances[l])
      
      
    }
    
  }
  
  
  # Vecteur pour stocker les labels des outliers
  outlier_labels <- rep(0, (t-1))

  phat <- rep(0,t-1)


  finboucle = t - 1

  for (i in ((1+depart):finboucle))
  {

    #print(i)

    #Somme des médianes géométriques
    S <- rep(0,ncol(Y))

    #Somme matrices

    Smatr <- matrix(0,ncol(Y),ncol(Y))

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


    moyennem = moyennem*(i-depart)/(i - depart +1) + 1/(i - depart +1)*m

    miter[i,] = moyennem





    for (j in 1 : k)
    {
      Smatr <- Smatr +  ((Y[(i - 1) * k + j, ] - moyennem) %*% t(Y[(i - 1) * k + j, ] - m)-V)/norm((Y[(i - 1) * k + j, ] - moyennem) %*% t(Y[(i - 1) * k + j, ] - m)-V,type = "F")

    }


    V <- V + gamma/k*k^(beta)*Smatr







    #Mise à jour de moyenneV

    moyenneV <- (i - depart)/(i - depart +1)*moyenneV + 1/(i - depart + 1)*V


    #Stockage des itérations de matrices dans un tableau
    VIter[,,i] <- moyenneV


    #Estimation des vecteurs propres de Vt et orthonormalisation
if(methode == "eigen"){
  VPropresV <- eigen(moyenneV)$vectors
  #VPropresV <- VPropresV %*% diag(1/sqrt(colSums(VPropresV^2)))
  valPV <- eigen(moyenneV)$values  
  #print(VPropresV)
  
  #lambdaResultat <- RobbinsMC2(niterRMon,vp=valPV,samp=1:sampsize,init = lambdatilde,initbarre = lambda,ctilde = sampsize*(i-1),cbarre =sampsize*(i-1),slog=sum((log(1:((sampsize*(i-1))+1))^w)))
  lambdaResultat <- RobbinsMC2(niterRMon,vp=valPV,samp=1:sampsize,init = lambdatilde,initbarre = lambda,ctilde = sampsize*(i-1),cbarre =sampsize*(i-1))
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
  #S <- 0
  
  # for(j in (1:ncol(Y)))
  
  #{
  
  # S<- S + 1/(lambdaIter[i,j])*sum(t(Y[i+1,] - moyennem)*(VP[,j]))^2
  
  #}
  #S <- 0
  #for(j in (1:ncol(Y)))
  #{
  # S <- S + 1/(lambdaIter[i,j])*(sum((Y[i,] - moyennem)*(eig_init$vectors)[,j]))^2
  #}
  
  
  
  #print(solve(Sigma[i,,]))
  distances[i] <- as.numeric(Y[i,] - m) %*% solve(Sigma[i,,]) %*% (as.numeric(t(Y[i,] - m)))
  #distances[i] <- S
  
}
    
    else{
    for (l in (1:ncol(Y)))
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
    lambdaResultat <- RobbinsMC2(niterRMon,vp=vpMCM[i,],samp=1:sampsize,init = lambdatilde,initbarre = lambda,ctilde = sampsize*(i-1),cbarre = sampsize*(i-1),slog=sum((log(1:((sampsize*(i-1))+1))^w)))
    lambda <- lambdaResultat$vp
    lambdatilde <- lambdaResultat$lambda
    #print(lambda)
    lambdaIter[i,] <- lambdatilde
    
    
    
      }
    #Rédiger partie simulation
    #vp  = lambdaEstim (précédente valeur de lambda) sample 100

    #Récupération des résultats de la fonction Outlier
    #resultatOutlier <- Outlier(donnee = Y[i, ],seuil_p_value = 0.05,VP = U[i,,],m = moyennem,lambdatilde)
    #resultatOutlier <- Outlier(donnee = Y[i, ],seuil_p_value = 0.05,VP = eigen(V)$vectors,m = moyennem,lambdatilde)

    # Vérifier si l'entrée est un outlier
    #if (resultatOutlier$outlier_label) {

      #Estimation des valeurs propres selon une procédure de Robbins-Monro

      #outlier_labels[i] <- 1


#    }
 #   stat[i] <- resultatOutlier$S

  #  phat[i] <- resultatOutlier$phat
    #print(phat[i])


  }

  
  miter[t,] = miter[(t-1),]
  VIter[,,t] = VIter[,,t-1]
  Sigma[t,,] = Sigma[(t - 1),,]
  
  
  return (list(m=m,V=V,lambdatilde = lambdatilde,lambdaIter = lambdaIter,moyennem=moyennem,moyenneV=moyenneV,miter = miter,VIter = VIter,U = U,vpMCM = vpMCM,Sigma = Sigma))
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

