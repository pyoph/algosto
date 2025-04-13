#' @title Online robust estimation of the variance
#'
#' @description Given a data set X, this function estimates online and
#' robustly variance-covariance matrix via the geometric median and the Median Covariation Matrix.
#'
#' @param X Matrix of size N*p corresponding to the data. The rows are the observations.
#' @param c Stochastic gradient paramater (\eqn{\sqrt{ncol(X)}}) by default). Must be positive.
#' @param gamma Stochastic gradient parameter (0.75 by default). Must belong to (0,1).
#' @param batch Size of the batch. Default is {ncol(X)}.
#' @param w A parameter of the function RobbinsMC.
#' @param cMC A parameter of the function RobbinsMC.
#' @param gammaMC A parameter of the function RobbinsMC.
#' @param epsilon Stopping criterion: the algorithm stops when the difference between two iterations is less than epsilon, by default 1e-08
#' @param mc_sample_size Number of data generated in the Monte Carlo procedure to estimate the eigenvalues of the variance. Default is batch.
#' @param nitermax Maximum number of iterations for first estimation of the median and the Median Covariation Matrix with the help of Weiszfeld algorithm. Default is 100.
#' @param model Type of distribution: can be "Gaussian" (by default), "Student" or "Laplace"
#'
#' @return A list containing ??? TO DO Antoine
#'
#' @md
#'
#' @examples
#' N <- 1000
#' p <- 10
#' X <- matrix(rnorm(N*p),ncol=p)
#' StreamingRobustVariance(X=X)
#'
#' @export
#'
#'
StreamingRobustVariance <- function(X,c = sqrt(ncol(X)), gamma = 0.75,w=2, batch = ncol(X), cMC=ncol(X),gammaMC=0.75,
                                    Ninit = 100,mc_sample_size = batch,nitermax=100, model='Gaussian', df=3)
{
  checkmate::assertChoice(model, c("Gaussian", "Student", "Laplace"))
  checkmate::checkTRUE(nrow(X) > Ninit)
  vpMCM = matrix(0,ncol(X),ncol(X))
  lambdaInit =  rep(1,ncol(X))
  niterr=0
  lambdatilde = rep(1,ncol(X)) #initialization of the estimates of the eigenvalues of the Variance
  lambdaIter = matrix(0,nrow(X),ncol(X)) # intialization of a matrix containing all the estimates
  Sigma = array(0, dim = c(nrow(X),ncol(X),ncol(X))) # intialization of a matrix containing all the estimates of the variance
  U = array(0, dim = c(nrow(X), ncol(X),ncol(X)))  #Initialisation of an array contaiing all the estimates of the eigenvectors
  if (Ninit > 0)
  {
    minit=WeiszfeldMedian(X[1:Ninit,],nitermax = nitermax)$median
    init_cov = robustbase::covComed(X[1:Ninit,])$cov
    Vinit <- WeiszfeldMedianCovariance(X[1:Ninit,],minit,init_cov = init_cov,nitermax = 100)$median
    eig_init = eigen(Vinit,symmetric = TRUE)
    valPV <- eig_init$values
    valPV = apply(cbind(valPV,rep(10^(-4),length(valPV))),MARGIN=1, FUN=max)
    lambdaInit = valPV
    lambdatilde = valPV
    if(model == 'Gaussian'){
      UU = matrix(rnorm(ncol(X)*Ninit*10),ncol=ncol(X))
    }
    else if(model == 'Student'){
      UU <- matrix(rnorm(ncol(X)*Ninit*10)/sqrt(rchisq(1,df=df))*sqrt(df-2), ncol=ncol(X))
    }
    else{ # model == 'Laplace'
      UU <- LaplacesDemon::rmvl(Ninit*10,mu=rep(0,ncol(X)),Sigma=diag(ncol(X)))
    }
    lambdaResultat <- robbinsMC( U = UU ,  delta=valPV,  c = cMC,w=w , gamma = gammaMC)
    #ctilde = sampsize*(i-1),cbarre =sampsize*(i-1)
    lambda <- c(lambdaResultat$vp)
    lambdatilde <- c(lambdaResultat$vp)
    lambdaInit <- lambda
    VPropresV <- eig_init$vectors
    
    m <- minit  #initialisation of the estimates of the median
    V <- Vinit # initialisation of the estimates of the Median Covariation Matrix
    moyennem <- minit # initialisation of the averaged estimates
    moyenneV <- Vinit # initialisation of the averaged estimates
    VPropresV <- eig_init$vectors
    #    VP <- VPropresV %*% diag(1/sqrt(colSums(VPropresV^2)))
    VP <- normalize_columnsRcpp(VPropresV)   #print(lambdaIter[i,])
    #   varian= VP %*% diag(lambda) %*% t(VP)
    varianc = reconstruct_covarianceRcpp(VP,lambda)
    
    for (l in 1 :Ninit) ####stocking the estimates of the variance, eigenvectors and calculating the Mahalonib distance of the first data
    {
      U[l,,] <- eig_init$vectors
      Sigma[l,,] <- varianc
      lambdaIter[l,] = lambda
    }
  }
  
  #Stockage des itérations de matrices dans un tableau
  VIter <- array(0, dim = c(nrow(X),ncol(X), ncol(X)))
  slog <- 1 #To calculate the weights in RobbinsMC
  sslog = 1 # To calculate the weighs in RobbinsMC
  
  #Stocking the estimates of the medien
  miter = matrix(0,nrow(X),ncol(X))
  miter[1:Ninit,] = matrix(rep(m,Ninit),byrow=T,ncol=ncol(X))
  
  
  for (i  in (Ninit):(nrow(X)-1))
  {
    niterr=niterr+1
    if (niterr*batch+Ninit > nrow(X)){
      break
    }
    gamman = c*sqrt(batch)/(i)^(gamma)
    upd= update_median_covarianceRcpp(X,  m,  moyennem,  V,  moyenneV,  Ninit,  niterr,  batch,  gamman,  w,  sslog)
    m=upd$m
    moyennem = upd$moyennem
    V=upd$V
    moyenneV=upd$moyenneV
    for (l in 1:batch)
    {
      miter[Ninit + (niterr-1)*batch + l,] = moyennem
      VIter[Ninit + (niterr-1)*batch + l,,] = moyenneV
    }
    reseig=eigen(moyenneV,symmetric = TRUE)
    vpMCM=reseig$vectors
    VPropresV <- reseig$vectors
    valPV <- reseig$values
    valPV = apply(cbind(valPV,rep(10^(-4),length(valPV))),MARGIN=1, FUN=max)
    if(model == 'Gaussian'){
      UU = matrix(rnorm(mc_sample_size*ncol(X)),ncol=ncol(X))
    }
    else if(model == 'Student'){
      UU <- matrix(rnorm(mc_sample_size*ncol(X))/sqrt(rchisq(1,df=df))*sqrt(df-2), ncol=ncol(X))
    }
    else{ # model == 'Laplace'
      UU <- LaplacesDemon::rmvl(mc_sample_size,mu=rep(0,ncol(X)),Sigma=diag(ncol(X)))
    }
    lambdaResultat <-   robbinsMC(U=UU, c=cMC, gamma = gammaMC, w=w,delta=valPV,init = lambdatilde,
                                  init_bar = lambda,c_tilde = mc_sample_size*(niterr-1),
                                  c_bar =mc_sample_size*(niterr-1), sumlog=sum((log(1:((mc_sample_size*(niterr-1))+1))^w)))
    #ctilde = sampsize*(i-1),cbarre =sampsize*(i-1)
    lambda <- c(lambdaResultat$vp)
    lambda = apply(cbind(lambda,rep(10^(-4),length(lambda))),MARGIN=1, FUN=max)
    lambdatilde <- c(lambdaResultat$lambda)
    lambdatilde = apply(cbind(lambdatilde,rep(10^(-4),length(lambda))),MARGIN=1, FUN=max)
    lambdaIter[Ninit + (niterr-1)*batch + 1:batch,] <- matrix(rep(lambda,batch),byrow=T,nrow=batch)
    #    VP <- VPropresV %*% diag(1/sqrt(colSums(VPropresV^2)))
    VP <- normalize_columnsRcpp(VPropresV)   #print(lambdaIter[i,])
    #   varian= VP %*% diag(lambda) %*% t(VP)
    varian = reconstruct_covarianceRcpp(VP,lambda)
    for (l in 1:batch)
    {
      Sigma[ Ninit + (niterr-1)*batch + l,,] <- varian
    }
    
  }
  return (list(m=m,V=V,lambdatilde = lambdatilde,lambdaIter = lambdaIter,moyennem=moyennem,moyenneV=moyenneV,miter = miter,VIter = VIter,U = U,vpMCM = vpMCM,Sigma = Sigma,niter=niterr,VP=VP))
}




#' @title Online detection of outliers
#'
#' @description Given a data set X, this function detects in a streaming way
#' the outliers.
#'
#' @param X Matrix of size N*p corresponding to the data. The rows are the observations.
#' @param c Stochastic gradient paramater (\eqn{\sqrt{ncol(X)}}) by default). Must be positive.
#' @param gamma Stochastic gradient parameter (0.75 by default). Must belong to (0,1).
#' @param batch Size of the batch. Default is {ncol(X)}.
#' @param w A parameter of the function RobbinsMC.
#' @param cMC A parameter of the function RobbinsMC.
#' @param gammaMC A parameter of the function RobbinsMC.
#' @param cutoff Threshold at which a data item is considered contaminated. Defult is \eqn{qchisq(p = 0.95, df = ncol(X))}.
#' @param epsilon Stopping criterion: the algorithm stops when the difference between two iterations is less than epsilon, by default 1e-08
#' @param mc_sample_size Number of data generated in the Monte Carlo procedure to estimate the eigenvalues of the variance. Default is batch.
#' @param nitermax Maximum number of iterations for first estimation of the median and the Median Covariation Matrix with the help of Weiszfeld algorithm. Default is 100.
#'
#' @return A list containing ???
#'
#' @md
#'
#' @examples
#' N <- 1000
#' p <- 10
#' X <- matrix(rnorm(N*p),ncol=p)
#' StreamingOutlierDetection(X=X)
#'
#' @export
#'
#'
StreamingOutlierDetection <- function(X,c = sqrt(ncol(X)), gamma = 0.75,w=2, batch = ncol(X), cMC=ncol(X),gammaMC=0.75,
                                      Ninit = 100,cutoff =qchisq(p = 0.95, df = ncol(X)),mc_sample_size = batch,
                                      nitermax=100)
{
  checkmate::checkTRUE(nrow(X) > Ninit)
  lambdaInit =  rep(1,ncol(X))
  niterr=0
  lambdatilde = rep(1,ncol(X)) #initialization of the estimates of the eigenvalues of the Variance
  lambdaIter = matrix(0,nrow(X),ncol(X)) # intialization of a matrix containing all the estimates
  Sigma = array(0, dim = c(nrow(X),ncol(X),ncol(X))) # intialization of a matrix containing all the estimates of the variance
  distances <- rep(0,nrow(X)) #initialization of a vector containing all the Mahalanobis distances
  outlier_labels = rep(0,nrow(X)) # initialization of a vector containing all the outlier labels (1 if outlier, 0 else).
  U = array(0, dim = c(nrow(X), ncol(X),ncol(X)))  #Initialisation of an array contaiing all the estimates of the eigenvectors
  if (Ninit > 0)
  {
    minit=WeiszfeldMedian(X[1:Ninit,],nitermax = nitermax)$median
    init_cov = robustbase::covComed(X[1:Ninit,])$cov
    Vinit <- WeiszfeldMedianCovariance(X[1:Ninit,],minit,init_cov = init_cov,nitermax = 100)$median
    eig_init = eigen(Vinit,symmetric = TRUE)
    valPV <- eig_init$values
    valPV = apply(cbind(valPV,rep(10^(-4),length(valPV))),MARGIN=1, FUN=max)
    lambdaInit = valPV
    lambdatilde = valPV
    UU=matrix(rnorm(ncol(X)*Ninit*10),ncol=ncol(X))
    lambdaResultat <- robbinsMC( U = UU ,  delta=valPV,  c = cMC,gamma=gammaMC, w=w)
    #ctilde = sampsize*(i-1),cbarre =sampsize*(i-1)
    lambda <- c(lambdaResultat$vp)
    lambdatilde <- c(lambdaResultat$vp)
    lambdaInit <- lambda
    VPropresV <- eig_init$vectors
    vpMCM=eig_init$vectors
    m <- minit  #initialisation of the estimates of the median
    V <- Vinit # initialisation of the estimates of the Median Covariation Matrix
    moyennem <- minit # initialisation of the averaged estimates
    moyenneV <- Vinit # initialisation of the averaged estimates
    VPropresV <- eig_init$vectors
    #    VP <- VPropresV %*% diag(1/sqrt(colSums(VPropresV^2)))
    VP <- normalize_columnsRcpp(VPropresV)   #print(lambdaIter[i,])
    #   varian= VP %*% diag(lambda) %*% t(VP)
    varianc = reconstruct_covarianceRcpp(VP,lambda)
    
    for (l in 1 :Ninit) ####stocking the estimates of the variance, eigenvectors and calculating the Mahalonib distance of the first data
    {
      U[l,,] <- eig_init$vectors
      Sigma[l,,] <- varianc
      lambdaIter[l,] = lambda
      S <-  mahalanobis_generalizedRcpp(X[l,],  moyennem, vpMCM,  lambda)
      distances[l] <- S
      if (S > cutoff) {outlier_labels[l] = 1}
      
      
    }
  }
  
  #Stockage des itérations de matrices dans un tableau
  VIter <- array(0, dim = c(nrow(X),ncol(X), ncol(X)))
  slog <- 1 #To calculate the weights in RobbinsMC
  sslog = 1 # To calculate the weighs in RobbinsMC
  
  #Stocking the estimates of the medien
  miter = matrix(0,nrow(X),ncol(X))
  miter[1:Ninit,] = matrix(rep(m,Ninit),byrow=T,ncol=ncol(X))
  
  for (i  in (Ninit):(nrow(X)-1))
  {
    niterr=niterr+1
    if (niterr*batch+Ninit > nrow(X)){
      break
    }
    gamman = c*sqrt(batch)/(i)^(gamma)
    upd= update_median_covarianceRcpp(X,  m,  moyennem,  V,  moyenneV,  Ninit,  niterr,  batch,  gamman,  w,  sslog)
    m=upd$m
    moyennem = upd$moyennem
    V=upd$V
    moyenneV=upd$moyenneV
    for (l in 1:batch)
    {
      miter[Ninit + (niterr-1)*batch + l,] = moyennem
      VIter[Ninit + (niterr-1)*batch + l,,] = moyenneV
    }
    reseig=eigen(moyenneV,symmetric = TRUE)
    vpMCM=reseig$vectors
    VPropresV <- reseig$vectors
    valPV <- reseig$values
    valPV = apply(cbind(valPV,rep(10^(-4),length(valPV))),MARGIN=1, FUN=max)
    UU=matrix(rnorm(mc_sample_size*ncol(X)),ncol=ncol(X))
    lambdaResultat <-   robbinsMC(U=UU, c=cMC,gamma=gammaMC, w=w,delta=valPV,init = lambdatilde,
                                  init_bar = lambda,c_tilde = mc_sample_size*(niterr-1),
                                  c_bar =mc_sample_size*(niterr-1), sumlog=sum((log(1:((mc_sample_size*(niterr-1))+1))^w)))
    #ctilde = sampsize*(i-1),cbarre =sampsize*(i-1)
    lambda <- c(lambdaResultat$vp)
    lambda = apply(cbind(lambda,rep(10^(-4),length(lambda))),MARGIN=1, FUN=max)
    lambdatilde <- c(lambdaResultat$lambda)
    lambdatilde = apply(cbind(lambdatilde,rep(10^(-4),length(lambda))),MARGIN=1, FUN=max)
    lambdaIter[Ninit + (niterr-1)*batch + 1:batch,] <- matrix(rep(lambda,batch),byrow=T,nrow=batch)
    #    VP <- VPropresV %*% diag(1/sqrt(colSums(VPropresV^2)))
    VP <- normalize_columnsRcpp(VPropresV)   #print(lambdaIter[i,])
    #   varian= VP %*% diag(lambda) %*% t(VP)
    varian = reconstruct_covarianceRcpp(VP,lambda)
    for (l in 1:batch)
    {        Sigma[ Ninit + (niterr-1)*batch + l,,] <- varian
    S <-  mahalanobis_generalizedRcpp(X[Ninit + (niterr-1)*batch + l,],  moyennem, vpMCM,  lambda)
    #iterations <- Ninit + (niterr-1)*batch + l
    #print(l)
    distances[Ninit + (niterr-1)*batch + l] <- S
    if (distances[Ninit + (niterr-1)*batch + l] > cutoff) {outlier_labels[Ninit + (niterr-1)*batch + l] = 1}
    
    }
    
  }
  return (list(m=m,V=V,lambdatilde = lambdatilde,lambdaIter = lambdaIter,moyennem=moyennem,moyenneV=moyenneV,miter = miter,VIter = VIter,U = U,vpMCM = vpMCM,outlier_labels = outlier_labels,distances = distances, Sigma = Sigma,niter=niterr,VP=VP))
}

#' @title Offline robust estimation of the variance
#'
#'
#' @description Given a data set X, this function estimates
#' robustly variance-covariance matrix via the geometric median and the Median Covariation Matrix.
#'
#' @param X Matrix of size N*p corresponding to the data. The rows are the observations.
#' @param methodMC Method to estimate the eigenvalues of the variance: can be "FixMC", "GradMC" and "RobbinsMC" (by default).
#' @param methodMCM Method for estimating the Median Covariation Matrix: can be "Weiszfeld" (by default) or "ASG".
#' @param model Type of distribution: can be "Gaussian" (by default), "Student" or "Laplace"
#' @param init A row vector for initializing the estimates of the median. Default is rep(0,ncol(X)).
#' @param epsMCM Stopping criterion for Weiszfeld algorithm: the algorithm stops when the difference between two iterations is less than epsilon, by default 1e-08
#' @param niterWeisz The maximum number of iteration for the Weiszfeld algorithm. Default is 50.
#' @param gammaMCM A positive constant if methodMCM='ASG'.
#' @param alphaMCM A positive constant if methodMCM='ASG'. Must belong to (0,1).
#' @param nstart Number of run if methodMCM='ASG'.
#' @param mc_sample_size Number of data generated in the Monte Carlo procedure to estimate the eigenvalues of the variance. Default is batch.
#' @param df Number of degress of freedom if %pdem='Student'.
#' @param epsilon Stopping criterion for MethodMC algorithm: the algorithm stops when the difference between two iterations is less than epsilon, by default 1e-08
#' @param niterMC Maximum number of iterations is MethodMC='FixMC' or 'GradMC'.
#' @param cMC A parameter of the function RobbinsMC.
#' @param gammaMC A parameter of the function RobbinsMC.
#' @param w A parameter of the function RobbinsMC.
#'
#' @return A list containing ??? TO DO Antoine
#'
#' @md
#'
#' @examples
#' N <- 1000
#' p <- 10
#' X <- matrix(rnorm(N*p),ncol=p)
#' OfflineRobustVariance(X=X)
#'
#' @export
#'
OfflineRobustVariance <- function(X,methodMC='RobbinsMC', methodMCM='Weiszfeld', model='Gaussian',
                                  init=rep(0,ncol(X)),
                                  init_cov=diag(ncol(X)),
                                  epsMCM=1e-08,niterWeisz=50,
                                  gammaMCM=2, alphaMCM=0.75, nstart=1,
                                  mc_sample_size=1000,
                                  df=3, epsilon=1e-08,
                                  niterMC=50,
                                  cMC=2, gammaMC=0.75, w=2)
{
  checkmate::assertChoice(methodMC, c("RobbinsMC", "GradMC", "FixMC"))
  checkmate::assertChoice(methodMCM, c("Weiszfeld",  "ASG"))
  checkmate::assertChoice(model, c("Gaussian", "Student", "Laplace"))
  d <- ncol(X)
  if(methodMCM=='Weiszfeld'){
    median <- WeiszfeldMedian(X=X, init=init, epsilon=epsMCM, nitermax=niterWeisz)
    medianCov <- WeiszfeldMedianCovariance(X=X, median_est=median$median, init_cov=init_cov, epsilon=epsMCM, nitermax=niterWeisz)
  }
  else if(methodMCM=='ASG')
  {
    median <- ASGMedian(X=X, init=init, gamma=gammaMCM, alpha=alphaMCM, nstart=nstart, epsilon=epsMCM)
    medianCov <- ASGMedianCovariance(X=X, median_est=median$median, init_cov=init_cov, gamma=gammaMCM, alpha=alphaMCM, nstart=nstart)
  }
  eig <- eigen(medianCov$median,symmetric = TRUE)
  eigenvec <- eig$vectors
  eigenval <- eig$values
  
  # Generate U for eigen values computation
  p <- length(eigenval)
  
  if(model == 'Gaussian'){
    U = matrix(rnorm(mc_sample_size*p),ncol=p)
  }
  else if(model == 'Student'){
    U <- matrix(rnorm(mc_sample_size*p)/sqrt(rchisq(1,df=df))*sqrt(df-2), ncol=p)
  }
  else{ # model == 'Laplace'
    U <- LaplacesDemon::rmvl(mc_sample_size,mu=rep(0,p),Sigma=diag(p))
  }
  
  # Compute eigen values of variance-covariance matrix using eigen values of median-covariance matrix
  if (methodMC=="RobbinsMC"){
    eigen_est <- robbinsMC(U=U,delta=eigenval,gamma=gammaMC,c=cMC,w=w,epsilon=epsilon)
  }
  else if (methodMC=="GradMC"){
    eigen_est <- gradMC(U=U,delta=eigenval,niter=niterMC,epsilon=epsilon,step=cMC*(1:niterMC)^(-0.5))
  }
  else if (methodMC=="FixMC"){
    eigen_est <- fixMC(U=U,delta=eigenval,niter=niterMC,epsilon=epsilon)
  }
  
  lambda <- c(eigen_est$vp)
  variance <- t(matrix(eigenvec,ncol=d,byrow=TRUE))%*%diag(lambda)%*%(matrix(eigenvec,ncol=d,byrow=TRUE))
  resultat <- list(median=median$median,variance=variance,covmedian=medianCov$median)
  return(resultat)
}


#' @title Offline detection of outliers
#'
#'
#' @description Given a data set X, this function estimates
#' robustly variance-covariance matrix via the geometric median and the Median Covariation Matrix.
#'
#' @param X Matrix of size N*p corresponding to the data. The rows are the observations.
#' @param methodMC Method to estimate the eigenvalues of the variance: can be "FixMC", "GradMC" and "RobbinsMC" (by default).
#' @param methodMCM Method for estimating the Median Covariation Matrix: can be "Weiszfeld" (by default) or "ASG".
#' @param cutoff Threshold at which a data item is considered contaminated. Defult is \eqn{qchisq(p = 0.95, df = ncol(X))}.
#' @param model Type of distribution: can be "Gaussian" (by default), "Student" or "Laplace"
#' @param init A row vector for initializing the estimates of the median. Default is rep(0,ncol(X)).
#' @param epsMCM Stopping criterion for Weiszfeld algorithm: the algorithm stops when the difference between two iterations is less than epsilon, by default 1e-08
#' @param niterWeisz The maximum number of iteration for the Weiszfeld algorithm. Default is 50.
#' @param gammaMCM A positive constant if methodMCM='ASG'.
#' @param alphaMCM A positive constant if methodMCM='ASG'. Must belong to (0,1).
#' @param nstart Number of run if methodMCM='ASG'.
#' @param mc_sample_size Number of data generated in the Monte Carlo procedure to estimate the eigenvalues of the variance. Default is batch.
#' @param df Number of degress of freedom if %pdem='Student'.
#' @param epsilon Stopping criterion for MethodMC algorithm: the algorithm stops when the difference between two iterations is less than epsilon, by default 1e-08
#' @param niterMC Maximum number of iterations is MethodMC='FixMC' or 'GradMC'.
#' @param cMC A parameter of the function RobbinsMC.
#' @param gammaMC A parameter of the function RobbinsMC.
#' @param w A parameter of the function RobbinsMC.
#'
#' @return A list containing ??? TO DO Antoine
#'
#' @md
#'
#' @examples
#' N <- 1000
#' p <- 10
#' X <- matrix(rnorm(N*p),ncol=p)
#' OfflineOutlierDetection(X=X)
#'
#' @export
#'
OfflineOutlierDetection <- function(X,methodMC='RobbinsMC', methodMCM='Weiszfeld',cutoff =qchisq(p = 0.95, df = ncol(X)),
                                    init=rep(0,ncol(X)),
                                    init_cov=diag(ncol(X)),
                                    epsMCM=1e-08,niterWeisz=50,
                                    gammaMCM=2, alphaMCM=0.75, nstart=1,
                                    mc_sample_size=1000,
                                    df=3, epsilon=1e-08,
                                    niterMC=50,
                                    cMC=2, gammaMC=0.75, w=2)
{
  checkmate::assertChoice(methodMC, c("RobbinsMC", "GradMC", "FixMC"))
  checkmate::assertChoice(methodMCM, c("Weiszfeld",  "ASG"))
  d <- ncol(X)
  distances=rep(0,nrow(X))
  outlier_labels=rep(0,nrow(X))
  if(methodMCM=='Weiszfeld'){
    median <- WeiszfeldMedian(X=X, init=init, epsilon=epsMCM, nitermax=niterWeisz)
    medianCov <- WeiszfeldMedianCovariance(X=X, median_est=median$median, init_cov=init_cov, epsilon=epsMCM, nitermax=niterWeisz)
  }
  else if(methodMCM=='ASG')
  {
    median <- ASGMedian(X=X, init=init, gamma=gammaMCM, alpha=alphaMCM, nstart=nstart, epsilon=epsMCM)
    medianCov <- ASGMedianCovariance(X=X, median_est=median$median, init_cov=init_cov, gamma=gammaMCM, alpha=alphaMCM, nstart=nstart)
  }
  eig <- eigen(medianCov$median,symmetric = TRUE)
  eigenvec <- eig$vectors
  eigenval <- eig$values
  
  # Generate U for eigen values computation
  p <- length(eigenval)
  U = matrix(rnorm(mc_sample_size*p),ncol=p)
  # Compute eigen values of variance-covariance matrix using eigen values of median-covariance matrix
  if (methodMC=="RobbinsMC"){
    eigen_est <- robbinsMC(U=U,delta=eigenval,gamma=gammaMC,c=cMC,w=w,epsilon=epsilon)
  }
  else if (methodMC=="GradMC"){
    eigen_est <- gradMC(U=U,delta=eigenval,niter=niterMC,epsilon=epsilon,step=cMC*(1:niterMC)^(-0.5))
  }
  else if (methodMC=="FixMC"){
    eigen_est <- fixMC(U=U,delta=eigenval,niter=niterMC,epsilon=epsilon)
  }
  lambda <- c(eigen_est$vp)
  variance <- t(matrix(eigenvec,ncol=d,byrow=TRUE))%*%diag(lambda)%*%(matrix(eigenvec,ncol=d,byrow=TRUE))
  #### outlier detection
  #Calcul des distances
  for (i in 1:nrow(X))
    
  {
    S <-  mahalanobis_generalizedRcpp(X[i,], median$median, eigenvec,  lambda)
    distances[i] <- S
    if (distances[i] > cutoff) {outlier_labels[i] = 1}
    
  }
  return(list(median=median$median,variance=variance,covmedian=medianCov$median,outlier_labels=outlier_labels,distances=distances))
  
}












