#' @title Online robust estimation of the variance
#'
#' @description Given a data set X, this function estimates online and
#' robustly variance-covariance matrix via the geometric median and the Median Covariation Matrix.
#'
#' @param X Matrix of size n*d corresponding to the data. The rows are the observations.
#' @param nDataInit Number of data used for initializing the algorithm. Default is 100.
#' @param model Type of distribution: can be "Gaussian" (by default), "Student" or "Laplace"
#' @param df Number of degrees of freedom if model='Student'.
#' @param computeOutliers Flag for outliers computation, FALSE by default. If TRUE, only Gaussian model is possible.
#' @param cutoff Threshold at which a data item is considered contaminated (if computeOutliers is TRUE). Default is 0.95.
#' @param cutoff_method The cutoff method used to flag as outliers. Can be "Chi-square" if the distances are considered to follow a Chi-square distribution, or "quantile" if the estimated quantile of order "cutoff" of the distance is used.
#' @param cutinit A parameter giving the order of the quantile used for initializing the estimate of quantile of the distance if "quantile" method is used. Default is 0.6.
#' @param c_m The constant in the stepsequence of Rbobbins-Monro algorithm for estimating the quantile of the distances. Default is 2.
#' @param batch Size of the batch. Default is ncol(X).
#' @param niterWeisz Maximum number of iterations for first estimation of the median and the Median Covariation Matrix with the help of Weiszfeld algorithm. Default is 100.
#' @param sizeMC Number of data generated in the Monte Carlo procedure to estimate the eigenvalues of the variance. Default is batch.
#' @param epsilon Stopping criterion for RobbinsMC and Weisfeld algorithms: the algorithms stop when the difference between two iterations is less than epsilon, by default 1e-08
#' @param eps_eigval Minimal value for eigen values (1e-04 by default)
#' @param eps_lambda Minimal value for estimated eigen values (1e-04 by default)
#' @param gamma Stochastic gradient parameter (0.75 by default). Must belong to (0,1).
#' @param c Stochastic gradient parameter (\eqn{\sqrt{ncol(X)}}) by default). Must be positive.
#' @param gammaRMC A parameter of the function RobbinsMC.
#' @param cRMC A parameter of the function RobbinsMC.
#' @param wRMC A parameter of the function RobbinsMC.

#'
#' @return A list with the following components:
#' \describe{
#'   \item{variance}{A robust estimate of the variance.}
#'   \item{outlier_labels}{If required, a binary vector of length \code{nrow(X)} indicating outliers (1 = outlier, 0 = inlier). Outlier detection is based on Mahalanobis distance.}
#'   \item{distances}{If required, a numeric vector of length \code{nrow(X)} giving the Mahalanobis distances of each observation from the estimated center.}
#'  }
#'
#' @md
#'
#' @examples
#' N <- 1000
#' p <- 10
#' X <- matrix(rnorm(N*p),ncol=p)
#' onlineRobustVariance(X=X)
#'
#' @importFrom robustbase covComed
#'
#' @export
onlineRobustVariance <- function(X,
                                 nDataInit=100,
                                 model='Gaussian',df=NULL,
                                 computeOutliers=FALSE,
                                 cutoff=0.95,cutoff_method="Chi-square",
                                 cutinit=0.6,c_m=2,
                                 batch=ncol(X),
                                 niterWeisz=100,
                                 sizeMC=batch,
                                 epsilon=1e-8,
                                 eps_eigval=1e-04,
                                 eps_lambda=1e-04,
                                 gamma=0.75,c=sqrt(ncol(X)),
                                 gammaRMC=0.75,cRMC=ncol(X),wRMC=2
)
{
  # Parameters check
  checkmate::assertChoice(model, c("Gaussian", "Student", "Laplace"))
  checkmate::assertInt(nDataInit, lower=0, upper=nrow(X))
  checkmate::assertInt(batch, lower=1, upper=ncol(X))
  checkmate::assertInt(niterWeisz, lower=1)
  checkmate::assertInt(sizeMC, lower=1)
  checkmate::assertNumber(epsilon, lower=0, finite=TRUE)
  checkmate::assertNumber(eps_eigval, lower=0, finite=TRUE)
  checkmate::assertNumber(eps_lambda, lower=0, finite=TRUE)
  if (model == "Student"){
    checkmate::assertInt(df, lower = 3)
  }
  if (computeOutliers){
    checkmate::assertChoice(model, c("Gaussian")) ###pas utile? ### Voir avec Daphné
    checkmate::assertChoice(cutoff_method, c("Chi-square","quantile" ))
    if (cutoff_method=="Chi-square")
      {
      cutoff= stats::qchisq(p=cutoff,df=ncol(X))
      }
  }

  # Initialization
  n <- nrow(X)
  d <- ncol(X)

  # Initialisation (in case nDataInit = 0)
  lambdatilde <- rep(1,d) # initialization of the estimates of the eigenvalues of the Variance
  m <- rep(0,d)
  averaged_m <- rep(0,d)
  V <- diag(1,d)
  averaged_V <- diag(1,d)

  distances <- NULL
  outlier_labels <- NULL
  varian <- NULL

  if (computeOutliers){
    distances <- rep(0,n) #initialization of a vector containing all the Mahalanobis distances
    outlier_labels <- rep(0,n) # initialization of a vector containing all the outlier labels (1 if outlier, 0 else).
  }

  if (nDataInit > 0){
    init_median <- WeiszfeldMedian(X[1:nDataInit,],nitermax=niterWeisz,epsilon=epsilon)
    init_median_cov <- robustbase::covComed(X[1:nDataInit,])$cov
    init_V <- WeiszfeldMedianCovariance(X[1:nDataInit,],median_est=init_median,init_median_cov=init_median_cov,nitermax=niterWeisz,epsilon=epsilon)

    init_eigen <- eigen(init_V,symmetric = TRUE)
    eigval_V <- init_eigen$values
    eigval_V <- pmax(eigval_V,eps_eigval)
    eigvec_V <- init_eigen$vectors

    if(model == 'Gaussian'){
      UU <- matrix(rnorm(d*nDataInit*10),ncol=d)
    }
    else if(model == 'Student'){
      UU <- matrix(rnorm(d*nDataInit*10)/sqrt(rchisq(1,df=df))*sqrt(df-2), ncol=d)
    }
    else{ # model == 'Laplace'
      UU <- LaplacesDemon::rmvl(nDataInit*10,mu=rep(0,d),Sigma=diag(d))
    }
    lambdaResultat <- robbinsMCList(U=UU,delta=eigval_V,gamma=gammaRMC,c=cRMC,w=wRMC,epsilon=epsilon)

    lambda <- lambdaResultat$lambda_bar
    lambdatilde <- lambdaResultat$lambda
    lambda=pmax(lambda,eps_lambda)
    lambdatilde=pmax(lambdatilde,eps_lambda)
    m <- init_median  # estimate of the median
    V <- init_V # estimate of the Median Covariation Matrix
    averaged_m <- init_median # averaged estimates
    averaged_V <- init_V # averaged estimates

    VP <- normalizeColumnsRcpp(eigvec_V)

    # Outliers detection
    if (computeOutliers){
      for (i in 1 :nDataInit) {
        distances[i] <- mahalanobisGeneralizedRcpp(X[i,],averaged_m,eigvec_V,lambda)
      }

      meddist <- stats::median(distances[1:nDataInit])
      if (cutoff_method=="Chi-square")
      {
        cutoffcor <- cutoff*meddist/stats::qchisq(.5,df = d)
      }
      if (cutoff_method=="quantile")
      {
        cutoffcor=quantile(distances[1:nDataInit],probs=cutinit)
        c_med=median(distances[1:nDataInit])
        for (i in 1 :nDataInit) {
        cutoffcor <- cutoffcor - c_m*c_med*(i+1)^(-0.75)*( as.numeric(distances[i] <= cutoffcor) - cutoff)
        }
      }

      for (i in 1:nDataInit){
        if (distances[i] > cutoffcor) {
          outlier_labels[i] <- 1
        }
      }
    }
  }

  slog <- 1 # To calculate the weights in RobbinsMC
  sslog <- 1 # To calculate the weighs in RobbinsMC

  niterr <- 0

  # Main loop
  for (i in (nDataInit+1):n){
    niterr <- niterr+1

    # If the current batch is the last one
    if ((niterr*batch+nDataInit >= n)  ){
      last_batch_size <- n - (nDataInit+(niterr-1)*batch)
      if (last_batch_size <= 0) break

      # Compute outliers for the last (incomplete) batch
      if (computeOutliers){
        for (j in 1:last_batch_size){
          S <-  mahalanobisGeneralizedRcpp(X[nDataInit + (niterr-1)*batch + j,],averaged_m,eigvec_V,lambda)

          distances[nDataInit + (niterr-1)*batch + j] <- S
          meddist <- stats::median(distances[1:nDataInit])
          if (cutoff_method=="Chi-square")
          {
            meddist <- meddist - (nDataInit + (niterr-1)*batch + j)^(-0.75)*( as.numeric(S <= meddist) - 0.5)
            cutoffcor <- cutoff*meddist/stats::qchisq(.5,df = d)
          }
          if (cutoff_method=="quantile")
          {
            cutoffcor <- cutoffcor - c_m*c_med*(nDataInit + (niterr-1)*batch + j)^(-0.75)*( as.numeric(distances[i] <= cutoffcor) - cutoff)
          }
          if (distances[nDataInit + (niterr-1)*batch + j] > cutoffcor){
            outlier_labels[nDataInit + (niterr-1)*batch + j] = 1
          }
        }
      }

      break
    }

    gamman <- c*sqrt(batch)/(i)^(gamma)

    upd <- updateMedianCovarianceRcpp(X, m, averaged_m, V, averaged_V, nDataInit, niterr, batch, gamman, wRMC, sslog)
    m <- upd$m
    averaged_m <- upd$averaged_m
    V <- upd$V
    averaged_V <- upd$averaged_V

    reseig <- eigen(averaged_V,symmetric = TRUE)
    eigvec_V <- reseig$vectors
    eigval_V <- reseig$values
    eigval_V <- pmax(eigval_V, eps_eigval)

    if(model == 'Gaussian'){
      UU <- matrix(rnorm(sizeMC*d),ncol=d)
    }
    else if(model == 'Student'){
      UU <- matrix(rnorm(sizeMC*d)/sqrt(rchisq(1,df=df))*sqrt(df-2), ncol=d)
    }
    else{ # model == 'Laplace'
      UU <- LaplacesDemon::rmvl(sizeMC,mu=rep(0,d),Sigma=diag(d))
    }
    lambdaResultat <- robbinsMCList(U=UU,delta=eigval_V,init=lambdatilde,init_bar=lambda,
                                    gamma=gammaRMC,c=cRMC,w=wRMC,
                                    c_bar=sizeMC*(niterr-1),c_tilde=sizeMC*(niterr-1),
                                    sumlog=sum((log(1:((sizeMC*(niterr-1))+1))^wRMC)),
                                    epsilon=epsilon)

    lambda <- lambdaResultat$lambda_bar
    lambda <- pmax(lambda,eps_lambda)
    lambdatilde <- lambdaResultat$lambda
    lambdatilde <- pmax(lambdatilde,eps_lambda)

    VP <- normalizeColumnsRcpp(eigvec_V)
    varian <- buildCovarianceRcpp(VP,lambda)

    if (computeOutliers){
      for (j in 1:batch){
        S <-  mahalanobisGeneralizedRcpp(X[nDataInit + (niterr-1)*batch + j,],averaged_m,eigvec_V,lambda)

        distances[nDataInit + (niterr-1)*batch + j] <- S
        if (cutoff_method=="Chi-square")
        {
          meddist <- meddist - (nDataInit + (niterr-1)*batch + j)^(-0.75)*( as.numeric(S <= meddist) - 0.5)
          cutoffcor <- cutoff*meddist/stats::qchisq(.5,df = d)
        }
        if (cutoff_method=="quantile")
        {
          cutoffcor <- cutoffcor - c_m*c_med*(nDataInit + (niterr-1)*batch + j)^(-0.75)*( as.numeric(distances[i] <= cutoffcor) - cutoff)
        }
        if (distances[nDataInit + (niterr-1)*batch + j] > cutoffcor){
          outlier_labels[nDataInit + (niterr-1)*batch + j] = 1
        }
      }
    }
  }

  return (list(variance=varian,outlier_labels=outlier_labels,distances=distances))
}


#' @title Offline robust estimation of the variance
#'
#' @description Given a data set X, this function estimates
#' robustly variance-covariance matrix via the geometric median and the Median Covariation Matrix.
#'
#' @param X Matrix of size n*p corresponding to the data. The rows are the observations.
#' @param methodMC Method to estimate the eigenvalues of the variance: can be "FixMC", "GradMC" and "RobbinsMC" (by default).
#' @param methodMCM Method for estimating the Median Covariation Matrix: can be "Weiszfeld" (by default) or "ASG".
#' @param model Type of distribution: can be "Gaussian" (by default), "Student" or "Laplace"
#' @param df Number of degrees of freedom if model='Student'.
#' @param computeOutliers Flag for outliers computation, FALSE by default.
#' @param cutoff Threshold at which a data item is considered contaminated (if computeOutliers is TRUE). Default is 0.95.
#' @param cutoff_method The cutoff method used to flag as outliers. Can be "Chi-square" if the distances are considered to follow a Chi-square distribution, or "quantile" if the estimated quantile of order "cutoff" of the distance is used.
#' @param cutinit A parameter giving the order of the quantile used for initializing the estimate of quantile of the distance if "quantile" method is used. Default is 0.6.
#' @param c_m The constant in the stepsequence of Rbobbins-Monro algorithm for estimating the quantile of the distances. Default is 2.
#' @param init_median A row vector for initializing the estimates of the median. Default is rep(0,ncol(X)).
#' @param init_median_cov A matrix for initializing the estimates of the MCM. Default is Identity.
#' @param niterWeisz The maximum number of iteration for the Weiszfeld algorithm. Default is 50.
#' @param niterMC Maximum number of iterations is MethodMC='FixMC' or 'GradMC'.
#' @param sizeMC Number of data generated in the Monte Carlo procedure to estimate the eigenvalues of the variance. Default is batch.
#' @param epsMCM Stopping criterion for Weiszfeld algorithm: the algorithm stops when the difference between two iterations is less than epsilon, by default 1e-08
#' @param epsMC Stopping criterion for MethodMC algorithm: the algorithm stops when the difference between two iterations is less than epsilon, by default 1e-08
#' @param eps_eigval Minimal value for eigen values (1e-04 by default)
#' @param eps_lambda Minimal value for estimated eigen values (1e-04 by default)
#' @param gammaMCM A positive constant if methodMCM='ASG'.
#' @param alphaMCM A positive constant if methodMCM='ASG'. Must belong to (0,1).
#' @param nstart Number of run if methodMCM='ASG'.
#' @param gammaRMC Robbins-Monro parameter for eigen values estimation, if methodMCM="RobbinsMC".
#' @param cRMC Robbins-Monro parameter for eigen values estimation, if methodMCM="RobbinsMC".
#' @param wRMC Robbins-Monro parameter for eigen values estimation, if methodMCM="RobbinsMC".
#'
#' @return A list with the following components:
#' \describe{
#'   \item{variance}{A robust estimate of the variance.}
#'   \item{outlier_labels}{If required, a binary vector of length \code{nrow(X)} indicating outliers (1 = outlier, 0 = inlier). Outlier detection is based on Mahalanobis distance.}
#'   \item{distances}{If required, a numeric vector of length \code{nrow(X)} giving the Mahalanobis distances of each observation from the estimated center.}
#' }
#'
#' @md
#'
#' @examples
#' N <- 1000
#' p <- 10
#' X <- matrix(rnorm(N*p),ncol=p)
#' offlineRobustVariance(X=X)
#'
#' @export
offlineRobustVariance <- function(X,methodMC='RobbinsMC', methodMCM='Weiszfeld',
                                  model='Gaussian',df=NULL,
                                  computeOutliers=FALSE,
                                  cutoff=0.95,cutoff_method="Chi-square",
                                  cutinit=0.6,c_m=2,
                                  init_median=rep(0,ncol(X)),init_median_cov=diag(ncol(X)),
                                  niterWeisz=50,niterMC=50,
                                  sizeMC=1000,
                                  epsMCM=1e-08,epsMC=1e-08,
                                  eps_eigval=1e-04,eps_lambda=1e-04,
                                  gammaMCM=2, alphaMCM=0.75,nstart=1,
                                  gammaRMC=0.75,cRMC=2,wRMC=2)
{
  checkmate::assertChoice(methodMC, c("RobbinsMC", "GradMC", "FixMC"))
  checkmate::assertChoice(methodMCM, c("Weiszfeld",  "ASG"))
  checkmate::assertChoice(model, c("Gaussian", "Student", "Laplace"))
  if (computeOutliers){
    checkmate::assertChoice(model, c("Gaussian")) ###pas utile? ### Voir avec Daphné
    checkmate::assertChoice(cutoff_method, c("Chi-square","quantile" ))
    if (cutoff_method=="Chi-square")
    {
      cutoff= stats::qchisq(p=cutoff,df=ncol(X))
    }
  }

  n <- nrow(X)
  p <- ncol(X)

  if(methodMCM=='Weiszfeld'){
    median <- WeiszfeldMedian(X=X, init_median=init_median, epsilon=epsMCM, nitermax=niterWeisz)
    medianCov <- WeiszfeldMedianCovariance(X=X, median_est=median, init_median_cov=init_median_cov, epsilon=epsMCM, nitermax=niterWeisz)
  }
  else if(methodMCM=='ASG')
  {
    median <- ASGMedian(X=X, init_median=init_median, gamma=gammaMCM, alpha=alphaMCM, nstart=nstart, epsilon=epsMCM)
    medianCov <- ASGMedianCovariance(X=X, median_est=median, init_median_cov=init_median_cov, gamma=gammaMCM, alpha=alphaMCM, nstart=nstart)
  }
  eig <- eigen(medianCov,symmetric = TRUE)
  eigenvec <- eig$vectors
  eigenval <- eig$values
  eigenval <- pmax(eigenval,eps_eigval)
  # Generate U for eigen values computation
  if(model == 'Gaussian'){
    U = matrix(rnorm(sizeMC*p),ncol=p)
  }
  else if(model == 'Student'){
    U <- matrix(rnorm(sizeMC*p)/sqrt(rchisq(1,df=df))*sqrt(df-2), ncol=p)
  }
  else{ # model == 'Laplace'
    U <- LaplacesDemon::rmvl(sizeMC,mu=rep(0,p),Sigma=diag(p))
  }

  # Compute eigen values of variance-covariance matrix using eigen values of median-covariance matrix
  if (methodMC=="RobbinsMC"){
    eigen_est <- robbinsMC(U=U,delta=eigenval,gamma=gammaRMC,c=cRMC,w=wRMC,epsilon=epsMC)
  }
  else if (methodMC=="GradMC"){
    eigen_est <- gradMC(U=U,delta=eigenval,niter=niterMC,epsilon=epsMC,step=cRMC*(1:niterMC)^(-0.5))
  }
  else if (methodMC=="FixMC"){
    eigen_est <- fixMC(U=U,delta=eigenval,niter=niterMC,epsilon=epsMC)
  }
  eigen_est = pmax(eigen_est,eps_lambda)
  variance <- t(matrix(eigenvec,ncol=p,byrow=TRUE))%*%diag(c(eigen_est))%*%(matrix(eigenvec,ncol=p,byrow=TRUE))

  distances = NULL
  outlier_labels = NULL

  if (computeOutliers){
    distances <- rep(0,n)
    outlier_labels <- rep(0,n)

    # Distances
    for (i in 1:n){
      S <-  mahalanobisGeneralizedRcpp(X[i,], median, eigenvec,eigen_est)
      distances[i] <- S
    }

    # Outliers detection
    if (cutoff_method=="Chi-square")
    {
      cutoffcor = cutoff*stats::median(distances)/stats::qchisq(.5,df = p)
    }
    if (cutoff_method=="quantile")
    {
      cutoffcor=quantile(distances,probs=cutinit)
      c_med=median(distances)
      for (i in 1 :n) {
        cutoffcor <- cutoffcor - c_m*c_med*(i+1)^(-0.75)*( as.numeric(distances[i] <= cutoffcor) - cutoff)
      }
    }
    for (i in 1:n){
      if (distances[i] > cutoffcor) {
        outlier_labels[i] = 1
      }
    }
  }

  return(list(variance=variance,outlier_labels=outlier_labels,distances=distances))
}
