#' Robust Mixture Model
#'
#' Robust EM algorithm for Mixture Model for a fixed number of clusters K
#'
#' @param X A numerical matrix giving the data of size nxd
#' @param K Integer giving the number of clusters
#' @param ninit The number of random initializations (10 by default)
#' @param initClust Initial clustering, typically an output of `mclust::hclass` with a single number of clusters or `genieclust::genie` (`NULL` by default)
#' @param methodMC Method to compute the eigen values and vectors of variance-covariance matrix starting from median covariation matrix, can be "FixMC", "GradMC" or "RobbinsMC" (by default)
#' @param methodMCM Method to compute the geometric median and the median covariation matrix, can be "Weiszfeld" (by default) or "ASG"
#' @param model Type of distribution: can be "Gaussian" (by default), "Student" or "Laplace"
#' @param df Degrees of freedom used when model is 'Student' (`NULL` by default). When used, it must be greater than 3
#' @param niterWiesz Maximum number of iterations in the geometric median computation using the Weiszfeld algorithm (50 by default)
#' @param niterEM Maximum number of updates in the EM algorithm (50 by default)
#' @param niterMC Maximum number of updates in the gradient descent algorithm in the estimation of the eigen values of the variance-covariance matrix (50 by default)
#' @param sizeMC Size of the sample for the Monte Carlo method in the estimation of the eigen values of the variance-covariance matrix (1000 by default)
#' @param logLikeInit Initial value of log likelihood (-1e10 by default).
#' @param eps_dist Distance between old and new centers after each iteration has to be greater than eps_dist (1e-04 by default)
#' @param eps_eigval Minimal value for eigen values (1e-10 by default)
#' @param eps_MCM_MC Stopping criterion used in the Weiszfeld geometric median computation and in the estimation of the eigen values of the variance-covariance matrix
#' @param eps_tau Small correction to avoid 0 and 1 values of tau (1e-04 by default)
#' @param logPhiTol Tolerance for small values of log densities (-20 by default). This must be a negative value, since the log densities are first adjusted by substracting their max.
#' @param outlierThreshold Threshold for outliers detection, on a log density criterion (-20 by default)
#' @param gammaRMC Robbins-Monro parameter for eigen values estimation, if methodMC="RobbinsMC"
#' @param cRMC Robbins-Monro parameter for eigen values estimation, if methodMC="RobbinsMC". When parameter `scale=TRUE`, `cRMC=1` should be a good choice.
#' @param wRMC Robbins-Monro parameter for eigen values estimation, if methodMC="RobbinsMC"
#' @param scale Boolean value for scaling data before running the algorithm (`FALSE` by default)
#'
#' @return A list with the following components:
#' \describe{
#'   \item{centers}{A numeric matrix of size \code{K × ncol(X)}, where each row corresponds to the center of a class.}
#'   \item{Sigma}{A numeric array of size \code{K × ncol(X) × ncol(X)} containing the covariance matrices for each class.}
#'   \item{LogLike}{The log-likelihood of the selected model.}
#'   \item{entropy}{The entropy of the selected model.}
#'   \item{tau}{A numeric matrix of size \code{n × K}, where each row contains the posterior probabilities of class membership for each observation.}
#'   \item{clusters}{An integer vector of length \code{nrow(X)} indicating the assigned class for each observation.}
#'   \item{niter}{The total number of iterations performed.}
#'   \item{initEM}{The initialization used by the EM algorithm.}
#'   \item{prop}{A numeric vector of length \code{K} giving the estimated proportion of each class.}
#'   \item{outliers}{A binary vector of length \code{nrow(X)} where a value of 1 indicates that the corresponding observation is considered as an outlier. An observation is flagged as an outlier if its log-likelihood under all classes is below \code{outlierThreshold}.}
#' }
#'
#' @seealso \link{multipleRobustMM}
#'
#' @references \url{https://cnrs.hal.science/hal-03853744/}
#' @references Cardot, H., Cenac, P. and Zitt, P-A. (2013). Efficient and fast estimation of the geometric median in Hilbert spaces with an averaged stochastic gradient algorithm. \emph{Bernoulli}, 19, 18-43.
#' @references Cardot, H. and Godichon-Baggioni, A. (2017). Fast Estimation of the Median Covariation Matrix with Application to Online Robust Principal Components Analysis.  \emph{Test}, 26(3), 461-480
#' @references Vardi, Y. and Zhang, C.-H. (2000). The multivariate L1-median and associated data depth. \emph{Proc. Natl. Acad. Sci. USA}, 97(4):1423-1426.
#'
#' @export
#' @md
robustMM <- function(X,K=2,ninit=10,initClust=NULL,
                  methodMC="RobbinsMC",methodMCM="Weiszfeld",model="Gaussian",df=NULL,
                  niterWiesz=50,niterEM=50,niterMC=50,
                  sizeMC=1000,
                  logLikeInit=-1e10,eps_dist=1e-04,eps_eigval=1e-10,eps_MCM_MC=1e-03,eps_tau=1e-04,logPhiTol=-20,outlierThreshold=-20,
                  gammaRMC=0.75,cRMC=ncol(X),wRMC=2,
                  scale=FALSE
                  )
{
  # Arguments check
  checkmate::assertChoice(methodMC, c("RobbinsMC", "GradMC", "FixMC"))
  checkmate::assertChoice(methodMCM, c("Weiszfeld", "Weiszfeld_init", "ASG", "ASG_init"))
  checkmate::assert_logical(scale)
  checkmate::assertChoice(model, c("Gaussian", "Student", "Laplace"))
  if (model == "Student"){
    checkmate::assertInt(df, lower = 3)
  }
  if (ninit <= 0 && is.null(initClust))
    stop("Either ninit>0 or initClust != NULL")

  # Robust scaling with median and MAD (median absolute deviation) if asked
  if (scale==TRUE){
    X <- DescTools::RobScale(X)
    robCenters <- attr(X,"scaled:center")
    robScales <- attr(X,"scaled:scale")
    cRMC <- 1
  }

  # Initialisation
  d <- ncol(X)
  n <- nrow(X)

  finalcenters <- matrix(0,nrow=K,ncol=d)
  finalcluster <- matrix(0,nrow=n,ncol=K)
  logPhi <- matrix(0,nrow=n,ncol=K) # log densities \log \phi_ik
  finalPhiik <- matrix(0,nrow=n,ncol=K)
  finalvar <- matrix(0,nrow=d,ncol=K*d)
  finalprop <- rep(0,K) # pi_k
  finalniter <- 0
  finaloutliers <- c()

  currentLogLike <- logLikeInit

  # If an initial clustering was given
  if (!is.null(initClust)){
    # initial proportion of individuals in each cluster
    prop <- sapply(1:K, function(k) sum(initClust==k))
    prop <- prop/sum(prop)

    # tau is a probability but we start with a one hot encoding of the given clustering
    # (probability 1 of being in a cluster and 0 to be in all the others)
    tau <- matrix(0,nrow=n,ncol=K)
    for (k in 1:K){
      I <- which(initClust==k)
      tau[I,k] <- 1
    }

    centers <- matrix(0,ncol=d,nrow=K) # pourquoi pas aléatoire dans X comme plus bas ? parce que j'ai une init et que j'ai des premiers tau que je vais donner au calcul de mediane

    # Initialisation of Sigma and V
    Sigma <- matrix(rep(diag(d),K),nrow=d)
    V <- matrix(rep(diag(d),K),nrow=d)

    l <- 0

    dist <- 2*eps_dist
    epsout <- -20

    while (l < niterEM && dist > eps_dist){
      l <- l+1
      old_centers <- centers

      # Update centers, variances and weights
      res <- update(X,K,tau,centers,Sigma,V,
                    methodMC,methodMCM,model,df,
                    niterWiesz,niterMC,sizeMC,
                    eps_MCM_MC,eps_eigval,
                    gammaRMC,cRMC,wRMC)

      #res <- update(K , d, tau, methodMCM, X, centers, Sigma, V, niterWiesz, eps_MCM_MC, methodMC, sizeMC, gammaRMC, wRMC, cRMC, niterMC, model, df, eps_eigval)
      centers <- res$centers
      Sigma <- res$Sigma
      V <- res$V
      prop <- res$prop

      dist <- sum((old_centers-centers)^2)/(K*d)

      # Compute logdensity logPhi
      for (k in 1 : K){
        centers_k <- centers[k,,drop=FALSE]
        var_k <- Sigma[,((k-1)*d+1):(k*d)]

        if (any(eigen(var_k)$values<0))
          cat('var_k has negative eigen values : l=',l,'\n')

        if (model=="Gaussian")
          logPhi[,k] <- mvtnorm::dmvnorm(X,mean=centers_k,sigma=var_k,log=TRUE)
        else if(model=="Student")
          logPhi[,k] <- mvtnorm::dmvt(X,delta=centers_k,sigma=(df-2)/df*var_k,df=df,log=TRUE)
        else # model=="Laplace"
          logPhi[,k] <- LaplacesDemon::dmvl(X,mu=centers_k,Sigma=var_k,log=TRUE)
      }

      if (any(is.infinite(logPhi))){
        cat('logPhi has infinity values : l=',l,'\n')
      }

      # Update tau with proper adjustments
      tau <- logPhi

      # Avoid exp to be zero or NaN for at least one element
      # This does not change final tau since it is equivalent to multiply and divide by the same constant (tau will be later divided by rowSums(tau))
      # tau - max(tau, axis=1)
      tau <- tau - apply(tau,1,max)

      # Avoid too small values (change on final tau, but only for values which contribute less than exp(logPhiTol) in the sum)
      tau[tau<logPhiTol] <- logPhiTol

      # Effective tau update
      for (k in 1:K)
        tau[,k] <- prop[k]*exp(tau[,k] ) # Probabilities of belonging to cluster k
      tau <- tau/rowSums(tau) # tau is a probability

      # Check for NaN values
      if (any(is.na(tau))){
        cat('tau has NaN values : l=',l,'\n')
      }

      # Avoid tau to be 0 or 1, to avoid tau getting stuck
      tau <- tau+eps_tau
      tau <- tau/rowSums(tau)
    } # end of while (l < niterEM && dist > eps_dist)

    # Detect outliers
    outliers <- c()
    max_logPhi <- apply(logPhi,1,max)
    outliers <- which(max_logPhi<=outlierThreshold) # this criteria could be improved (PhD Paul)

    # Avoid log densities to be too small
    logPhi[logPhi<logPhiTol] <- logPhiTol # TODO garder ou pas

    entropy <- - sum(tau*log(tau + (tau==0)))
    LogLikeEM <- sum(tau*logPhi)+ sum(tau*log(prop)) + entropy

    if (all(!is.na(LogLikeEM)) && LogLikeEM > currentLogLike){
      finalcenters <- centers
      finalcluster <- tau
      finalvar <- Sigma
      finalPhiik <- logPhi
      currentLogLike <- LogLikeEM
      finalniter <- l
      finalprop <- prop
      finaloutliers <- outliers
    }
  } # end of if (!is.null(initClust))

  if (ninit > 0){
    for (o in 1:ninit){
      # Initialisation
      prop <- rep(1/K,K)

      # Tau initialisation
      tau <- matrix(1,nrow=n,ncol=K)

      # Centers, Sigma and V initialisation
      centers <- as.matrix(X[sample(1:n,K),, drop=FALSE],ncol=d) # centers initialisation with random points of X
      Sigma <- matrix(rep(diag(d),K),nrow=d)
      V <- matrix(rep(diag(d),K),nrow=d)

      l <- 0
      dist <- 2*eps_dist

      # EM
      while (l < niterEM && dist > eps_dist)
      {
        l <- l+1

        # Update probabilities
        for (k in 1 : K){
          centers_k <- centers[k,, drop=FALSE]
          var_k <- Sigma[,((k-1)*d+1):(k*d)]

          if (any(eigen(var_k)$values<0))
            cat('var_k has negative eigen values : l=',l,'k=', k,'\n')

          if (model=="Gaussian")
            logPhi[,k] <- mvtnorm::dmvnorm(X,mean=centers_k,sigma=var_k,log=TRUE)
          else if(model=="Student")
            logPhi[,k] <- mvtnorm::dmvt(X,delta=centers_k,sigma=(df-2)/df*var_k,df=df,log=TRUE)
          else # model=="Laplace"
            logPhi[,k] <- LaplacesDemon::dmvl(X,mu=centers_k,Sigma=var_k,log=TRUE)
        }

        if (any(is.infinite(logPhi))){
          cat('logPhi has infinity values : l=',l,'\n')
        }

        # Update tau, with proper adjustments
        tau <- logPhi

        # Avoid exp to be zero or NaN for at least one element
        tau <- tau - apply(tau,1,max)

        # Avoid too small values
        tau[tau<logPhiTol] <- logPhiTol

        for (k in 1:K)
          tau[,k] <- prop[k]*exp(tau[,k] )
        tau <- tau/rowSums(tau)

        tau <- tau+eps_tau
        tau <- tau/rowSums(tau)

        old_centers <- centers

        # Check for NaN values
        if (any(is.na(tau))){
          cat('tau has NaN values : l=',l,'\n')
        }

        # Update centers, variances and weights
        res <- update(X,K,tau,centers,Sigma,V,
        methodMC,methodMCM,model,df,
        niterWiesz,niterMC,sizeMC,
        eps_MCM_MC,eps_eigval,
        gammaRMC,cRMC,wRMC)

        #res <- update(K , d, tau, methodMCM, X, centers, Sigma, V, niterWiesz, eps_MCM_MC, methodMC, sizeMC, gammaRMC, wRMC, cRMC, niterMC, model, df, eps_eigval)
        centers <- res$centers
        Sigma <- res$Sigma
        V <- res$V
        prop <- res$prop

        dist <- sum((old_centers-centers)^2)/(K*d)
      } # end of while (l < niterEM && dist > eps_dist)

      # Detect outliers
      outliers <- c()
      max_logPhi <- apply(logPhi,1,max)
      outliers <- which(max_logPhi<=outlierThreshold)

      # Avoid log densities to be too small
      logPhi[logPhi<logPhiTol] <- logPhiTol

      entropy <- - sum(tau*log(tau + (tau==0)))
      LogLikeEM <- sum(tau*logPhi)+ sum(tau*log(prop)) + entropy

      # Update if better criteria
      if (all(!is.na(LogLikeEM)) && LogLikeEM > currentLogLike){
        finalcenters <- centers
        finalcluster <- tau
        finalvar <- Sigma
        currentLogLike <- LogLikeEM
        finalniter <- l
        finalPhiik <- logPhi
        finalprop <- prop
        finaloutliers <- outliers
      }
    } # end of for (o in 1:ninit)
  } # end of if (ninit > 0)

  # convert back with scaling (unscale output)
  if (scale==TRUE){
    for (k in 1:K){
      finalvar[,((k-1)*d+1):(k*d)] <- diag(robScales)%*%finalvar[,((k-1)*d+1):(k*d)]%*%diag(robScales)
      finalcenters[k,] <- finalcenters[k,]%*%diag(robScales)+ robCenters
    }
  }

  clusters <- apply(finalcluster,1,which.max)

  return(list(centers=finalcenters,Sigma=finalvar,Loglike=currentLogLike,entropy=entropy, tau=finalcluster,clusters=clusters,niter=finalniter,initEM=initClust,prop=finalprop,outliers=finaloutliers))
}


#' @keywords internal
update <- function(X,K,tau,centers,Sigma,V,
                   methodMC,methodMCM,model,df,
                   niterWiesz,niterMC,sizeMC,
                   eps_MCM_MC,eps_eigval,
                   gammaRMC,cRMC,wRMC){
  if (any(is.na(tau)))
    cat('tau has NaN values : K=',K,'\n')

  #print(paste0("K=", K))
  d <- ncol(X)
  prop <- rep(0, K)

  for (k in 1:K){
    #print(paste0("k=", k))
    if(any(tau[,k]>0)){
      prop[k] <- mean(tau[,k])
      if (methodMCM=="Weiszfeld_init"){
        med <- WeiszfeldMedian(X, init_median=centers[k,,drop=FALSE], weights=(tau[,k])/sum(tau[,k]), epsilon = eps_MCM_MC, nitermax = niterWiesz)
        medcovmat <- WeiszfeldMedianCovariance(X, median_est = med, init_median_cov=V[,((k-1)*d+1):(k*d)], weights=(tau[,k])/sum(tau[,k]), epsilon = eps_MCM_MC, nitermax = niterWiesz)
      }
      else if (methodMCM=="Weiszfeld"){
        med <- WeiszfeldMedian(X, weights=(tau[,k])/sum(tau[,k]), epsilon = eps_MCM_MC, nitermax = niterWiesz)
        medcovmat <- WeiszfeldMedianCovariance(X, median_est = med, weights=(tau[,k])/sum(tau[,k]), epsilon = eps_MCM_MC, nitermax = niterWiesz)
      }
      else if (methodMCM=="ASG_init"){
        med <- ASGMedian(X, init_median=centers[k,,drop=FALSE], weights=(tau[,k]))
        medcovmat <- ASGMedianCovariance(X, median_est = med, init_median_cov=V[,((k-1)*d+1):(k*d)], weights=(tau[,k]))
      }
      else if (methodMCM=="ASG"){
        med <- ASGMedian(X, weights=(tau[,k]))
        medcovmat <- ASGMedianCovariance(X, median_est = med, weights=(tau[,k]))
      }

      if (all(!is.na(medcovmat)) && all(is.finite(medcovmat))){
        centers[k,] <- med
        eig <- eigen(medcovmat)
        eigenvec <- eig$vectors  # vecteur propre de \hat V du papier
        eigenval <- eig$values  # delta du papier

        # Generate U for eigen values computation
        p <- length(eigenval)

        if(model == 'Gaussian'){
          U = matrix(rnorm(sizeMC*p),ncol=p)
        }
        else if(model == 'Student'){
          U <- matrix(rnorm(sizeMC*p)/sqrt(rchisq(1,df=df))*sqrt(df-2), ncol=p)
        }
        else{ # model == 'Laplace'
          U <- LaplacesDemon::rmvl(sizeMC,mu=rep(0,p),Sigma=diag(p))
        }

        # Compute eigen values of variance-covariance matrix using eigen values of median-covariation matrix by Monte Carlo
        if (methodMC=="RobbinsMC"){
          eigval_est <- robbinsMC(U=U,delta=eigenval,init=eigenval,gamma=gammaRMC,c=cRMC,w=wRMC,epsilon=eps_MCM_MC)
        }
        else if (methodMC=="GradMC"){
          eigval_est <- gradMC(U=U,delta=eigenval,niter=niterMC,epsilon=eps_MCM_MC,step=cRMC*(1:niterMC)^(-0.5))
        }
        else if (methodMC=="FixMC"){
          eigval_est <- fixMC(U=U,delta=eigenval,niter=niterMC,epsilon=eps_MCM_MC)
        }

        # Warn user if small eigen values adjustment
        if(any(eigval_est<eps_eigval))
          cat('Small eigen values are being replaced by ',eps_eigval,'\n')
        lambdak <- pmax(eigval_est,eps_eigval)

        # Update Sigma_k and V_k
        V[,((k-1)*d+1):(k*d)] <- medcovmat
        Sigma[,((k-1)*d+1):(k*d)] <- t(matrix(eigenvec,ncol=d,byrow=T))%*%diag(c(lambdak))%*%(matrix(eigenvec,ncol=d,byrow=T)) # calcul de \psi(\hat V)

      }
    }
  }
  return(list(centers=centers, Sigma=Sigma, V=V, prop=prop))
}

#' Robust Mixture Model
#'
#' Robust EM algorithm for Mixture Model for a number of clusters in a range
#'
#' @param X A numerical matrix giving the data of size nxd
#' @param nclust Range of number of clusters
#' @param ninit The number of random initializations (10 by default)
#' @param init Kind of initialization, can be 'Mclust' (default) or 'genie'
#' @param methodMC Method to compute the eigen values and vectors of variance-covariance matrix starting from median covariation matrix, can be "FixMC", "GradMC" or "RobbinsMC" (by default)
#' @param methodMCM Method to compute the geometric median and the median covariation matrix, can be "Weiszfeld" (by default) or "ASG"
#' @param model Type of distribution: can be "Gaussian" (by default), "Student" or "Laplace"
#' @param df Degrees of freedom used when model is 'Student' (NULL by default). When used, it must be greater than 3
#' @param niterWiesz Maximum number of iterations in the geometric median computation using the Weiszfeld algorithm
#' @param niterEM Maximum number of updates in the EM algorithm
#' @param niterMC Maximum number of updates in the gradient descent algorithm in the estimation of the eigen values of the variance-covariance matrix
#' @param sizeMC Size of the sample for the Monte Carlo method in the estimation of the eigen values of the variance-covariance matrix
#' @param logLikeInit Initial value of log likelihood
#' @param eps_dist Distance between old and new centers after each iteration has to be greater than eps_dist
#' @param eps_eigval Minimal value for eigen values (1e-10 by default)
#' @param eps_MCM_MC Stopping criterion used in the Weiszfeld geometric median computation and in the estimation of the eigen values of the variance-covariance matrix
#' @param eps_tau Small correction to avoid 0 and 1 values in tau computation (1e-04 by default)
#' @param logPhiTol Tolerance for small values of log densities (-20 by default)
#' @param outlierThreshold Threshold for outliers detection (-20 by default)
#' @param gammaRMC Robbins-Monro parameter for eigen values estimation, if methodMC="RobbinsMC"
#' @param cRMC Robbins-Monro parameter for eigen values estimation, if methodMC="RobbinsMC". When parameter scale=TRUE, cRMC=1 should be a good choice.
#' @param wRMC Robbins-Monro parameter for eigen values estimation, if methodMC="RobbinsMC"
#' @param scale Boolean value for scaling data before running the algorithm (FALSE by default)
#' @param criterion Criterion for choice of best clustering, can be 'ICL'  or 'BIC' (default)
#' @param nb_cores Number of cores in parallel computation, by default parallel::detectCores()
#'
#' @return An object of class \code{RMM}, which is a list with the following components:
#' \describe{
#'   \item{allresults}{A list containing the results from the \code{RobustMM} function for each number of classes specified in \code{nclust}.}
#'   \item{bestresult}{A list corresponding to the output of \code{RobustMM} for the best model selected using the BIC or ICL criterion.}
#'   \item{ICL}{A numeric vector of length \code{length(nclust)} containing the Integrated Complete-data Likelihood values for the different models.}
#'   \item{BIC}{A numeric vector of length \code{length(nclust)} containing the Bayesian Information Criterion values for the different models.}
#'   \item{data}{The original data matrix \code{X}.}
#'   \item{nclust}{An integer vector specifying the number of classes (models) considered.}
#'   \item{Kopt}{The optimal number of classes selected based on the BIC or ICL criterion.}
#' }
#'
#' @seealso \link{robustMM}
#'
#' @references \url{https://cnrs.hal.science/hal-03853744/}
#' @references Cardot, H., Cenac, P. and Zitt, P-A. (2013). Efficient and fast estimation of the geometric median in Hilbert spaces with an averaged stochastic gradient algorithm. \emph{Bernoulli}, 19, 18-43.
#' @references Cardot, H. and Godichon-Baggioni, A. (2017). Fast Estimation of the Median Covariation Matrix with Application to Online Robust Principal Components Analysis.  \emph{Test}, 26(3), 461-480
#' @references Vardi, Y. and Zhang, C.-H. (2000). The multivariate L1-median and associated data depth. \emph{Proc. Natl. Acad. Sci. USA}, 97(4):1423-1426.
#'
#' @export
#' @md
multipleRobustMM <- function(X,
                             nclust=2:5,
                             ninit=10,
                             init='Mclust',
                             methodMC="RobbinsMC",methodMCM="Weiszfeld",model="Gaussian",df=NULL,
                             niterWiesz=50,niterEM=50,niterMC=50,
                             sizeMC=1000,
                             logLikeInit=-1e10,eps_dist=1e-04,eps_eigval=1e-10,eps_MCM_MC=1e-03,eps_tau=1e-04,logPhiTol=-20,outlierThreshold=-20,
                             gammaRMC=0.75,cRMC=ncol(X),wRMC=2,
                             scale=FALSE,
                             criterion='BIC',
                             nb_cores = parallel::detectCores())
{
  # Arguments check
  checkmate::assertIntegerish(nclust)
  checkmate::assertChoice(methodMC, c("RobbinsMC", "GradMC", "FixMC"))
  checkmate::assertChoice(methodMCM, c("Weiszfeld", "Weiszfeld_init", "ASG", "ASG_init"))
  checkmate::assert_logical(scale)
  checkmate::assertChoice(model, c("Gaussian", "Student", "Laplace"))
  if (model == "Student"){
    checkmate::assertInt(df, lower = 3)
  }

  checkmate::assertChoice(init, c('Mclust', 'genie'), null.ok=TRUE)
  checkmate::assertChoice(criterion, c('ICL', 'BIC'))
  checkmate::assertInt(nb_cores, upper = parallel::detectCores())

  initClust <- NULL

  ICL <- c()
  BIC <- c()

  if (length(nclust)==1)
  {
    Kopt <- nclust

    if (init=='Mclust'){
      clas <- mclust::hcVVV(data=X)
      initClust <- mclust::hclass(clas,nclust)
    }
    # else init=='genie'
    else {
      initClust <- genieclust::genie(X,k=nclust)
    }

    resultat <- robustMM(X,K=nclust,ninit=ninit,initClust=initClust,
                         methodMC=methodMC,methodMCM=methodMCM,model=model,df=df,
                         niterWiesz=niterWiesz,niterEM=niterEM,niterMC=niterMC,
                         sizeMC=sizeMC,
                         logLikeInit=logLikeInit,eps_dist=eps_dist,eps_eigval=eps_eigval,eps_MCM_MC=eps_MCM_MC,eps_tau=eps_tau,logPhiTol=logPhiTol,outlierThreshold=outlierThreshold,
                         gammaRMC=gammaRMC,cRMC=cRMC,wRMC=wRMC,
                         scale=scale)

    bestresult <- resultat

    penality <- 0.5*log(nrow(X))*(nclust-1+ nclust*ncol(X) + nclust*ncol(X)*(ncol(X)+1)/2)
    ICL <- c(ICL, bestresult$Loglike - penality - bestresult$entropy)
    BIC <- c(BIC, bestresult$Loglike - penality)
  } # end of if (length(nclust)==1)

  else if (length(nclust)>1)
  {
    future::plan(future::multisession, workers = nb_cores)

    KK = NULL
    resultat <- foreach::foreach(KK=nclust, .options.future = list(seed = TRUE))  %dofuture%
    {
      if (init=='Mclust'){
        clas <- mclust::hcVVV(data=X)
        initClust <- mclust::hclass(clas,KK)
      } #else if (init=='genie')
      else {
        initClust <- genieclust::genie(X,k=KK)
      }
      #         cat('Running for : K=',K,'\n') TODO faire marcher ce print (si pas trop long)
      #          cat("Running K =", min(nclust), "...\n")

      resultatk <- robustMM(X,K=KK,ninit=ninit,initClust=initClust,
                           methodMC=methodMC,methodMCM=methodMCM,model=model,df=df,
                           niterWiesz=niterWiesz,niterEM=niterEM,niterMC=niterMC,
                           sizeMC=sizeMC,
                           logLikeInit=logLikeInit,eps_dist=eps_dist,eps_eigval=eps_eigval,eps_MCM_MC=eps_MCM_MC,eps_tau=eps_tau,logPhiTol=logPhiTol,outlierThreshold=outlierThreshold,
                           gammaRMC=gammaRMC,cRMC=cRMC,wRMC=wRMC,
                           scale=scale)
      return(resultatk)
    } # end of foreach::foreach(...) %dofuture%

    for (k in 1:length(nclust)){
      penality <- 0.5*log(nrow(X))*(nclust[k]-1+ nclust[k]*ncol(X) + nclust[k]*ncol(X)*(ncol(X)+1)/2)
      ICL <- c(ICL,resultat[[k]]$Loglike - penality - resultat[[k]]$entropy)
      BIC <- c(BIC,resultat[[k]]$Loglike - penality)
    }
    if (criterion=='ICL'){
      k <- which.max(ICL)
      Kopt <- nclust[k]
    }
    else if (criterion=='BIC'){
      k <- which.max(BIC)
      Kopt <- nclust[k]
    }
    bestresult <- resultat[[k]]
  }

  RMM_res <- list(allresults=resultat,bestresult=bestresult,ICL=ICL,BIC=BIC,data=X,nclust=nclust,Kopt=Kopt)

  class(RMM_res) <- "RMM"
  return(RMM_res)
}
