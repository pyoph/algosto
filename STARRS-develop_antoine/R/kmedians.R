#' K-medians algorithms
#'
#' @description Perform k-medians clustering on a data matrix.
#'
#' @param X A numerical matrix giving the data
#' @param nclust A vector of positive integers giving the possible numbers of clusters. Default is \code{1:15}.
#' @param ninit A non negative integer giving the number of random initializations. Default is \code{0}.
#' @param niter A positive integer giving the number of iterations for the EM algorithms. Default is \code{20}.
#' @param method The selected method for the K-medians algorithm. Can be
#' \itemize{
#'   \item \code{Offline}: Default. Medians are computed with Weizsfeld algorithm.
#'   \item \code{Semi-Online}: Medians are computed with ASG algorithm.
#'   \item \code{Online}: Medians are computed with Online algorithm.
#' }
#' @param init A logical argument telling if the function \code{'genie'} is used for initializing the algorithm. Default is \code{TRUE}.
#' @param nb_cores Number of cores in parallel computation, by default parallel::detectCores()
#'
#' @return List with the following components :
#'
#' \itemize{
#'   \item \code{bestresults}: A list giving all the results for the clustering selected by \code{'capushe'}.
#'   \itemize{
#'    \item \code{cluster} A vector of positive integers giving the clustering.
#'    \item \code{centers} A numerical matrix giving the centers of the clusters.
#'    \item \code{SE} An integer giving the Sum of Errors.
#'   }
#'   \item \code{allresults}: A list containing all the results.
#'   \item \code{SE}: A vector giving the Sum of Errors for each considered number of clusters.
#'   \item \code{cap}: The results given by the function \code{'capushe'} if \code{nclust} is  of length larger than \code{10}.
#'   \item \code{Kopt}: An integer giving the number of clusters selected by \code{'capushe'} if \code{nclust} is  of length larger than \code{10}.
#' }
#'
#' @references Godichon-Baggioni, A. and Surendran, S. A penalized criterion for selecting the number of clusters for K-medians. \emph{arxiv.org/abs/2209.03597}
#'
#' @examples
#' \donttest{
#'   n <- 1000
#'   K <- 3
#'   pcont <- 0.2
#'   sphereGM <- rSphereGM(n=n,K=K,pcont=pcont)
#'   X <-sphereGM$X
#'   res <- kmedians(X, nb_cores=1)
#'   plot(X, col=res$bestresult$cluster)
#' }
#'
#' @export
kmedians <- function(X,nclust=1:15,ninit=0,niter=20,method='Offline', init=TRUE, nb_cores=parallel::detectCores())
{
  checkmate::assertChoice(method,c("Offline", "Online", "Semi-Online"))

  if (length(nclust)==1)
  {
    resultat <- kmedian(X,K=nclust,ninit=ninit,niter=niter,init=init,method=method)
    return(list(bestresult=resultat,
                allresults=resultat,
                SE=resultat$SE,
                cap=NULL,
                Kopt=nclust))
  }

  future::plan(future::multisession, workers = nb_cores)

  KK = NULL
  resultat <- foreach::foreach(KK=nclust, .options.future = list(seed = TRUE))  %dofuture%
  {
    resultatk <- kmedian(X,K=KK,ninit=ninit,niter=niter,init=init,method=method)
    return(resultatk)
  }

  SE <- c()
  for (i in 1:length(nclust))
  {
    SE <- c(SE,resultat[[i]]$SE)
  }

  if (length(nclust)<= 10)
  {
    message('Capushe cannot be used, nclust must be of length larger that 10. No model has been selected.')
    return(list(bestresult=NULL,
                allresults=resultat,
                SE=SE,
                cap=NULL,
                Kopt=NULL))
  }

  cap <- capushe::capushe(cbind(nclust,sqrt(nclust),nclust,SE))
  Kopt <- as.numeric(cap@DDSE@model)
  I <- which(nclust==Kopt)
  bestresult <- resultat[[I]]

  return(list(bestresult=bestresult,
              allresults=resultat,
              SE=SE,
              cap=cap,
              Kopt=Kopt))
}

#' K-medians algorithm for a fixed number of clusters
#'
#' @keywords internal
kmedian <- function(X,K=3,ninit=0,niter=20,init=TRUE,method='Offline')
{
  checkmate::assertInt(K, lower=1, upper=nrow(X))

  if (method=='Online')
  {
    resultat <- kmedianOnline(X,K=K,ninit=ninit,niter=niter,init=init)
  }
  else
  {
    resultat <- kmedianOffline(X,K=K,ninit=ninit,niter=niter,init=init,method=method)
  }

  return(resultat)
}

#' K-medians algorithm for a fixed number of clusters
#' This function uses the Online algorithm
#' for computing the median of a set of points.
#'
#' @keywords internal
kmedianOnline <- function(X,K=3,ninit=0,niter=20,init=TRUE)
{
  d <- ncol(X)
  n <- nrow(X)

  if (K==1)
  {
    centers <- ASGMedian(X)
    finalcenters <- centers
    finalcluster <- rep(1,n)
    finaldist <- sqrt(rowSums((X-matrix(rep(centers,n),nrow=n,byrow=TRUE))^2))
    finalSE <- sum(finaldist)

    return(list(cluster=finalcluster, centers=finalcenters, SE=finalSE))
  }

  centers <- matrix(0,nrow=K,ncol=d)
  finalcenters <- centers
  dist <- matrix(0,nrow=n,ncol=K)
  finaldist <- dist
  cluster <- rep(0,n)
  finalcluster <- cluster
  finalSE <- 10^10

  if (init==TRUE)
  {
    cluster <- genieclust::genie(X,k=K)
    for (k in 1:K)
    {
      I <- which(cluster==k)
      if (length(I)>1)
      {
        centers[k,] <- WeiszfeldMedian(X[I,])
      }
      if (length(I) == 1)
      {
        centers[k,] <- X[I,]
      }
    }
    resultat <- kmedianClustering(X,ncenters=centers)

    SE <- 0
    for (k in 1:K)
    {
      dist[,k] <- sqrt(rowSums((X-matrix(rep(resultat$centers[k,],n),nrow=n,byrow=T))^2))
      I <- which(resultat$cluster==k)
      SE <- SE+sum(dist[I,k])
    }

    finalcluster <- resultat$cluster
    finalcenters <- resultat$centers
    finaldist <- dist
    finalSE <- SE
  }

  if (ninit>0)
  {
    for (o in 1:ninit)
    {
      dist <- matrix(0,nrow=n,ncol=K)
      centers <- X[sample(1:n,K),]
      resultat <- kmedianClustering(X,ncenters=centers)

      SE <- 0
      for (k in 1:K)
      {
        dist[,k] <- sqrt(rowSums((X-matrix(rep(resultat$centers[k,],n),nrow=n,byrow=T))^2))
        I <- which(resultat$cluster==k)
        SE <- SE+sum(dist[I,k])
      }
      if (SE < finalSE)
      {
        finalcenters <- resultat$centers
        finalcluster <- resultat$cluster
        finaldist <- resultat$withinsrs
        finalSE <- sum(resultat$withinsrs)
      }
    }
  }

  return(list(cluster=finalcluster, centers=finalcenters, dist=finaldist, SE=finalSE))
}

#' K-medians algorithm for a fixed number of clusters
#' This function uses the Weiszfeld algorithm (Offline) or the ASG algorithm (Semi-Online)
#' for computing the median of a set of points.
#'
#' @keywords internal
kmedianOffline <- function(X,K=3,ninit=0,niter=20,init=TRUE, method='Semi-Online')
{
  d <- ncol(X)
  n <- nrow(X)

  if (K==1)
  {
    if (method=='Semi-Online')
      centers <- ASGMedian(X)
    else
      centers <- WeiszfeldMedian(X)

    finalcenters <- centers
    finalcluster <- rep(1,n)
    finaldist <- sqrt(rowSums((X-matrix(rep(centers,n),nrow=n,byrow=TRUE))^2))
    finalSE <- sum(finaldist)

    return(list(cluster=finalcluster, centers=finalcenters, SE=finalSE))
  }

  centers <- matrix(0,nrow=K,ncol=d)
  finalcenters <- centers
  dist <- matrix(0,nrow=n,ncol=K)
  finaldist <- dist
  cluster <- rep(0,n)
  finalcluster <- cluster
  finalSE <- 10^10

  if (init==TRUE)
  {
    cluster <- genieclust::genie(X,k=K)
    for (k in 1:K)
    {
      I <- which(cluster==k)
      if (length(I)>1)
      {
        if (method=='Semi-Online')
          centers[k,] <- ASGMedian(X[I,])
        else
          centers[k,] <- WeiszfeldMedian(X[I,])
      }
      if (length(I) == 1)
      {
        centers[k,] <- X[I,]
      }
    }

    l <- 0
    distclust <- 1
    while (l < niter && distclust > 0)
    {
      cluster0 <- cluster
      l <- l+1
      for (k in 1:K)
      {
        dist[,k] <- sqrt(rowSums((X-matrix(rep(centers[k,],n),nrow=n,byrow=TRUE))^2))
      }
      cluster <- max.col(-dist)
      SE <- 0
      for (k in 1:K)
      {
        I <- which(cluster==k)
        SE <- SE+sum(dist[I,k])
        if (length(I)>1)
        {
          if (method=='Semi-Online')
            centers[k,] <- ASGMedian(X[I,])
          else
            centers[k,] <- WeiszfeldMedian(X[I,])
        }
        if (length(I) == 1)
        {
          centers[k,] <- X[I,]
        }
      }
      distclust <- sum((cluster-cluster0)^2)

      if (SE < finalSE)
      {
        finalcenters <- centers
        finalcluster <- cluster
        finaldist <- dist
        finalSE <- SE
      }
    }
  }

  if (ninit>0)
  {
    for (o in 1:ninit)
    {
      dist <- matrix(0,nrow=n,ncol=K)
      centers <- X[sample(1:n,K),]

      l <- 0
      distclust <- 1
      while (l < niter && distclust > 0)
      {
        cluster0 <- cluster
        l <- l+1
        for (k in 1:K)
        {
          dist[,k] <- sqrt(rowSums((X-matrix(rep(centers[k,],n),nrow=n,byrow=TRUE))^2))
        }
        cluster <- max.col(-dist)
        distclust <- sum((cluster-cluster0)^2)
        SE <- 0
        for (k in 1:K)
        {
          I <- which(cluster==k)
          SE <- SE+sum(dist[I,k])
          if (length(I)>1)
          {
            if (method=='Semi-Online')
              centers[k,] <- ASGMedian(X[I,])
            else
              centers[k,] <- WeiszfeldMedian(X[I,])
          }
          if (length(I) == 1)
          {
            centers[k,] <- X[I,]
          }
        }

        if (SE < finalSE)
        {
          finalcenters <- centers
          finalcluster <- cluster
          finaldist <- dist
          finalSE <- SE
        }
      }
    }
  }

  return(list(cluster=finalcluster,centers=finalcenters,dist=finaldist,SE=finalSE))
}

#' @keywords internal
kmedianClustering <- function (X, ncenters=2, gamma=1, alpha=0.75, nstart=10, nstartkmeans=10, iter.max=20)
{
  X <- as.matrix(X)
  p <- ncol(X)
  n <- nrow(X)

  ### Initialization
  if (is.matrix(ncenters)==FALSE) {
    k <- ncenters
    centers <- kmeans(X, ncenters,nstart=nstartkmeans,iter.max=iter.max,algorithm="MacQueen")$centers
  }
  else {
    k <- nrow(ncenters)
    centers <- ncenters
  }

  Z <- stoKmed_rcpp(X, X, centers, gamma=gamma, alpha=alpha)
  best <- sum(Z$wss)

  if (nstart >= 2) {
    for (i in 2:nstart) {
      ind.init <- sample(c(1:n),k)
      centers <- X[ind.init, ]
      newZ <-  stoKmed_rcpp(X, X, centers, gamma=gamma, alpha=alpha)
      if ((z <- sum(newZ$wss)) < best) {
        Z <- newZ
        best <- z
      }
    }
  }
  centers <- matrix(Z$centers, k)
  dimnames(centers) <- list(1L:k, dimnames(X)[[2L]])
  cluster <- Z$cl

  return(list(cluster=cluster, centers=centers, withinsrs=Z$wss, size=Z$nc))
}
