#' @title Robust PCA using Weiszfeld algorithm
#'
#' @description Robust PCA using Weiszfeld algorithm for the median and the median covariation computation
#'
#' @param X Matrix of size N*p corresponding to the data. The rows are the observations.
#' @param init_median Initial value of the median, by default equal to 0
#' @param init_median_cov Initial value of the median covariation matrix, by default equal to the identity matrix
#' @param weights Vector of size N of the weights, by default equal to 1
#' @param scores Number of scores to compute, by default 2
#' @param epsilon Stopping criterion: the algorithm stops when the difference between two iterations is less than epsilon, by default 1e-08
#' @param nitermax Maximum number of iterations for the Weiszfeld algorithm, by default 100
#'
#' @return A list containing :
#' - `median` : the estimated median
#' - `covmedian` : the estimated median covariation matrix
#' - `scores` : the scores
#' - `vectors`: the eigen vectors
#'
#' @md
#'
#' @examples
#' N <- 1000
#' p <- 4
#' X <- matrix(rnorm(N*p),ncol=p)
#' WeiszfeldRobustPCA(X)
#'
#' @export
WeiszfeldRobustPCA <- function(X, init_median=rep(0,ncol(X)), init_median_cov=(diag(ncol(X))), weights=rep(1,nrow(X)), scores=2, epsilon=1e-08, nitermax=100)
{
  med <- WeiszfeldMedian(X, init_median=init_median, weights=weights, epsilon=epsilon, nitermax=nitermax)
  medcovmat <- WeiszfeldMedianCovariance(X, init_median_cov=init_median_cov, weights=weights, median_est=med, epsilon=epsilon, nitermax=nitermax)

  vectors <- RSpectra::eigs_sym(medcovmat, scores)$vectors
  vscores <- sweep(X,2,med)%*%vectors

  return(list(median=med, covmedian=medcovmat, scores=vscores, vectors=vectors))
}

#' @title Robust PCA using ASG algorithm
#'
#' @description Robust PCA using ASG algorithm for the median and the median covariation computation
#'
#' @param X Matrix of size N*p corresponding to the data. The rows are the observations.
#' @param init_median Initial value of the median, by default equal to 0
#' @param init_median_cov Initial value of the median covariation matrix, by default equal to the identity matrix
#' @param weights Vector of size N of the weights, by default equal to 1
#' @param epsilon Stopping criterion: the algorithm stops when the difference between two iterations is less than epsilon, by default 1e-08
#' @param scores Number of scores to compute, by default 2
#' @param gamma ASG parameter (2 by default)
#' @param alpha ASG parameter (0.75 by default)
#' @param nstart Number of starts for the ASG algorithm, by default 1
#'
#' @return A list containing :
#'
#' - `median` : the estimated median
#' - `covmedian` : the estimated median covariation matrix
#' - `scores` : the scores
#' - `vectors`: the eigen vectors
#'
#' @md
#'
#' @examples
#' N <- 1000
#' p <- 4
#' X <- matrix(rnorm(N*p),ncol=p)
#' ASGRobustPCA(X)
#'
#' @export
ASGRobustPCA <- function(X, init_median=rep(0,ncol(X)), init_median_cov=diag(ncol(X)), weights=rep(1,nrow(X)), epsilon=1e-08, scores=2, gamma=2, alpha=0.75, nstart=1)
{
  med <- ASGMedian(X, init_median=init_median, weights=weights, gamma=gamma, alpha=alpha, nstart=nstart, epsilon=epsilon)
  medcovmat <- ASGMedianCovariance(X, init_median_cov=init_median_cov, median_est=med, weights=weights, gamma=gamma, alpha=alpha, nstart=nstart)

  vectors <- RSpectra::eigs_sym(medcovmat, scores)$vectors
  vscores <- sweep(X,2,med)%*%vectors

  return(list(median=med,covmedian=medcovmat,scores=scores,vectors=vectors))
}
