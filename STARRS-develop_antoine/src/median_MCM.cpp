#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' @title Weiszfeld geometric median computation
//'
//' @description Compute the median of a matrix X using the Weiszfeld algorithm
//'
//' \eqn{m_{t+1} = \frac{\sum_{k=1}{^n}X_k/\| X_k-m_t\|}{\sum_{k=1}{^n}1/\| X_k-m_t\|}}
//'
//' @param X A numerical matrix corresponding to the data of size \eqn{n\times p}. The rows are \eqn{n} observations in \eqn{\mathbf R^p}.
//' @param init_median Initial value of the median, by default equal to 0 (internally handled)
//' @param weights Vector of size \eqn{n} of the weights, by default equal to 1 (internally handled)
//' @param epsilon Stopping criterion: the algorithm stops when the difference between two iterations is less than epsilon, by default 1e-08
//' @param nitermax Maximum number of iterations for the Weiszfeld algorithm, by default 100
//'
//' @return The estimated median
//'
//' @md
//'
//' @examples
//' N <- 1000
//' p <- 4
//' X <- matrix(rnorm(N*p),ncol=p)
//' WeiszfeldMedian(X)
//'
//' @export
// [[Rcpp::export]]
arma::rowvec WeiszfeldMedian(const arma::mat& X,
                             const Rcpp::Nullable<arma::rowvec>& init_median = R_NilValue,
                             const Rcpp::Nullable<arma::rowvec>& weights = R_NilValue,
                             double epsilon = 1e-08,
                             int nitermax = 100){
  // Inputs
  const int n = X.n_rows ;
  const int p = X.n_cols ;

  // Containers and initialisation of the algorithm
  arma::rowvec meanvec(p);
  arma::rowvec medvec(p);
  arma::rowvec W(n);
  arma::rowvec w(n);

  // Set default values if needed
  if(init_median.isNotNull()){
    meanvec = Rcpp::as<arma::rowvec>(init_median);
    medvec = Rcpp::as<arma::rowvec>(init_median);
  }
  else{
    meanvec = arma::zeros<arma::rowvec>(p);
    medvec = arma::zeros<arma::rowvec>(p);
  }

  if(weights.isNotNull()){
    W = Rcpp::as<arma::rowvec>(weights);
  }
  else{
    W = arma::ones<arma::rowvec>(n);
  }

  double diffxn, normxm = 10;
  int iter = 0;

  /* Weiszfeld loop */
  while (iter < nitermax and normxm > epsilon)
  {
    for (int it=0 ; it < n ; it++){
      diffxn = norm(X.row(it)-meanvec);
      w(it) = (diffxn > 0) ? W(it)/diffxn : 0;
    }
    w = w/sum(w); // normalisation
    medvec = w*X; // geometric median update, element by element product
    normxm = norm(medvec-meanvec)/sqrt(double(p)); // stopping criterion
    meanvec = medvec;
    iter++;
  }
  // Returns ;
  return medvec;
}

//' @title Weiszfeld geometric median covariation matrix computation
//'
//' @description Compute the median covariation matrix of a matrix X using the Weiszfeld algorithm
//'
//' TODO comparer avec notations papier --> ajouter les poids w_k comme ds med
//'
//' \eqn{V_{t+1} = \frac{\sum_{k=1}{^n}\| (X_k-m^*)(X_k-m^*)^T-V_t\|_F^{-1}(X_k-m^*)(X_k-m^*)^T}{\sum_{k=1}{^n}\| (X_k-m^*)(X_k-m^*)^T-V_t\|_F^{-1}}}
//'
//' @param X A numerical matrix corresponding to the data of size \eqn{n\times p}. The rows are \eqn{n} observations in \eqn{\mathbf R^p}.
//' @param median_est Estimation of the median \eqn{m^*}, by default \code{WeiszfeldMedian(X)} (internally handled)
//' @param init_median_cov Initial value of the median covariation matrix, by default equal to the identity matrix (internally handled)
//' @param weights Vector of size $n$ of the weights, by default equal to 1 (internally handled)
//' @param epsilon Stopping criterion: the algorithm stops when the difference between two iterations is less than epsilon, by default 1e-08
//' @param nitermax Maximum number of iterations for the Weiszfeld algorithm, by default 100
//'
//'
//' @return The estimated median covariation matrix
//'
//' @examples
//' N <- 1000
//' p <- 4
//' X <- matrix(rnorm(N*p),ncol=p)
//' med_est <- WeiszfeldMedian(X)
//' WeiszfeldMedianCovariance(X,median_est=med_est)
//'
//' @export
// [[Rcpp::export]]
arma::mat WeiszfeldMedianCovariance(const arma::mat& X,
                                    const Rcpp::Nullable<arma::rowvec>& median_est = R_NilValue,
                                    const Rcpp::Nullable<arma::mat>& init_median_cov = R_NilValue,
                                    const Rcpp::Nullable<arma::rowvec>& weights = R_NilValue,
                                    double epsilon=1e-08,
                                    int nitermax=100)
{
  // X : n * p  matrix
  // Inputs
  const int n = X.n_rows ;
  const int p = X.n_cols ;

  arma::rowvec median;
  arma::mat median_cov;
  arma::rowvec W;

  // Set default values if needed
  if(median_est.isNotNull()){
    median = Rcpp::as<arma::rowvec>(median_est);
  }
  else{
    median = WeiszfeldMedian(X);
  }

  if(init_median_cov.isNotNull()){
    median_cov = Rcpp::as<arma::mat>(init_median_cov);
  }
  else{
    median_cov = arma::eye(p,p);
  }

  if(weights.isNotNull()){
    W = Rcpp::as<arma::rowvec>(weights);
  }
  else{
    W = arma::ones<arma::rowvec>(n);
  }

  // Containers
  arma::mat Xcent(n,p);
  for (int it=0 ; it < n ; it++)
  {
    Xcent.row(it) = X.row(it)-median;
  }

  arma::mat medcovest(p,p);

  arma::rowvec w(n);
  double diffxn, normxm = 1;
  int iter = 0;

  while (iter < nitermax and normxm > epsilon )
  {
    medcovest.fill(0.0);
    for (int it=0 ; it < n ; it++){
      diffxn = arma::norm(arma::trans(Xcent.row(it))*Xcent.row(it)-median_cov,"fro");
      w(it) = (diffxn > 0) ? W(it)/diffxn : 0;
    }
    w = w/sum(w);
    for (int it=0 ; it < n ; it++)
    {
      medcovest += w(it)*arma::trans(Xcent.row(it))*Xcent.row(it);
    }
    normxm = norm(medcovest-median_cov,"fro")/p;
    median_cov = medcovest;
    iter++;
  }
  // Returns ;
  return medcovest;
}


//' @title Averaged Stochastic Gradient geometric median computation
//'
//' @description Compute the median of a matrix X using the Averaged Stochastic Gradient algorithm
//'
//' \eqn{m_{k+1} = m_k + \gamma_{k+1}\frac{X_{k+1} − m_k}{\| X_{k+1} − m_k\|}}
//'
//' \eqn{\bar m_{k+1} = \bar m_k + \frac 1 {k+1} (m_{k+1} − \bar m_k )}
//'
//' @param X  A numerical matrix corresponding to the data of size \eqn{n\times p}. The rows are n observations in \eqn{\mathbf R^p}.
//' @param init_median  Initial value of the median, by default equal to 0 (internally handled)
//' @param weights Vector of size n of the weights, by default equal to 1 (internally handled)
//' @param gamma ASG parameter (2 by default)
//' @param alpha ASG parameter (0.75 by default)
//' @param nstart Number of starts for the ASG algorithm, by default 1
//' @param epsilon Stopping criterion: the algorithm stops when the difference between two iterations is less than epsilon, by default 1e-08
//'
//' @return The estimated median
//'
//' @examples
//' N <- 1000
//' p <- 4
//' X <- matrix(rnorm(N*p),ncol=p)
//' ASGMedian(X)
//'
//' @export
// [[Rcpp::export]]
arma::rowvec ASGMedian(const arma::mat& X,
                        const Rcpp::Nullable<arma::rowvec>& init_median = R_NilValue,
                        const Rcpp::Nullable<arma::rowvec>& weights = R_NilValue,
                        double gamma = 2,
                        double alpha = 0.75,
                        int nstart = 1,
                        double epsilon = 1e-8)
{
  try {
    if (alpha < 0.5 || alpha > 1) {
      throw std::range_error("alpha must be in [0.5, 1]");
    }
    if (gamma < 0 ) {
      throw std::range_error("gamma must be positive");
    }
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }

  // X : n * p  matrix
  // Inputs
  const int n = X.n_rows;
  const int p = X.n_cols;

  arma::rowvec medvec;
  arma::rowvec medrm;
  arma::rowvec W;

  // Set default values if needed
  if(init_median.isNotNull()){
    medvec = Rcpp::as<arma::rowvec>(init_median);
    medrm = Rcpp::as<arma::rowvec>(init_median);
  }
  else{
    medvec = arma::zeros<arma::rowvec>(p);
    medrm = arma::zeros<arma::rowvec>(p);
  }

  if(weights.isNotNull()){
    W = Rcpp::as<arma::rowvec>(weights);
  }
  else{
    W = arma::ones<arma::rowvec>(n);
  }

  // Containers and initialisation of the algorithm
  double w, normxm ;

  // Number of replications of the algorithm
  for (int nbcomp = 0 ; nbcomp < nstart ; nbcomp++){
    // Stochastic gradient algorithms
    for (int it = 1 ; it < n ; it++)
    {
      normxm = arma::norm(X.row(it)-medrm);
      if (normxm > epsilon) {
        w = W(it)*sqrt(double(p))*gamma*pow(double(it+1),-alpha)/normxm;
        medrm += w * (X.row(it)-medrm) ;
      }
      medvec += (medrm-medvec)/(it+1);
    }
  }
  // Returns ;
  return medvec;
}

//' @title Averaged Stochastic Gradient geometric median covariation matrix computation
//'
//' @description Compute the median covariation matrix of a matrix X using the Averaged Stochastic Gradient algorithm
//'
//' TODO -V_n ??
//' \eqn{V_{k+1} = V_k + \gamma_{k+1}\frac{(X_{k+1} − \bar m_k)(X_{k+1} − \bar m_k)^T-V_n }{\| (X_{k+1} − \bar m_k)(X_{k+1} − \bar m_k)^T-V_n \|_F}}
//'
//' \eqn{\bar V_{k+1} = \bar V_k + \frac 1 {k+1} (V_{k+1} − \bar V_k )}
//'
//' @param X  A numerical matrix corresponding to the data of size \eqn{n\times p}. The rows are n observations in \eqn{\mathbf R^p}.
//' @param median_est Estimation of the median \eqn{m^*}, typically an output of \code{ASGMedian(X)} (internally handled)
//' @param init_median_cov Initial value of the median covariation matrix, by default equal to the identity matrix (internally handled)
//' @param weights Vector of size $n$ of the weights, by default equal to 1 (internally handled)
//' @param gamma ASG parameter (2 by default)
//' @param alpha ASG parameter (0.75 by default)
//' @param nstart Number of starts for the ASG algorithm, by default 1
//'
//' @return The estimated median covariation matrix
//'
//' @examples
//' N <- 1000
//' p <- 4
//' X <- matrix(rnorm(N*p),ncol=p)
//' med <- ASGMedian(X)
//' ASGMedianCovariance(X, median_est=med)
//'
//' @export
// [[Rcpp::export]]
arma::mat ASGMedianCovariance(const arma::mat&  X,
                              const Rcpp::Nullable<arma::rowvec>& median_est = R_NilValue,
                              const Rcpp::Nullable<arma::mat>& init_median_cov = R_NilValue,
                              const Rcpp::Nullable<arma::rowvec>& weights = R_NilValue,
                              double gamma = 2,
                              double alpha = 0.75,
                              int nstart = 1){
  // X : n * p  matrix
  // Inputs
  const int n = X.n_rows ;
  const int p = X.n_cols ;

  arma::rowvec median;
  arma::mat median_cov;
  arma::rowvec W;

  // Set default values if needed
  if(median_est.isNotNull()){
    median = Rcpp::as<arma::rowvec>(median_est);
  }
  else{
    median = ASGMedian(X);
  }

  if(init_median_cov.isNotNull()){
    median_cov = Rcpp::as<arma::mat>(init_median_cov);
  }
  else{
    median_cov = arma::eye(p,p);
  }

  if(weights.isNotNull()){
    W = Rcpp::as<arma::rowvec>(weights);
  }
  else{
    W = arma::ones<arma::rowvec>(n);
  }

  // Containers
  arma::mat  diffmat(p,p), medrm(p,p),diffmed(p,p)  ;
  double nrmrm, w ;

  // Initialization of the algorithm
  medrm = median_cov ;

  // Number of replications of the algorithm
  for (int nbcomp = 0 ; nbcomp < nstart ; nbcomp++){
    // Stochastic gradient algorithms
    for (int it = 1 ; it < n ; it++)
    {
      diffmed = X.row(it)-median;
      diffmat = arma::trans(diffmed)*diffmed;
      diffmat -= medrm ;
      nrmrm = arma::norm(diffmat,"fro") ; // Frobenius norm divided by the dimension
      w = p *W(it) * gamma * pow(double(it+1),-alpha) * pow(nrmrm,-1);
      medrm += w*diffmat ;
      median_cov += (medrm-median_cov)/(it+1);
    }
  }
  return median_cov;
}
