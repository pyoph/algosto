#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List update_mean_cov(const arma::vec& Xn1, 
                     const arma::vec& mean_n, 
                     const arma::mat& cov_n, 
                     int n) {
  // Mise à jour de la moyenne
  arma::vec mean_n1 = mean_n + (1.0 / (n + 1)) * (Xn1 - mean_n);
  
  // Mise à jour de la covariance non centrée (matrice des moments)
  arma::mat cov_n1 = cov_n + (1.0 / (n + 1)) * (Xn1 * Xn1.t() - cov_n);
  
  return List::create(
    Named("mean") = mean_n1,
    Named("cov") = cov_n1
  );
}
