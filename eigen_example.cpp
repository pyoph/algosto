#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List eigen_cpp(arma::mat M) {
  arma::vec eigval;
  arma::mat eigvec;
  
  // Calcul des valeurs propres et vecteurs propres
  arma::eig_sym(eigval, eigvec, M);
  
  return Rcpp::List::create(Rcpp::Named("values") = eigval,
                            Rcpp::Named("vectors") = eigvec);
}