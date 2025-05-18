#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List update_mean_Sigma2(const arma::mat& X) {
  int n_obs = X.n_rows;
  int d = X.n_cols;
  //Rcout << "X.n_rows = " << n_obs << ", X.n_cols = " << d << "\n";
  
  arma::vec mean = arma::zeros(d);       // moyenne empirique
  arma::mat Sigma2 = arma::zeros(d, d);  // matrice des moments d'ordre 2
 // Rcout << "Sigma2 = \n" << Sigma2 << "\n\n";
 Rcout << "Sigma2 = \n" << Sigma2 << "\n\n";
  for (int n = 0; n < 5; ++n) {
    arma::vec x = X.row(n).t();  // vecteur colonne
   
    mean   = mean   + (1.0 / (n + 1)) * (x - mean);
    Sigma2 = Sigma2  + (1.0 / (n + 1)) * (x * x.t() - Sigma2 );
    //Rcout << "Étape " << n+1 << " :\n";
    //Rcout << "x = " << x.t();
    //Rcout << "mean = " << mean.t();
    //Rcout << "Sigma2 = \n" << Sigma2 << "\n\n";
    Rcout << "Sigma2 = \n" << Sigma2 << "\n\n";
    
  }
  
  return List::create(
    Named("mean") = mean,
    Named("Sigma2") = Sigma2
  );
}