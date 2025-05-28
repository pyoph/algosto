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
 //Rcout << "Sigma2 = \n" << Sigma2 << "\n\n";
  for (int n = 0; n < 5; ++n) {
    arma::vec x = X.row(n).t();  // vecteur colonne
   
    mean   = mean   + (1.0 / (n + 1)) * (x - mean);
    Sigma2 = Sigma2  + (1.0 / (n + 1)) * (x * x.t() - Sigma2 );
    //Rcout << "Étape " << n+1 << " :\n";
    //Rcout << "x = " << x.t();
    //Rcout << "mean = " << mean.t();
    //Rcout << "Sigma2 = \n" << Sigma2 << "\n\n";
    //Rcout << "Sigma2 = \n" << Sigma2 << "\n\n";
    
  }
  
  return List::create(
    Named("mean") = mean,
    Named("Sigma2") = Sigma2
  );
}



// [[Rcpp::export]]
  
double mahalanobis_generalizedRcpp(const arma::rowvec& x,
                                   const arma::rowvec& moyennem,
                                   const arma::mat& eigvecs,
                                   const arma::rowvec& lambdaInit) {
  arma::rowvec centered = x - moyennem;
  arma::vec proj = eigvecs.t() * centered.t();  // Projeter x - moyennem sur les vecteurs propres
  arma::vec scaled = arma::square(proj) / lambdaInit.t();  // lambdaInit est un rowvec → transpose en vec
  return arma::accu(scaled);
}

// [[Rcpp::export]]


List update_mean_Sigma2Trimmed(const arma::mat& X, double cutoff) {
  int n_obs = X.n_rows;
  int d = X.n_cols;
  double dist ;
  //Rcout << "X.n_rows = " << n_obs << ", X.n_cols = " << d << "\n";
  
  arma::vec mean = arma::ones(d);       // moyenne empirique
  arma::vec meanOld = arma::ones(d);
  arma::mat Sigma2 = arma::zeros(d, d);  // matrice des moments d'ordre 2
  arma::mat Sigma2Old = arma::zeros(d, d);
  arma::mat eigvecs = arma::eye(d, d);
  arma::vec eigvals = arma::ones(d);
  
  // Rcout << "Sigma2 = \n" << Sigma2 << "\n\n";
  //Rcout << "Sigma2 = \n" << Sigma2 << "\n\n";
  for (int n = 0; n < 5; ++n) {
    arma::vec x = X.row(n).t();  // vecteur colonne
    
    mean   = mean   + (1.0 / (n + 1)) * (x - mean);
    Sigma2 = Sigma2  + (1.0 / (n + 1)) * (x * x.t() - Sigma2 );
    //Rcout << "Étape " << n+1 << " :\n";
    //Rcout << "x = " << x.t();
    //Rcout << "mean = " << mean.t();
    //Rcout << "Sigma2 = \n" << Sigma2 << "\n\n";
    //Rcout << "Sigma2 = \n" << Sigma2 << "\n\n";
    
     arma::eig_sym(eigvals,eigvecs,Sigma2);
    //dist = mahalanobis_generalizedRcpp(x,mean,eigvecs,eigvals);
    //if (dist > cutoff){Sigma2= Sigma2Old;
      //mean = meanOld;}
    //else {Sigma2Old = Sigma2;
      //meanOld = mean;}
      
      Rcout << "eigvals = " << eigvals;
      
      Rcout << "eigvecs = " << eigvecs;
  }
  
  return List::create(
    Named("mean") = mean,
    Named("Sigma2") = Sigma2
  );
}