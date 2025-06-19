#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


 
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


