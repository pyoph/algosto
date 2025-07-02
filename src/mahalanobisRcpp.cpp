#include <RcppArmadillo.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

double mahalanobis_generalizedRcpp(const arma::rowvec& x,
                                   const arma::rowvec& moyennem,
                                   const arma::mat& eigvecs,
                                   const arma::rowvec& lambdaInit) {
  arma::rowvec centered = x - moyennem;
  arma::vec proj = eigvecs.t() * centered.t();  // Projeter x - moyennem sur les vecteurs propres
  arma::vec scaled = arma::square(proj) / lambdaInit.t();  // lambdaInit est un rowvec → transpose en vec
  return arma::accu(scaled);
}




// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
*/
