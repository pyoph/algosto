#include <Rcpp.h>
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

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
*/
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

// Fonction pour calculer les valeurs propres et les vecteurs propres
// [[Rcpp::export]]
Rcpp::List calculValeursEtVecteursPropres(const Eigen::MatrixXd& mat) {
  // Créer un solveur pour la matrice
  Eigen::EigenSolver<Eigen::MatrixXd> solver(mat);
  
  // Extraire les valeurs propres (sous forme de vecteur complexe)
  Eigen::VectorXcd valeurs_propres = solver.eigenvalues();
  
  // Extraire les vecteurs propres (sous forme de matrice complexe)
  Eigen::MatrixXcd vecteurs_propres = solver.eigenvectors();
  
  // Retourner les résultats sous forme de liste R
  return Rcpp::List::create(
    Rcpp::Named("valeurs_propres") = valeurs_propres,
    Rcpp::Named("vecteurs_propres") = vecteurs_propres
  );
}