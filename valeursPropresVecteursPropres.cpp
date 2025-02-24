#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Fonction pour calculer les valeurs propres et les vecteurs propres
// [[Rcpp::export]]
Rcpp::List calculValeursEtVecteursPropresArmadillo(const arma::mat& mat) {
  // Vecteur pour les valeurs propres (complexes)
  arma::cx_vec valeurs_propres;
  
  // Matrice pour les vecteurs propres (complexes)
  arma::cx_mat vecteurs_propres;
  
  // Calcul des valeurs propres et des vecteurs propres
  // Utilisation correcte de eig_gen
  if (!arma::eig_gen(valeurs_propres, vecteurs_propres, mat)) {
    Rcpp::stop("Le calcul des valeurs propres a échoué.");
  }
  
  // Retourner les résultats sous forme de liste R
  return Rcpp::List::create(
    Rcpp::Named("valeurs_propres") = valeurs_propres,
    Rcpp::Named("vecteurs_propres") = vecteurs_propres
  );
}