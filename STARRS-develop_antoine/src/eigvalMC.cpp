#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' @title Robbins-Monro method for the estimation of the eigen values of the variance-covariance matrix
//'
//' @description Given the eigen values \eqn{\delta} of the median covariation matrix, this function estimates
//' the eigen values \eqn{\lambda} of a variance-covariance matrix using the Robbins-Monro method and Monte Carlo approximation.
//'
//' @param U Matrix of size N*p corresponding to \eqn{\Sigma^{-1/2}(X-\mu)}. The rows are the observations.
//'
//' - In a gaussian model typically `U = matrix(rnorm(N*p),ncol=p)`
//' - In a Student model `U <- matrix(rnorm(N*p)/sqrt(rchisq(1,df=df))*sqrt((df-2)), ncol=p))`
//' - In a Laplace model `U <- LaplacesDemon::rmvl(N,mu=rep(0,p),Sigma=diag(p))`
//' @param delta Vector of size p of the eigen values of the median covariation matrix
//' @param init Initial value of the vector of eigen values of the variance-covariance matrix, by default equal to delta (internally handled)
//' @param gamma Robbins-Monro parameter (0.75 by default)
//' @param c Robbins-Monro parameter (2 by default)
//' @param w Robbins-Monro parameter (2 by default)
//' @param c_bar Hyperparameter for fine tuning. Do not change unless you are an experienced user (0 by default)
//' @param c_tilde Hyperparameter for fine tuning. Do not change unless you are an experienced user (0 by default)
//' @param sumlog Hyperparameter for fine tuning. Do not change unless you are an experienced user (1 by default)
//' @param epsilon Stopping criterion: the algorithm stops when the difference between two iterations is less than epsilon, by default 1e-08
//'
//' @return The vector of the estimated eigen values \eqn{\lambda}
//'
//' @examples
//' delta <- c(1, 1)
//' N <- 1000
//' p <- length(delta)
//' U <- matrix(rnorm(N*p),ncol=p)
//' robbinsMC(U=U,delta=delta)
//'
//' @export
// [[Rcpp::export]]
arma::vec robbinsMC(const arma::mat &U,
                    const arma::vec &delta,
                    Rcpp::Nullable<arma::vec> init = R_NilValue,
                    double gamma = 0.75,
                    double c = 2,
                    double w = 2,
                    double c_bar = 0,
                    double c_tilde = 0,
                    double sumlog = 1,
                    double epsilon = 1e-8)
{
  int N = U.n_rows;

  arma::vec lambda;
  arma::vec lambda_bar;

  // Set default values if needed
  if (init.isNotNull()){
    lambda = Rcpp::as<arma::vec>(init);
    lambda_bar = Rcpp::as<arma::vec>(init);
  }
  else {
    lambda = delta;
    lambda_bar = delta;
  }

  for (int i = 0; i < N; ++i)
  {
    arma::rowvec U_i = U.row(i);
    arma::vec U_i_sq = trans(U_i % U_i); // U_i^2

    arma::vec lambda_Ui2 = lambda % U_i_sq;
    arma::vec delta_minus_lambda_Ui2 = delta - lambda_Ui2;

    double E2 = std::pow(
        arma::dot(delta_minus_lambda_Ui2, delta_minus_lambda_Ui2) + arma::as_scalar(lambda_Ui2.t() * lambda_Ui2) - arma::accu(lambda_Ui2 % lambda_Ui2),
        -0.5);

    arma::vec E1 = U_i_sq * E2;

    arma::vec lambda_old = lambda;

    double i_gamma = std::pow(i + 1 + c_tilde, -gamma);
    lambda = lambda - c * i_gamma * lambda % E1 + c * i_gamma * delta * E2;

    // Weighted averaged version of lambda
    sumlog = sumlog + std::pow(std::log(i + 2 + c_bar), w);
    lambda_bar = lambda_bar + std::pow(std::log(i + 2 + c_bar), w) / sumlog * (lambda - lambda_bar);

    double eps = std::sqrt(arma::accu(arma::square(lambda - lambda_old)));

    if (eps < epsilon)
      break;
  }

  return lambda_bar;
}

//' @title Fix point method for the estimation of the eigen values of the variance-covariance matrix
//'
//' @description Given the eigen values \eqn{\delta} of the median covariation matrix, this function estimates
//' the eigen values \eqn{\lambda} of a variance-covariance matrix using the fix point method and Monte Carlo approximation.
//'
//' @param U Matrix of size N*p corresponding to \eqn{\Sigma^{-1/2}(X-\mu)}. The rows are the observations.
//'
//' - In a gaussian model typically `U = matrix(rnorm(N*p),ncol=p)`
//' - In a Student model `U <- matrix(rnorm(N*p)/sqrt(rchisq(1,df=df))*sqrt((df-2)), ncol=p))`
//' - In a Laplace model `U <- LaplacesDemon::rmvl(N,mu=rep(0,p),Sigma=diag(p))`
//' @param delta Vector of size p of the eigen values of the median covariation matrix
//' @param init Initial value of the vector of eigen values of the variance-covariance matrix, by default equal to delta (internally handled)
//' @param niter Maximum number of iterations for the gradient descent algorithm, by default 10
//' @param epsilon Stopping criterion: the algorithm stops when the difference between two iterations is less than epsilon, by default 1e-08
//' @param out_index Indexes of the iterations for which we want to output the values of the estimates, default is the last iteration
//'
//' @return The vector of the estimated eigen values \eqn{\lambda}
//'
//' @examples
//' delta <- c(1, 1)
//' N <- 1000
//' p <- length(delta)
//' U <- matrix(rnorm(N*p),ncol=p)
//' gradMC(U=U,delta=delta)
//'
//' @export
// [[Rcpp::export]]
arma::vec fixMC(const arma::mat &U,
                const arma::vec &delta,
                Rcpp::Nullable<arma::vec> init = R_NilValue,
                int niter = 10,
                double epsilon = 1e-8,
                IntegerVector out_index = IntegerVector::create())
{
  int N = U.n_rows;
  int p = U.n_cols;

  arma::vec lambda;

  // Set default values if needed
  if (init.isNotNull())
    lambda = Rcpp::as<arma::vec>(init);
  else
    lambda = delta;

  for (int l = 0; l < niter; ++l)
  {
    arma::vec E1(p, arma::fill::zeros);
    double E2 = 0.0;

    for (int i = 0; i < N; ++i)
    {
      arma::rowvec U_i = U.row(i);
      arma::vec U_i_sq = trans(U_i % U_i); // element-wise square

      arma::vec lambda_Ui2 = lambda % U_i_sq;
      arma::vec delta_minus_lambda_Ui2 = delta - lambda_Ui2;

      double term1 = arma::dot(delta_minus_lambda_Ui2, delta_minus_lambda_Ui2);
      double term2 = arma::as_scalar(lambda_Ui2.t() * lambda_Ui2);
      double term3 = arma::accu(lambda_Ui2 % lambda_Ui2);

      double Ei = std::pow(term1 + term2 - term3, -0.5);

      E1 += U_i_sq * Ei;
      E2 += Ei;
    }

    arma::vec lambda_old = lambda;
    lambda = E2 * delta / E1;

    double eps = std::sqrt(arma::accu(arma::square(lambda - lambda_old)));

    if (eps < epsilon)
      break;
  }

  return lambda;
}

//' @title Gradient descent for the estimation of the eigen values of the variance-covariance matrix
//'
//' @description Given the eigen values \eqn{\delta} of the median covariation matrix, this function estimates
//' the eigen values \eqn{\lambda} of a variance-covariance matrix using the gradient descent method and Monte Carlo approximation.
//'
//' @param U Matrix of size N*p corresponding to \eqn{\Sigma^{-1/2}(X-\mu)}. The rows are the observations.
//'
//' - In a gaussian model typically `U = matrix(rnorm(N*p),ncol=p)`
//' - In a Student model `U <- matrix(rnorm(N*p)/sqrt(rchisq(1,df=df))*sqrt((df-2)), ncol=p))`
//' - In a Laplace model `U <- LaplacesDemon::rmvl(N,mu=rep(0,p),Sigma=diag(p))`
//' @param delta Vector of size p of the eigen values of the median covariation matrix
//' @param init Initial value of the vector of eigen values of the variance-covariance matrix, by default equal to delta (internally handled)
//' @param niter Maximum number of iterations for the gradient descent algorithm, by default 10
//' @param epsilon Stopping criterion: the algorithm stops when the difference between two iterations is less than epsilon, by default 1e-08
//' @param step Step of the gradient descent, by default 1
//' @param out_index Indexes of the iterations for which we want to output the values of the estimates, default is the last iteration
//'
//' @return A list containing the estimated eigen values \eqn{\lambda}
//' and the values of the estimates for different iterations
//'
//' @examples
//' delta <- c(1, 1)
//' N <- 1000
//' p <- length(delta)
//' U <- matrix(rnorm(N*p),ncol=p)
//' gradMC(U=U,delta=delta)
//'
//' @export
// [[Rcpp::export]]
arma::vec gradMC(const arma::mat &U,
                 const arma::vec &delta,
                 Rcpp::Nullable<arma::vec> init = R_NilValue,
                 int niter = 10,
                 double epsilon = 1e-8,
                 NumericVector step = NumericVector::create(),
                 IntegerVector out_index = IntegerVector::create())
{
  int N = U.n_rows;
  int p = U.n_cols;

  arma::vec lambda;

  // Set default values if needed
  if (init.isNotNull())
    lambda = Rcpp::as<arma::vec>(init);
  else
    lambda = delta;

  arma::mat lambdalist(p, 0); // Each column will be a lambda

  // Default steps if none provided
  if (step.size() == 0)
  {
    step = rep(1.0, niter);
  }

  for (int l = 0; l < niter; ++l)
  {
    arma::vec E1(p, arma::fill::zeros);
    double E2 = 0.0;

    for (int i = 0; i < N; ++i)
    {
      arma::rowvec U_i = U.row(i);
      arma::vec U_i_sq = trans(U_i % U_i); // U_i^2

      arma::vec lambda_Ui2 = lambda % U_i_sq;
      arma::vec delta_minus_lambda_Ui2 = delta - lambda_Ui2;

      double term1 = arma::dot(delta_minus_lambda_Ui2, delta_minus_lambda_Ui2);
      double term2 = arma::as_scalar(lambda_Ui2.t() * lambda_Ui2);
      double term3 = arma::accu(lambda_Ui2 % lambda_Ui2);

      double Ei = std::pow(term1 + term2 - term3, -0.5);

      E1 += U_i_sq * Ei;
      E2 += Ei;
    }

    arma::vec lambda_old = lambda;

    double step_l = step[l];
    lambda = lambda - (step_l / N) * lambda % E1 + (step_l / N) * delta * E2;

    double eps = std::sqrt(arma::accu(arma::square(lambda - lambda_old)));

    // Save if iteration is in out_index
    if (std::find(out_index.begin(), out_index.end(), l + 1) != out_index.end())
    {
      lambdalist.insert_cols(lambdalist.n_cols, lambda);
    }

    // Optional stopping criterion:
    if (eps < epsilon)
      break;
  }

  return lambda;
}
