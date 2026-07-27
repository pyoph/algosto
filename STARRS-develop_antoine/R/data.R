#' Multivariate Gaussian data with Student contamination
#'
#' @description Synthetic data set
#'
#' A set of \eqn{n = 1500} observations with dimension \eqn{d = 5}, equally
#' spread into \eqn{K = 3}, each being sampled from a multivariate Gaussian distribution \eqn{\mathcal N _d (\mu_k , \Sigma_k )}
#' (with \eqn{k = 1, 2, 3}). The mean vectors are set to
#'
#' \eqn{\mu_1 = {0}_d , \mu_2 = 3 {1}_d , \mu_2 = −3 {1}_d}
#' (where \eqn{1_d} stands for the d-dimensional vector filled with ones) and the variance matrices to
#' \eqn{\Sigma_1 = 2 \; {I}_d, \qquad \Sigma_2 = \text{diag}([1, 2, \dots d]), \qquad \Sigma_3 = \text{diag}([1, 1/2, \dots 1/d])},
#' (\eqn{I_d} being the \eqn{d}-dimensional identity matrix).
#'
#' In addition, we consider that a fraction \eqn{f = 30\%} of the genuine observations have been replaced
#' with observations arising from a mixture with \eqn{K} components of student distribution,
#' with same respective location and scale parameters and with degree of freedom 1.
#'
#' @format Output of the function \code{GMMcontStudent}
#'
"GMMcontStudent"

#' Robust clustering output
#'
#' @description Output of robust clustering on synthetic data set GMMcontStudent
#'
#' @format Output of the call to the robust clustering function \code{multipleRobustMM(GMMcontStudent$X, nclust=1:5)}
#'
"robustGMMOutput"

#' Gaussian data with Student contamination
#'
#' @description TODO
#'
#' @format TODO
#'
"gaussStudentCont_n5000_d5"

#' Heteroscedastic setting gaussian GMM with Student contamination
#'
#' @description TODO
#'
#' @format TODO
#'
"heteroGMMStudentCont_K3_d5_n1500"

#' Homoscedastic setting gaussian GMM with Student contamination
#'
#' @description TODO
#'
#' @format TODO
#'
"homoGMMStudentCont_K3_d5_n1500"
