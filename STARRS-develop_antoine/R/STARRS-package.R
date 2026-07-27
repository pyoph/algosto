#' @useDynLib STARRS, .registration=TRUE
#'
#' @importFrom Rcpp sourceCpp evalCpp
#' @importFrom stats BIC kmeans
#' @importFrom DescTools  RobScale
#' @importFrom mvtnorm rmvnorm dmvnorm
#' @importFrom parallel detectCores
#' @importFrom foreach foreach "%do%"
#' @importFrom doFuture "%dofuture%"
#' @importFrom stats rnorm rt runif rchisq
#' @importFrom ggplot2 ggplot scale_x_continuous scale_y_continuous geom_point geom_line aes scale_color_manual facet_wrap geom_bar scale_size_continuous
#' @importFrom reshape2 melt
#' @importFrom capushe capushe
#' @importFrom knitr knit
#' @importFrom bookdown render_book
#'
NULL
#> NULL
