#' Generate a Random Positive Definite Correlation Matrix
#'
#' This function creates a random symmetric correlation matrix of dimension \code{p} × \code{p},
#' where off-diagonal elements are drawn uniformly at random between \code{lower} and \code{upper}, and the diagonal is set to 1.
#' The resulting matrix is adjusted to be the nearest positive definite matrix using Higham's algorithm.
#'
#' @param p An integer specifying the number of variables (i.e., the dimension of the correlation matrix).
#' @param lower A numeric value in (0, 1) giving the lower bound for the uniformly sampled off-diagonal correlations (default is 0.1).
#' @param upper A numeric value in (0, 1) giving the upper bound for the uniformly sampled off-diagonal correlations (default is 0.9).
#'
#' @return A \code{p} × \code{p} positive definite correlation matrix with unit diagonal and off-diagonal values between \code{lower} and \code{upper}.
#'
#' @details The function uses the \code{Matrix::nearPD()} function to ensure the matrix is valid as a correlation matrix.
#'
#' @examples
#' generate_positive_correlation_matrix(p = 5)
#'
#' @importFrom Matrix nearPD
#' @export
generate_positive_correlation_matrix <- function(p, lower = 0.1, upper = 0.9) {
  R <- matrix(runif(p * p, lower, upper), nrow = p)
  diag(R) <- 1
  R[lower.tri(R)] <- t(R)[lower.tri(R)]
  R <- as.matrix(nearPD(R, corr = TRUE)$mat)
  return(R)
}
