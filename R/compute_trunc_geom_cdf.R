#' Compute the Truncated Geometric Cumulative Distribution Function (CDF)
#'
#' This function calculates the cumulative distribution function (CDF) for a truncated geometric distribution with a maximum of \code{b} support points.
#' The geometric distribution is truncated to values from 0 to \code{b - 1}, then normalized.
#'
#' @param lambda A numeric value between 0 and 1 representing the success probability of the geometric distribution.
#' @param b An integer specifying the number of support points for truncation. The distribution is truncated at \code{b - 1}.
#'
#' @return A numeric vector of length \code{b} containing the cumulative probabilities of the truncated and normalized geometric distribution.
#'
#' @examples
#' compute_trunc_geom_cdf(lambda = 0.2, b = 10)
#'
#' @export
compute_trunc_geom_cdf <- function(lambda, b) {
  probs <- (1 - lambda)^(0:(b - 1)) * lambda
  probs <- probs / sum(probs)
  cumsum(probs)
}
