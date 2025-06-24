#' Map Uniform Samples to Truncated Geometric Categories (Vectorized)
#'
#' This function maps uniform random values to ordinal categories based on a truncated geometric cumulative distribution function (CDF).
#' It is useful for simulating categorical outcomes from a truncated geometric distribution.
#'
#' @param u A numeric vector of values between 0 and 1, typically drawn from a uniform distribution.
#' @param cdf_vals A numeric vector representing the cumulative probabilities (CDF) of a truncated and normalized geometric distribution.
#'                  Must be monotonically increasing and between 0 and 1.
#'
#' @return An integer vector of the same length as \code{u}, where each element is the corresponding ordinal category (starting at 1).
#'
#' @details Internally uses \code{findInterval()} to map values in \code{u} to bins defined by \code{cdf_vals}.
#'
#' @examples
#' cdf_vals <- compute_trunc_geom_cdf(lambda = 0.3, b = 5)
#' u <- runif(10)
#' map_to_trunc_geom_vectorized(u, cdf_vals)
#'
#' @export
map_to_trunc_geom_vectorized <- function(u, cdf_vals) {
  findInterval(u, c(0, cdf_vals)) # automatically returns ordinal categories
}
