#' Run Parallel Simulations of Cramér's V for Power Estimation
#'
#' This function runs multiple simulations to estimate the probability of detecting a significant deviation from uniformity
#' in ordinal data using Cramér's V. Each simulation draws samples from a truncated geometric distribution with optional correlation
#' between variables, then evaluates whether any variable exceeds the critical value of Cramér's V.
#'
#' @param simulations Number of independent simulation sets to run.
#' @param iterations Number of repetitions per simulation (total runs = \code{simulations × iterations}).
#' @param n Sample size for each simulation.
#' @param p Number of variables (features).
#' @param lambda Parameter of the truncated geometric distribution.
#' @param b Number of bins (categories) for discretizing the distribution.
#' @param cdf_vals A numeric vector giving the cumulative distribution function of the truncated geometric distribution
#'        (typically from \code{\link{compute_trunc_geom_cdf}}).
#' @param lowerCor Minimum correlation coefficient used to generate positive-definite correlation matrices (default: 0.1).
#' @param upperCor Maximum correlation coefficient used to generate positive-definite correlation matrices (default: 0.9).
#' @param chiSquareU The upper threshold for the chi-square statistic under the null hypothesis (used to compute critical Cramér's V).
#'
#' @return An integer count indicating the number of simulation runs (out of total) where at least one variable's Cramér's V exceeds its critical value.
#'
#' @details
#' For multivariate data (\code{p > 1}), the function generates positively correlated latent variables using Cholesky decomposition.
#' These are transformed to uniform [0,1] variables using the normal CDF, then mapped to ordinal categories using a truncated geometric CDF.
#' Cramér's V is computed for each variable, and if any exceed the critical threshold, the run is counted as a "positive."
#' The function uses parallel computation to speed up execution.
#'
#' @examples
#' cdf <- compute_trunc_geom_cdf(lambda = 0.3, b = 5)
#' run_cramersV_simulation(
#'   simulations = 2, iterations = 5, n = 100, p = 3, lambda = 0.3, b = 5,
#'   cdf_vals = cdf, lowerCor = 0.2, upperCor = 0.8, chiSquareU = qchisq(0.95, df = 4)
#' )
#'
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @import foreach
#' @export
run_cramersV_simulation <- function(simulations, iterations, n, p, lambda, b, cdf_vals,
                                    lowerCor = 0.1, upperCor = 0.9, chiSquareU) {
  total_runs <- simulations * iterations
  cramersVUpper <- sqrt(chiSquareU / n) / sqrt(b - 1)
  expected <- rep(n / b, b)
  cramersV <- numeric(p)

  if (p > 1) {
    positive_corr_matrix <- generate_positive_correlation_matrix(p, lower = lowerCor, upper = upperCor)
    chol_decomp <- chol(positive_corr_matrix)
  }

  # Setup parallel cluster
  cores <- parallel::detectCores(logical = FALSE) - 1
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)

  count <- foreach::foreach(run = 1:total_runs, .combine = "+",
                            .packages = c("stats"),
                            .export = c("map_to_trunc_geom_vectorized")) %dopar% {
                              if (p > 1) {
                                latent_matrix <- matrix(rnorm(n * p), nrow = n, ncol = p)
                                latent_vars <- latent_matrix %*% chol_decomp
                              } else {
                                latent_vars <- matrix(rnorm(n), ncol = 1)
                              }

                              latent_uniform_vars <- pnorm(latent_vars)

                              for (k in 1:p) {
                                ordinal_X <- map_to_trunc_geom_vectorized(latent_uniform_vars[, k], cdf_vals)
                                sample <- table(factor(ordinal_X, levels = 1:b))
                                chiSquare <- sum(((sample - expected)^2) / expected)
                                cramersV_k <- sqrt(chiSquare / n) / sqrt(b - 1)
                                if (cramersV_k > cramersVUpper) {
                                  return(1)
                                }
                              }

                              return(0)
                            }

  parallel::stopCluster(cl)
  return(count)
}
