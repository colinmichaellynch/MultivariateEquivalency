#' Multivariate Power Analysis Based on Cramér's V
#'
#' This function estimates power, required sample size, or optimal lambda for a multivariate bin-based test of equivalency.
#' The method uses simulations based on Cramér's V and a truncated geometric distribution to evaluate statistical performance across parameter values.
#' Exactly one of \code{powerThreshold}, \code{lambda}, or \code{n} must be \code{NULL}.
#'
#' @param powerThreshold Desired power level (e.g., 0.8). If \code{NULL}, this function will estimate power given \code{lambda} and \code{n}.
#' @param lambda The parameter of the truncated geometric distribution. If \code{NULL}, the function will solve for it using \code{lambdaVec}.
#' @param n Sample size. If \code{NULL}, the function will solve for it using binary search.
#' @param alpha Overall significance level (default is 0.01).
#' @param b Number of bins per variable (default is 8).
#' @param lowerCor Minimum pairwise correlation for the simulated correlation matrix (default is 0.1).
#' @param upperCor Maximum pairwise correlation for the simulated correlation matrix (default is 0.9).
#' @param lambdaVec Vector of candidate lambda values to search over (used when \code{lambda = NULL}).
#' @param varVec Vector of number of variables (features) to test across.
#' @param nVec Vector of sample sizes (used only when \code{lambda = NULL}).
#' @param simulations Number of simulations per configuration (default is 100).
#' @param iterations Number of repetitions within each simulation (default is 100).
#'
#' @return A \code{data.table} with results depending on the analysis mode:
#' \itemize{
#'   \item If estimating power: returns power for each p in \code{varVec}.
#'   \item If solving for \code{n}: returns optimal sample sizes per variable number.
#'   \item If solving for \code{lambda}: returns subset of combinations with equivalent sample size within 10\% of \code{n}.
#' }
#'
#' @import data.table
#' @export
#' @seealso \code{\link{compute_trunc_geom_cdf}}, \code{\link{run_cramersV_simulation}}
multivariatePowerAnalysis <- function(powerThreshold = NULL, lambda = NULL, n = NULL,
                                      alpha = 0.01, b = 8,
                                      lowerCor = 0.1, upperCor = 0.9,
                                      lambdaVec, varVec, nVec,
                                      simulations = 100, iterations = 100) {
  nLow <- 10
  nHigh <- 1000
  nulls <- sum(sapply(list(powerThreshold, lambda, n), is.null))

  if (nulls != 1) {
    stop("Exactly one of powerThreshold, lambda, or n must be NULL.")
  }

  ## POWER ESTIMATION MODE
  if (is.null(powerThreshold)) {
    cdf_vals <- compute_trunc_geom_cdf(lambda, b)
    powerVec <- numeric(length(varVec))

    for (u in seq_along(varVec)) {
      p <- varVec[u]
      alphaPrime <- 1 - (1 - alpha)^(1 / p)
      chiSquareU <- qchisq(p = alphaPrime, df = b - 1, lower.tail = FALSE)
      count <- run_cramersV_simulation(simulations, iterations, n, p, lambda, b, cdf_vals, lowerCor, upperCor, chiSquareU)
      powerVec[u] <- count / (iterations * simulations)
      message(sprintf("Progress: %.2f%%", 100 * u / length(varVec)))
    }

    return(data.table(variableNumber = varVec, Power = powerVec))
  }

  ## N-ESTIMATION MODE (FIXED LAMBDA)
  if (is.null(n)) {
    estimate_power <- function(n_trial, p, lambda, chiSquareU, cdf_vals) {
      count <- run_cramersV_simulation(simulations, iterations, n_trial, p, lambda, b, cdf_vals, lowerCor, upperCor, chiSquareU)
      count / (iterations * simulations)
    }

    binary_search_n <- function(p, lambda, chiSquareU, cdf_vals) {
      while (nHigh - nLow > 1) {
        nMid <- floor((nLow + nHigh) / 2)
        if (estimate_power(nMid, p, lambda, chiSquareU, cdf_vals) >= 0.8) {
          nHigh <- nMid
        } else {
          nLow <- nMid
        }
      }
      return(nHigh)
    }

    nstar <- numeric(length(varVec))
    cdf_vals <- compute_trunc_geom_cdf(lambda, b)

    for (u in seq_along(varVec)) {
      p <- varVec[u]
      alphaPrime <- 1 - (1 - alpha)^(1 / p)
      chiSquareU <- qchisq(p = alphaPrime, df = b - 1, lower.tail = FALSE)
      nstar[u] <- round(binary_search_n(p, lambda, chiSquareU, cdf_vals))
      message(sprintf("Progress: %.2f%%", 100 * u / length(varVec)))
    }

    return(data.table(variableNumber = varVec, n = nstar))
  }

  ## LAMBDA ESTIMATION MODE
  if (is.null(lambda)) {
    results <- list()

    for (q in seq_along(lambdaVec)) {
      lambda_try <- lambdaVec[q]
      cdf_vals <- compute_trunc_geom_cdf(lambda_try, b)
      nstar <- numeric(length(varVec))

      for (u in seq_along(varVec)) {
        p <- varVec[u]
        alphaPrime <- 1 - (1 - alpha)^(1 / p)
        chiSquareU <- qchisq(p = alphaPrime, df = b - 1, lower.tail = FALSE)
        nstar[u] <- round({
          estimate_power <- function(n_trial) {
            count <- run_cramersV_simulation(simulations, iterations, n_trial, p, lambda_try, b, cdf_vals, lowerCor, upperCor, chiSquareU)
            count / (iterations * simulations)
          }

          binary_search_n_local <- function() {
            low <- nLow
            high <- nHigh
            while (high - low > 1) {
              mid <- floor((low + high) / 2)
              if (estimate_power(mid) >= 0.8) {
                high <- mid
              } else {
                low <- mid
              }
            }
            return(high)
          }

          binary_search_n_local()
        })
      }

      results[[q]] <- data.table(variableNumber = varVec,
                                 nArtifact = nstar,
                                 lambda = lambda_try)
      message(sprintf("Progress: %.2f%%", 100 * q / length(lambdaVec)))
    }

    combined <- rbindlist(results)
    percent <- 0.10
    subset <- subset(combined, nArtifact > (n * (1 - percent)) & nArtifact < (n * (1 + percent)))
    return(subset)
  }
}

