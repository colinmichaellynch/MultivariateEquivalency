#' Multivariate Equivalency Test Based on Cramér's V
#'
#' This function performs an independent multivariate equivalency test between observed data and reference binning limits.
#' It evaluates whether each feature in a dataset deviates significantly from a uniform reference distribution
#' using Cramér's V and compares it against a Bonferroni-corrected upper bound derived from the chi-square distribution.
#'
#' @param referenceLimits A named list where each element is a numeric vector of bin edges for one variable (as produced by \code{binLimits()}).
#' @param data A numeric data frame containing the observed values, with columns matching the names in \code{referenceLimits}.
#' @param b An integer specifying the number of bins used for each variable.
#' @param alpha A numeric value for the overall significance level (default is 0.01). Bonferroni correction is applied across \code{p} variables.
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{\code{cramersV}}{The observed Cramér's V statistic for each variable.}
#'   \item{\code{cramersVUpper}}{The upper bound for Cramér's V under the null hypothesis, based on the chi-square distribution.}
#'   \item{\code{df}}{Degrees of freedom used in the chi-square test (\code{b - 1}).}
#'   \item{\code{equivalent}}{A character label indicating whether all variables are statistically equivalent to the reference (\code{"Equivalent"} or \code{"Not Equivalent"}).}
#' }
#'
#' @examples
#' set.seed(123)
#' df <- data.frame(A = rnorm(100), B = runif(100))
#' bins <- binLimits(df, b = 4)
#' multivariateEquivalency(referenceLimits = bins, data = df, b = 4)
#'
#' @export
multivariateEquivalency <- function(referenceLimits, data, b, alpha = 0.01) {
  n <- nrow(data)
  p <- length(referenceLimits)
  alphaPrime <- 1 - (1 - alpha)^(1 / p)
  chiSquareU <- qchisq(p = alphaPrime, df = b - 1, lower.tail = FALSE)
  cramersVUpper <- sqrt(chiSquareU / n) / sqrt(b - 1)
  names <- names(referenceLimits)
  cramersV <- c()

  for (i in 1:p) {
    name <- names[i]
    bins <- as.numeric(referenceLimits[[name]])
    observations <- data[, i]
    binnedData <- cut(observations, breaks = bins, include.lowest = TRUE, right = FALSE)
    observed <- table(binnedData)
    expected <- rep(n / b, b)
    chiSquare <- sum(((observed - expected)^2) / expected)
    cramersV[i] <- sqrt(chiSquare / n) / sqrt(b - 1)
  }

  if (any(cramersV > cramersVUpper)) {
    equivalent <- "Not Equivalent"
  } else {
    equivalent <- "Equivalent"
  }

  dataOutput <- data.frame(cramersV, cramersVUpper, df = b - 1, equivalent)
  return(dataOutput)
}
