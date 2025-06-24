#' Compute Percentile-Based Bin Limits for Each Variable
#'
#' This function divides each variable (column) of a numeric data frame into \code{b} bins based on percentile cutoffs.
#' It returns a list of bin edge vectors, one for each column, where the first and last edges are set to \code{-Inf} and \code{Inf}, respectively, to cover all values.
#'
#' @param data A data frame or list where each element (column) is a numeric vector to be binned.
#' @param b An integer specifying the number of bins (must be >= 1).
#'
#' @return A list of numeric vectors. Each vector contains \code{b + 1} bin edges for the corresponding column in \code{data}.
#'         The first edge is \code{-Inf} and the last is \code{Inf}, ensuring all values fall into some bin. Values must be real numbers.
#'
#' @examples
#' df <- data.frame(a = rnorm(100), b = runif(100))
#' binLimits(df, b = 4)
#'
#' @export
binLimits <- function(data, b) {
  # Define the percentile breakpoints
  percentileVec <- seq(0, 1, length.out = b + 1)

  # Compute bin edges for each column
  bin_limits_list <- lapply(data, function(col) {
    bins <- quantile(col, percentileVec, na.rm = TRUE)
    bins[1] <- -Inf
    bins[length(bins)] <- Inf
    return(bins)
  })

  return(bin_limits_list)
}
