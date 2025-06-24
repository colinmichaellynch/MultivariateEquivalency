#' Hotelling's T² Multivariate Control Chart for Process Stability
#'
#' This function evaluates the stability of a multivariate process using Hotelling's T² statistic.
#' It calculates the T² value for each observation, compares it to a control limit derived from the F-distribution,
#' and plots the results. Out-of-control observations are flagged and optionally displayed in the plot.
#'
#' @param data A numeric matrix or data frame where rows are observations and columns are variables.
#' @param alpha A numeric value representing the significance level for the control chart (default is 0.05).
#'
#' @return An integer indicating the number of observations that fall outside the control limit (i.e., potential out-of-control points).
#'
#' @details
#' The function computes the sample mean vector and covariance matrix, and then calculates the Hotelling's T² statistic for each observation:
#' \deqn{T^2 = (x_i - \bar{x})^\top S^{-1} (x_i - \bar{x})}
#' The control limit is calculated using the F-distribution:
#' \deqn{UCL = \frac{(n - 1)^2}{n} \cdot \frac{p}{n - p} \cdot F_{1 - \alpha}(p, n - p)}
#' where \eqn{n} is the number of observations and \eqn{p} is the number of variables.
#'
#' A plot is generated with the T² statistic for each observation and the control limit shown as a red dashed line.
#'
#' @examples
#' set.seed(123)
#' X <- matrix(rnorm(200), ncol = 2)
#' multivariateProcessStability(X)
#'
#' @export
multivariateProcessStability <- function(data, alpha = 0.05) {
  # Ensure it's a matrix
  data <- as.matrix(data)

  # Get dimensions
  n <- nrow(data)
  p <- ncol(data)

  # Sample mean and covariance matrix
  x_bar <- colMeans(data)
  S <- cov(data)

  # Inverse of covariance matrix
  S_inv <- solve(S)

  # Calculate T^2 for each observation
  T2_values <- apply(data, 1, function(x) {
    t(x - x_bar) %*% S_inv %*% (x - x_bar)
  })

  # Control limit using the F-distribution
  control_limit <- ((n - 1)^2 / n) * qf(1 - alpha, p, n - p) * p / (n - p)

  # Identify out-of-control points
  out_of_control <- which(T2_values > control_limit)
  num_out_of_control <- length(out_of_control)

  # Adjust plot limits
  y_max <- max(max(T2_values), control_limit) * 1.1

  # Plot
  plot(T2_values, type = "b", pch = 19,
       xlab = "Observation", ylab = expression(T^2),
       main = "Hotelling T^2 Control Chart",
       ylim = c(0, y_max))
  abline(h = control_limit, col = "red", lty = 2)
  legend("topright", legend = c("Control Limit"), col = "red", lty = 2)

  # Mark out-of-control points
  if (num_out_of_control > 0) {
    points(out_of_control, T2_values[out_of_control], col = "blue", pch = 4, cex = 1.5)
  }

  # Return number of out-of-control points
  return(num_out_of_control)
}
