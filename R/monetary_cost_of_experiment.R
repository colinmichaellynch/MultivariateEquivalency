#' Compute the Monetary Cost of an Experiment
#'
#' This function calculates the total monetary cost of an experiment as a function of the number of features, sample size, and per-unit costs.
#' The cost model assumes that the cost per sample increases linearly with the number of features selected.
#'
#' @param p An integer specifying the number of features used in the experiment.
#' @param nstar An integer representing the total number of experimental samples.
#' @param delta A numeric value representing the additional cost incurred per feature (relative to baseline).
#' @param omega A numeric value representing the base cost per sample (without added feature cost).
#'
#' @return A numeric value representing the total monetary cost of the experiment.
#'
#' @examples
#' monetary_cost_of_experiment(p = 5, nstar = 100, delta = 0.1, omega = 50)
#'
#' @export
monetary_cost_of_experiment <- function(p, nstar, delta, omega) {
  cost <- nstar * omega * (1 + p * delta)
  return(cost)
}
