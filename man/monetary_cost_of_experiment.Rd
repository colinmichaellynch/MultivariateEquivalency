% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/monetary_cost_of_experiment.R
\name{monetary_cost_of_experiment}
\alias{monetary_cost_of_experiment}
\title{Compute the Monetary Cost of an Experiment}
\usage{
monetary_cost_of_experiment(p, nstar, delta, omega)
}
\arguments{
\item{p}{An integer specifying the number of features used in the experiment.}

\item{nstar}{An integer representing the total number of experimental samples.}

\item{delta}{A numeric value representing the additional cost incurred per feature (relative to baseline).}

\item{omega}{A numeric value representing the base cost per sample (without added feature cost).}
}
\value{
A numeric value representing the total monetary cost of the experiment.
}
\description{
This function calculates the total monetary cost of an experiment as a function of the number of features, sample size, and per-unit costs.
The cost model assumes that the cost per sample increases linearly with the number of features selected.
}
\examples{
monetary_cost_of_experiment(p = 5, nstar = 100, delta = 0.1, omega = 50)

}
