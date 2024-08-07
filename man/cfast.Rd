% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cfast.R
\name{cfast}
\alias{cfast}
\title{Confirmatory factor analysis.}
\usage{
cfast(
  data,
  lambda,
  phi = NULL,
  psi = NULL,
  cor = "pearson",
  estimator = "uls",
  rotate = NULL,
  missing = "pairwise.complete.cases",
  nobs = NULL,
  group = NULL,
  invariance = "metric",
  control = NULL,
  positive = FALSE,
  random_starts = 1L,
  cores = 1L
)
}
\description{
Confirmatory factor analysis.
}
