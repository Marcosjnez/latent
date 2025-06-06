% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getmodel.R
\name{getmodel}
\alias{getmodel}
\title{Get the default model for Latent Class Analysis.}
\usage{
getmodel(data, model = rep("multinomial", ncol(data)), nclasses = 2L)
}
\arguments{
\item{data}{data.frame or matrix of response.}

\item{item}{Character vector with the model for each item.}

\item{nclasses}{Number of latent classes.}

\item{model}{List of parameter labels. See 'details' for more information.}

\item{constraints}{Should the model be checked for identification? Defaults to TRUE.}
}
\value{
List with the following objects:
\item{none}{.}
\item{none}{.}
}
\description{
Get the default model for Latent Class Analysis.
}
\details{
\code{getmodel} generates the model for the probability of belonging
to the classes and the conditional response probabilities. These models may
be modified by the user to set equality constraints or to fix parameters.
}
\references{
None yet.
}
