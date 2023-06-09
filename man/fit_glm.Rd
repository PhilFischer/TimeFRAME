% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_glm.R
\name{fit_glm}
\alias{fit_glm}
\title{Fit generalized linear model with spline bases}
\usage{
fit_glm(
  model,
  x,
  t,
  sd = data.frame(),
  eta = 0,
  M = 8,
  M.r = 4,
  iter = 2000,
  cores = 1,
  chains = 4
)
}
\arguments{
\item{model}{an object of type \code{FrameModel}.}

\item{x}{data frame with \eqn{d} columns corresponding to isotopic measurements.}

\item{t}{vector corresponding to scaled time points with same length as data frame \code{x}.}

\item{sd}{data frame corresponding to measurement errors with same shape as \code{x} (Default empty data.frame).}

\item{eta}{numeric or vector of length \eqn{d}, measurement noise magnitude overall or per isotopic measurement (Default 0).}

\item{M}{integer, number of degrees of freedom for source contributions (Default 8).}

\item{M.r}{integer, number of degrees of freedom for fractionation (Default 4).}

\item{iter}{integer, number of sampler iterations (Default 2000).}

\item{cores}{integer, number of cores to use (Default = 1).}

\item{chains}{integer, number of chains to use (Default = 4).}
}
\description{
Fit time series model using a generalzed linear model structure with time series points expanded in natural cubic spline bases.
}
\details{
See \link{frame_model} for description of model dimensions \eqn{d}, \eqn{K} and \eqn{L}.
}
\examples{
# Define sources and fit the model
model <- frame_model(n2o_sources, n2o_frac)
fit <- fit_glm(model, n2o_test[,-1], t = n2o_test$t, eta = 5)
coef(fit)

# Using only 2 sources and 2 isotopic measurements
model2 <- frame_model(n2o_sources[1:2,c(1:2,4:5)], n2o_frac[,c(1:2,4:5)])
fit2 <- fit_glm(model2, n2o_test2[,-1], t = n2o_test2$t, eta = 5)
coef(fit2)
}
