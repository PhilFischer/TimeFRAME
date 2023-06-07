#' Fit Gaussian process prior on measurements
#'
#' @description
#' Fit time series model using Gaussian process priors on measurements and otherwise fitting source contribution and fractionation weights independently.
#'
#' @export
#' @param model an object of type `FrameModel`.
#' @param x data frame with \eqn{d} columns corresponding to isotopic measurements.
#' @param t vector corresponding to scaled time points with same length as data frame `x`.
#' @param sd data frame corresponding to measurement errors with same shape as `x` (Default empty data.frame).
#' @param eta numeric or vector of length \eqn{d}, measurement noise magnitude overall or per isotopic measurement (Default 0).
#' @param rho numeric, scaled correlation length for Gaussian process (Default 0.5).
#' @param iter integer, number of sampler iterations (Default 2000).
#' @param cores integer, number of cores to use (Default = 1).
#' @param chains integer, number of chains to use (Default = 4).
#'
#' @details
#' See [TimeFRAME::frame_model] for description of model dimensions \eqn{d}, \eqn{K} and \eqn{L}.
#'
fit_gp <- function(model, x, t, sd = data.frame(), eta = 0, rho = 0.5, iter = 2000, cores = 1, chains = 4) {
  if(nrow(sd) == 0) sd <- matrix(0, nrow = nrow(x), ncol = ncol(x))
  else if(nrow(sd) != nrow(x) && ncol(sd) != ncol(x)) stop("Matrix of data and sd must have the same size, found ", dim(data), " and ", dim(sd))
  if(nrow(x) != length(t)) stop("Data and time points have different lengths.")
  if((length(eta) != 1) && (length(eta) != model$dims$d)) stop(sprintf("Argument eta has wrong length. Found %s, but %s expected.", length(eta), model$dims$d))

  if(length(eta) == 1) eta <- rep(eta, model$dims$d)
  if(any(sd + eta == 0)) stop("Standard deviations contain zeros, the model cannot be initialized.")

  fit <- sampling(stanmodels$gp, data = list(
    d = ncol(x),
    N = nrow(x),
    X = t(x),
    X_sd = t(sd),
    K = nrow(model$sources),
    b = t(model$sources),
    L = nrow(model$frac),
    frac = if(nrow(model$frac) > 0) t(model$frac) else data.frame(),
    sigma = 1,
    rho = rho
  ),
  iter = iter, refresh = 0, cores = cores, chains = chains, pars = c("f", "r", "S", "A", "mu"))

  model$stanfit <- fit
  model$data <- x
  model$t <- t
  model$model_name <- fit@model_name
  class(model) <- c("FrameFit", class(model))

  return(model)
}
