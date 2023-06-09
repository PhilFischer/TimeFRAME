#' Fit independent time step model
#'
#' @description
#' Fit a time series model using FRAME at each point in time independently.
#'
#' @export
#' @param model an object of type `FrameModel`.
#' @param x data frame with \eqn{d} columns corresponding to isotopic measurements.
#' @param t vector corresponding to scaled time points with same length as data frame `x` (Default NULL).
#' @param sd data frame corresponding to measurement errors with same shape as `x` (Default empty data.frame).
#' @param eta numeric or vector of length \eqn{d}, measurement noise magnitude overall or per isotopic measurement (Default 0).
#' @param iter integer, number of sampler iterations (Default 2000).
#' @param cores integer, number of cores to use (Default = 1).
#' @param chains integer, number of chains to use (Default = 4).
#'
#' @details
#' See [TimeFRAME::frame_model] for description of model dimensions \eqn{d}, \eqn{K} and \eqn{L}.
#'
#' @examples
#' # Define sources and fit the model
#' model <- frame_model(n2o_sources, n2o_frac)
#' fit <- fit_frame(model, n2o_test[,-1], eta = 5)
#' coef(fit)
#'
#' # Using only 2 sources and 2 isotopic measurements
#' model2 <- frame_model(n2o_sources[1:2,c(1:2,4:5)], n2o_frac[,c(1:2,4:5)])
#' fit2 <- fit_frame(model2, n2o_test2[,-1], eta = 5)
#' coef(fit2)
fit_frame <- function(model, x, t = NULL, sd = data.frame(), eta = 0, iter = 2000, cores = 1, chains = 4) {
  if(nrow(sd) == 0) sd <- matrix(0, nrow = nrow(x), ncol = ncol(x))
  else if(nrow(sd) != nrow(x) && ncol(sd) != ncol(x)) stop("Matrix of data and sd must have the same size, found ", dim(data), " and ", dim(sd))
  if((length(eta) != 1) && (length(eta) != model$dims$d)) stop(sprintf("Argument eta has wrong length. Found %s, but %s expected.", length(eta), model$dims$d))

  if(length(eta) == 1) eta <- rep(eta, model$dims$d)
  if(any(sd + eta == 0)) stop("Standard deviations contain zeros, the model cannot be initialized.")

  fit <- sampling(stanmodels$frame, data = list(
    d = ncol(x),
    N = nrow(x),
    X = t(x),
    X_sd = t(sd),
    K = nrow(model$sources),
    b = t(model$sources),
    L = nrow(model$frac),
    frac = if(nrow(model$frac) > 0) t(model$frac) else data.frame(),
    eta = eta,
    alpha = 1,
    alpha_r = 1
  ),
  iter = iter, refresh = 0, cores = cores, chains = chains, pars = c("f", "r", "S", "A", "mu"))

  model$stanfit <- fit
  model$data <- x
  model$t <- t
  model$model_name <- fit@model_name
  class(model) <- c("FrameFit", class(model))

  return(model)
}
