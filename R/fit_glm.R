#' Fit generalized linear model model with spline bases
#'
#' @export
#' @import splines
#' @param model Object of type FrameModel.
#' @param x Data frame corresponding to isotopic measurements.
#' @param t Vector corresponding to time steps.
#' @param sd Data frame corresponding to measurement errors (Default empty data.frame).
#' @param M Integer, number of degrees of freedom for source contributions (Default 8).
#' @param M.r Integer, number of degrees of freedom for fractionation (Default 4).
#' @param iter Numeric, number of sampler iterations.
fit_glm <- function(model, x, t, sd = data.frame(), M = 8, M.r = 4, iter = 2000) {
  if(nrow(sd) == 0) sd <- matrix(0, nrow = nrow(x), ncol = ncol(x))
  else if(nrow(sd) != nrow(x) && ncol(sd) != ncol(x)) stop("Matrix of data and sd must have the same size, found ", dim(data), " and ", dim(sd))
  if(nrow(x) != length(t)) stop("Data and time points have different lengths.")

  basis <- splines::ns(t, df = M, intercept = TRUE)
  basis_r <- splines::ns(t, df = M.r, intercept = TRUE)
  basis <- basis / mean(diag(basis %*% t(basis)))
  basis_r <- basis_r / mean(diag(basis_r %*% t(basis_r)))

  fit <- sampling(stanmodels$spline, data = list(
    d = ncol(x),
    N = nrow(x),
    X = t(x),
    X_sd = t(sd),
    K = nrow(model$sources),
    b = t(model$sources),
    L = nrow(model$frac),
    frac = if(nrow(model$frac) > 0) t(model$frac) else data.frame(),
    sigma = 1,
    M = M,
    M_r = M.r,
    basis = t(basis),
    basis_r = t(basis_r)
  ),
  iter = iter, refresh = 0, pars = c("f", "r", "S", "A", "mu"))

  model$stanfit <- fit
  model$data <- x
  model$t <- t
  model$model_name <- fit@model_name
  class(model) <- c("FrameFit", class(model))

  return(model)
}
