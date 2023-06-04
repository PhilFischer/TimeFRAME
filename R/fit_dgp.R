#' Fit Dirichlet-Gaussian process prior model
#'
#' @export
#' @param model Object of type FrameModel.
#' @param x Data frame corresponding to isotopic measurements.
#' @param t Vector corresponding to time points.
#' @param sd Data frame corresponding to measurement errors (Default empty data.frame).
#' @param rho Numeric, scaled correlation length for source contributions (Default 0.1).
#' @param rho.r Numeric, scaled correlation length for fractionation (Default 0.2).
#' @param iter Numeric, number of sampler iterations (Default 2000).
fit_dgp <- function(model, x, t, sd = data.frame(), rho = 0.1, rho.r = 0.2, iter = 2000) {
  if(nrow(sd) == 0) sd <- matrix(0, nrow = nrow(x), ncol = ncol(x))
  else if(nrow(sd) != nrow(x) && ncol(sd) != ncol(x)) stop("Matrix of data and sd must have the same size, found ", dim(data), " and ", dim(sd))

  sampling(stanmodels$dgp, data = list(
    d = ncol(x),
    N = nrow(x),
    X = t(x),
    X_sd = t(sd),
    K = nrow(model$sources),
    b = t(model$sources),
    L = nrow(model$frac),
    frac = if(nrow(model$frac) > 0) t(model$frac) else data.frame(),
    sigma = 1,
    rho = rho,
    rho_r = rho.r
  ),
  iter = iter, refresh = 0, pars = c("f", "r", "S", "A", "mu"))
}
