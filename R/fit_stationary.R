#' Fit stationary FRAME model
#'
#' @export
#' @param model Object of type FrameModel
#' @param x Data frame corresponding to isotopic measurements
#' @param iter Number of sampler iterations (Default 2000)
fit_stationary <- function(model, x, iter = 2000) {
  if(ncol(x) != model$dims$d) stop(sprintf("Argument x has the wrong number of columns. Found %s, but %s expected.", ncol(x), model$d))

  fit <- sampling(stanmodels$stationary, data = list(
    d = model$dims$d,
    N = nrow(x),
    x = t(x),
    K = model$dims$K,
    b = t(model$sources[,1:model$dims$d]),
    sdev = t(model$sources[,(model$dims$d+1):(2*model$dims$d)]),
    sigma = rep(1, model$dims$d),
    L = model$dims$L,
    frac = if(nrow(model$frac) > 0) t(model$frac[,1:model$dims$d]) else data.frame(),
    frac_sdev = if(nrow(model$frac) > 0) t(model$frac[,(model$dims$d+1):(2*model$dims$d)]) else data.frame()
  ), refresh = 0, iter = iter)
  model$stanfit <- fit
  model$data <- x
  model$model_name <- fit@model_name
  class(model) <- c("FrameFit", class(model))
  return(model)
}
