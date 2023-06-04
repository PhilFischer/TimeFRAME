#' Fit stationary FRAME model
#'
#' @export
#' @param model Object of type FrameModel
#' @param x Data frame corresponding to isotopic measurements
#' @param iter Number of sampler iterations (Default 2000)
fit_stationary <- function(model, x, iter = 2000) {
  sampling(stanmodels$stationary, data = list(
    N = nrow(x),
    x = t(x),
    K = nrow(model$sources),
    b = t(model$sources[,1:model$K]),
    sdev = t(model$sources[,(model$K+1):(2*model$K)]),
    sigma = c(1,1),
    L = nrow(model$frac),
    frac = if(nrow(model$frac) > 0) t(model$frac[,1:model$L]) else data.frame(),
    frac_sdev = if(nrow(model$frac) > 0) t(model$frac[,(model$K+1):(2*model$K)]) else data.frame()
  ), refresh = 0, iter = iter)
}
