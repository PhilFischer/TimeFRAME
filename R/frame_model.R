#' Create FRAME model specification
#'
#' @description Create FRAME model specification using source data and optionally fractionation data to be used for model fitting.
#' @param sources data frame with \eqn{d} columns corresponding to isotopic measurement dimensions and \eqn{2K} rows, the first \eqn{K} corresponding to source locations and the next \eqn{K} corresponding to source spreads.
#' @param frac data frame with \eqn{d} columns corresponding to isotopic measurement dimensions and \eqn{2L} rows,
#' the first \eqn{L} corresponding to fractionation factor means and the next \eqn{L} corresponding to their standard deviations.
#' @export
frame_model <- function(sources, frac = data.frame()) {
  if(ncol(sources) %% 2 != 0) stop("Parameter sources must have an even number of columns, the first half corresponding to locations and the rest to spreads.")
  d <- ncol(sources) / 2
  K <- nrow(sources)
  L <- nrow(frac)
  vnames <- rownames(sources)
  if(L > 0) vnames <- c(vnames, rownames(frac))
  if(d < K+L-1) warning(sprintf("Model has %s unknowns but only %s measurements.", K+L-1, d))
  model <- list(
    sources = as.data.frame(sources),
    frac = as.data.frame(frac),
    dims = list(d = d, K = K, L = L),
    dimnames = colnames(sources)[1:d],
    vnames = vnames
  )
  class(model) <- c("FrameModel", class(model))
  return(model)
}


#' Print FRAME model
#'
#' @description Print an object of type FrameModel.
#' @param x an object of type FrameModel.
#' @param ... further arguments passed to or from other methods.
#' @export
print.FrameModel <- function(x, ...) {
  cat("FRAME Model\n")
  cat("===========\n\n")
  cat(sprintf("Sources: %s\n", x$dims$K))
  cat(sprintf("Fr. Factors: %s\n", x$dims$L))
  cat(sprintf("Dimensions: %s\n", x$dims$d))
  if(x$dims$d > x$dims$K+x$dims$L-1) cat(sprintf("(%s extra dofs)", x$dims$d - x$dims$K - x$dims$L + 1))
  cat("\n")
}
