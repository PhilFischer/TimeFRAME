#' Create model specification
#'
#' @description
#' Create `FrameModel` specification using source data and optionally fractionation data to be used for model fitting.
#'
#' @export
#' @param sources data frame with \eqn{d} columns corresponding to isotopic measurement dimensions and \eqn{2K} rows, the first \eqn{K} corresponding to source locations and the next \eqn{K} corresponding to source spreads.
#' @param frac data frame with \eqn{d} columns corresponding to isotopic measurement dimensions and \eqn{2L} rows,
#' the first \eqn{L} corresponding to fractionation factor means and the next \eqn{L} corresponding to their standard deviations.
#'
#' @details
#' The following dimensions are defined and must be adhered to for functions using a `FrameModel`:
#' * \eqn{d} number of isotopic measurements
#' * \eqn{K} number of sources that contribute to the mixture
#' * \eqn{L} number of fractionation effects to consider
#'
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


#' Print `FrameModel`
#'
#' @export
#' @description Print an object of type FrameModel.
#' @param x an object of type FrameModel.
#' @param ... further arguments passed to or from other methods.
print.FrameModel <- function(x, ...) {
  cat("FRAME Model\n")
  cat("===========\n\n")
  cat(sprintf("Sources: %s\n", x$dims$K))
  cat(sprintf("Fr. Factors: %s\n", x$dims$L))
  cat(sprintf("Dimensions: %s\n", x$dims$d))
  if(x$dims$d > x$dims$K+x$dims$L-1) cat(sprintf("(%s extra dofs)", x$dims$d - x$dims$K - x$dims$L + 1))
  cat("\n")
}


#' Plot `FrameModel` using ggplot2
#'
#' @export
#' @description Plot an object of type `FrameModel` using ggplot2. Dual isotope plots will be generated and if more than 2 measurement dimensions are present a list of plots will be returned.
#' @param object an object of type FrameModel.
#' @param x data matrix of isotopic measurements (Default NULL).
#' @param ... other arguments.
autoplot.FrameModel <- function(object, x = NULL, ...) {
  comb <- combn(1:object$dims$K, 2)

  p <- apply(comb, 2, function(index) {
    names <- c(X = object$dimnames[index[1]], Y = object$dimnames[index[2]])
    sources <- rbind(
      cbind(type = "min", object$sources[,index] - object$sources[,object$dims$K + index] / 2) %>% tibble::rownames_to_column("name"),
      cbind(type = "max", object$sources[,index] + object$sources[,object$dims$K + index] / 2) %>% tibble::rownames_to_column("name")
    ) %>%
      dplyr::rename(all_of(names)) %>%
      tidyr::pivot_wider(id_cols = name, names_from = type, values_from = c(X, Y))

    p <- object$sources[,index] %>%
      as.data.frame() %>%
      tibble::rownames_to_column("name") %>%
      dplyr::rename(dplyr::all_of(names)) %>%
      ggplot2::ggplot() +
      ggplot2::geom_rect(data = sources, ggplot2::aes(xmin = X_min, xmax = X_max, ymin = Y_min, ymax = Y_max, fill = name)) +
      ggplot2::geom_text(ggplot2::aes(x = X, y = Y, label = name)) +
      xlab(names[1]) + ylab(names[2])

    if(!is.null(x)) {
      data <- as.data.frame(x)[,index]
      names(data) <- c("X", "Y")
      p <- p + ggplot2::geom_point(data = data, ggplot2::aes(x = X, y = Y))
    }

    return(p + guides(fill = "none") + theme_bw())
  })

  if(length(p) == 1) return(p[[1]])
  else return(p)
}

