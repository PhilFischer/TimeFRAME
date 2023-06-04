#' Computes measurement means from a vector or data frame of source contribution and fractionation parameters
#'
#' @param model Object of type FrameModel
#' @param f Data frame of `K+L` columns corresponding to source contributions and fractionation
#' @export
param_map <- function(model, f) {
  if(!inherits(model, "FrameModel")) stop("Parameter model is not of type FrameModel.")

  # Prepare input
  if(is.vector(f)) {
    if(length(f) != model$dims$K+model$dims$L) stop(sprintf("Parameter f has wrong length. Found %s, but should be %s.", length(f), model$K + model$L))

    ff <- f[1:model$dims$K]
    if(sum(ff) != 1) warning("Parameter f does not fulfill sum condition.")
    if(model$dims$L > 0) {
      rr <- f[(model$dims$K+1):(model$dims$K+model$dims$L)]
    }

  }
  else {
    if(ncol(f) != model$dims$K+model$dims$L) stop(sprintf("Parameter f has wrong number of columns. Found %s, but should be %s.", ncol(f), model$K + model$L))

    ff <- t(as.matrix(f[,1:model$dims$K]))
    if(any(colSums(ff) != 1)) warning("Parameter f does not fulfill sum condition.")
    if(model$dims$L > 0) {
      rr <- t(as.matrix(f[,(model$dims$K+1):(model$dims$K+model$dims$L)]))
    }

  }

  # Calculate mixing
  S <- t(as.matrix(model$sources[,1:model$dims$d,drop=FALSE]))
  mu <- S %*% ff

  # Calculate fractionation
  if(model$dims$L > 0) {
    A <- t(as.matrix(model$frac[,1:model$dims$d,drop=FALSE]))
    mu <- mu + A %*% log(rr)
  }

  # Post-process
  if(is.vector(f)) {
    mu <- as.vector(mu)
    names(mu) <- model$dimnames
  }
  else {
    mu <- as.matrix(t(mu))
    colnames(mu) <- model$dimnames
  }

  return(mu)
}



#' Isotope Mapping computes the source contribution and fractionation weights as the solution to a linear system of equations.
#'
#' @param model Object of type FrameModel
#' @param x Data frame of `d` columns corresponding to isotopic measurements
#' @export
isotope_map <- function(model, x) {
  if(!inherits(model, "FrameModel")) stop("Parameter model is not of type FrameModel.")

  # Prepare input
  S <- cbind(as.matrix(model$sources[,1:model$dims$d]), matrix(1, ncol = 1, nrow = model$dims$K))
  if(model$dims$L > 0) S <- rbind(S, cbind(-as.matrix(model$frac[,1:model$dims$d]), matrix(0, ncol = 1, nrow = model$dims$L)))

  # Solve with pseudo-inverse
  if(is.vector(x)) {
    if(length(x) != model$dims$d) stop(sprintf("Parameter x does not have the correct number of columns. Found %s, but should be %s.", length(x), model$d))
    Z <- as.matrix(c(x, 1))
    f <- solve(S %*% t(S), S %*% Z)
    if(model$dims$L > 0) f[(model$dims$K+1):length(f),] <- exp(-f[(model$dims$K+1):length(f),])
    f <- as.vector(f)
    names(f) <- model$vnames
  }
  else {
    if(ncol(x) != model$dims$d) stop(sprintf("Parameter x does not have the correct number of columns. Found %s, but should be %s.", length(x), model$d))
    Z <- t(as.matrix(cbind(x, 1)))
    f <- solve(S %*% t(S), S %*% Z)
    if(model$dims$L > 0) f[(model$dims$K+1):nrow(f),] <- exp(-f[(model$dims$K+1):nrow(f),])
    f <- t(f)
    colnames(f) <- model$vnames
  }

  if(any(f < 0) || any(f > 1)) warning("Estimates outside of domain [0,1] created.")
  return(f)
}
