#' Computes measurement means from source contribution and fractionation weights
#'
#' @description
#' Accepts a vector or data frame of parameters corresponding to source contribution and fractionation weights. The regular mixing equation is used to compute the mean measurement value over source and fractionation means.
#'
#' @export
#' @param model an object of type `FrameModel`.
#' @param f vector with \eqn{K+L} components or data frame with \eqn{K+L} columns corresponding to source contributions and fractionation.
#'
#' @details
#' See [TimeFRAME::frame_model] for description of model dimensions \eqn{d}, \eqn{K} and \eqn{L}.
#'
#' @examples
#' model <- frame_model(n2o_sources, n2o_frac)
#' df <- param_map(model, n2o_weights[,-1])
#'
#' model2 <- frame_model(n2o_sources[1:2,c(1:2,4:5)], n2o_frac[,c(1:2,4:5)])
#' df2 <- param_map(model2, n2o_weights2[,-1])
param_map <- function(model, f) {
  if(!inherits(model, "FrameModel")) stop("Parameter model is not of type FrameModel.")

  # Prepare input
  if(is.vector(f)) {
    if(length(f) != model$dims$K+model$dims$L) stop(sprintf("Parameter f has wrong length. Found %s, but should be %s.", length(f), model$dims$K + model$dims$L))

    ff <- f[1:model$dims$K]
    if(sum(ff) != 1) warning("Parameter f does not fulfill sum condition.")
    if(model$dims$L > 0) {
      rr <- f[(model$dims$K+1):(model$dims$K+model$dims$L)]
    }

  }
  else {
    if(ncol(f) != model$dims$K+model$dims$L) stop(sprintf("Parameter f has wrong number of columns. Found %s, but should be %s.", ncol(f), model$dims$K + model$dims$L))

    ff <- t(as.matrix(f[,1:model$dims$K]))
    if(!all.equal(colSums(ff), rep(1.0, ncol(ff)), tolerance = 1e-6)) warning("Parameter f does not fulfill sum condition.")
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


#' Sample isotopic measurements from true parameter values
#'
#' @description
#' Uses distributions of source isotopic signature and fractionation factor in object of type `FrameModel` to sample parameters and then simulates measurements from them using noise defined in `eta`.
#'
#' @export
#' @param model an object of type `FrameModel`.
#' @param f vector with \eqn{K+L} components or data frame with \eqn{K+L} columns containing true parameter values.
#' @param eta numeric or vector of length \eqn{d}, measurement noise overall or per isotopic measurement.
#'
#' @details
#' See [frame_model] for description of model dimensions \eqn{d}, \eqn{K} and \eqn{L}.
#'
#' @examples
#' model <- frame_model(n2o_sources, n2o_frac)
#' df <- sample_measurements(model, n2o_weights[,-1], 5)
#'
#' model2 <- frame_model(n2o_sources[1:2,c(1:2,4:5)], n2o_frac[,c(1:2,4:5)])
#' df2 <- sample_measurements(model2, n2o_weights2[,-1], 5)
sample_measurements <- function(model, f, eta) {
  if(!inherits(model, "FrameModel")) stop("Parameter model is not of type FrameModel.")

  # Prepare input
  if(is.vector(f)) {
    if(length(f) != model$dims$K+model$dims$L) stop(sprintf("Parameter f has wrong length. Found %s, but should be %s.", length(f), model$dims$K + model$dims$L))

    ff <- f[1:model$dims$K]
    if(sum(ff) != 1) warning("Parameter f does not fulfill sum condition.")
    if(model$dims$L > 0) {
      rr <- f[(model$dims$K+1):(model$dims$K+model$dims$L)]
    }

  }
  else {
    if(ncol(f) != model$dims$K+model$dims$L) stop(sprintf("Parameter f has wrong number of columns. Found %s, but should be %s.", ncol(f), model$dims$K + model$dims$L))

    ff <- t(as.matrix(f[,1:model$dims$K]))
    if(!all.equal(colSums(ff), rep(1.0, ncol(ff)), tolerance = 1e-6)) warning("Parameter f does not fulfill sum condition.")
    if(model$dims$L > 0) {
      rr <- t(as.matrix(f[,(model$dims$K+1):(model$dims$K+model$dims$L)]))
    }

  }

  # Calculate mixing
  S <- apply(model$sources, 1, function(x) {
    smin <- head(x, model$dims$d) - tail(x, model$dims$d) / 2
    smax <- head(x, model$dims$d) + tail(x, model$dims$d) / 2
    return(runif(model$dims$d, smin, smax))
  })
  mu <- S %*% ff

  # Calculate fractionation
  if(model$dims$L > 0) {
    A <- apply(model$frac, 1, function(x) {
      amean <- head(x, model$dims$d)
      asd <- tail(x, model$dims$d)
      return(rnorm(model$dims$d, amean, asd))
    })
    mu <- mu + A %*% log(rr)
  }

  # Post-process
  if(is.vector(f)) {
    mu <- as.vector(mu) + rnorm(model$dims$d, mean = 0, sd = eta)
    names(mu) <- model$dimnames
  }
  else {
    mu <- t(apply(mu, 2, function(x) x + rnorm(model$dims$d, mean = 0, sd = eta)))
    colnames(mu) <- model$dimnames
  }

  return(mu)
}


#' Solve for source contributions and fractionation via isotope mapping
#'
#' @description
#' Isotope Mapping computes the source contribution and fractionation weights as the solution to a linear system of equations. In certain cases where no solution exists due to geometric constraints the function will return unconstrained estimates and throw a warning and if no unconstrained solution exists due to collinearity there will be an error.
#'
#' @export
#' @param model an object of type `FrameModel`.
#' @param x vector with \eqn{d} components or data frame with \eqn{d} columns corresponding to isotopic measurements.
#'
#' @details
#' See [TimeFRAME::frame_model] for description of model dimensions \eqn{d}, \eqn{K} and \eqn{L}.
#'
#' @examples
#' model <- frame_model(n2o_sources, n2o_frac)
#' df <- isotope_map(model, n2o_test[,-1])
#'
#' model2 <- frame_model(n2o_sources[1:2,c(1:2,4:5)], n2o_frac[,c(1:2,4:5)])
#' df2 <- isotope_map(model2, n2o_test2[,-1])
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
