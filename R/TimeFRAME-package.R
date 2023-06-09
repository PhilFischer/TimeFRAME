#' The 'TimeFRAME' package.
#'
#' @description TimeFRAME is a data analysis package that can be used for Bayesian estimation of source contributions and fractionation of isotopic measurement time series. It uses Bayesian parameter estimation with Stan to estimate uncertainty and produce posterior samples. Additionally, the package provides utility functions to simulate measurements and isotopic signature data for the study of nitrous oxide.
#'
#' @docType package
#' @name TimeFRAME-package
#' @aliases TimeFRAME
#' @useDynLib TimeFRAME, .registration = TRUE
#' @importFrom ggplot2 autoplot
#' @import methods
#' @import Rcpp
#' @import rstan
#' @import tibble
#' @import dplyr
#' @import stringr
#' @import readr
#' @rawNamespace import(tidyr, except = c(extract))
#'
#' @references Stan Development Team (2023). RStan: the R interface to Stan. R package version 2.21.8. https://mc-stan.org
#' @references L. Yu, E. Harris, D. Lewicka-Szczebak, et al., “What can we learn from N2O isotope data? – analytics, processes and modelling,” Rapid Communications in Mass Spectrometry, vol. 34, no. 20, Aug. 2020. \url{https://doi.org/10.1002/rcm.8858}

NULL
