#' Extract summary statistics of `FrameFit` object as data frame.
#'
#' @description
#' Collects summary statistics from the `FrameFit` object and arranges them in a data frame.
#'
#' @export
#' @param x an object of type `FrameFit`.
#' @param row.names `NULL` or a character vector giving the row names for the data frame. Missing values are not allowed.
#' @param optional logical. If `TRUE`, setting row names and converting column names (to syntactic names: see `make.names`) is optional. Note that all of *R*'s *base* package `as.data.frame()` methods use `optional` only for column names treatment, basically with the meaning of `data.frame(*, check.names = !optional)`. See also the `make.names` argument of the `matrix` method.
#' @param ... other arguments.
#'
#' @returns Data frame having a column `variable` for the variable description names, `index` if no time points were supplied or `t` with time points, `mean` for the posterior means, `se_mean` indicating estimation error based on effective sample size, `sd` for posterior standard deviation, a list of quantiles including the median 50% quantile, `n_eff` with effective sample sizes according to the Stan documentation, and `Rhat` the Gelman-Rubin statistic according to the Stan documentation.
#'
as.data.frame.FrameFit <- function(x, row.names = NULL, optional = FALSE, ...) {
  if(class(x)[1] != "FrameFit") stop("Object must be of type FrameFit.")

  if(x$model_name == "stationary") return(as.data.frame.FrameFit.stationary(x, row.names, ...))
  if(x$model_name %in% c("frame", "gp")) pattern <- c(f = "\\1\\3-\\2", r = "\\1\\2-\\3")
  if(x$model_name %in% c("spline", "dgp", "hdgp")) pattern <- c(f = "\\1\\2-\\3", r = "\\1\\2-\\3")

  res <- rstan::summary(x$stanfit, probs = c(0.025, 0.16, 0.5, 0.84, 0.975))$summary %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "variable") %>%
    dplyr::filter(stringr::str_starts(variable, "f") | stringr::str_starts(variable, "r")) %>%
    dplyr::mutate(pattern = pattern[stringr::str_sub(variable, 1, 1)]) %>%
    dplyr::mutate(variable = stringr::str_replace(variable, "([a-zA-Z]+)\\[([0-9]+),([0-9]+)\\]", pattern)) %>%
    tidyr::separate(variable, c("variable", "index"), sep = "-") %>%
    dplyr::mutate(index = readr::parse_integer(index)) %>%
    dplyr::select(-pattern) %>%
    dplyr::arrange(index)

  if(!is.null(x$vnames)) res <- res %>% group_by(index) %>% mutate(source = factor(x$vnames, levels = x$vnames), .after = variable)
  if(!is.null(x$t)) res <- res %>% group_by(variable) %>% mutate(t = x$t, .after = index) %>% select(-index)

  return(as.data.frame(res, row.names = row.names, optional = optional))
}


as.data.frame.FrameFit.stationary <- function(x, row.names, optional, ...) {
  res <- rstan::summary(x$stanfit, probs = c(0.025, 0.16, 0.5, 0.84, 0.975))$summary %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "variable") %>%
    dplyr::filter(stringr::str_starts(variable, "f") | stringr::str_starts(variable, "r")) %>%
    dplyr::mutate(variable = stringr::str_replace(variable, "([a-zA-Z]+)\\[([0-9]+)\\]", "\\1\\2")) %>%
    tidyr::separate(variable, c("variable"), sep = "-")

  if(!is.null(x$vnames)) res <- res %>% mutate(source = factor(x$vnames, levels = x$vnames), .after = variable)

  return(res)
}


