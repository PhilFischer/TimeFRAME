#' Extract summary statistics of `FrameFit` object as data frame.
#'
#' @description
#' Collects summary statistics from the `FrameFit` object and arranges them in a data frame.
#'
#' @export
#' @param object an object of type `FrameFit`.
#' @param ... other arguments.
#'
#' @returns Data frame having a column `variable` for the variable description names, `index` if no time points were supplied or `t` with time points, `mean` for the posterior means, `se_mean` indicating estimation error based on effective sample size, `sd` for posterior standard deviation, a list of quantiles including the median 50% quantile, `n_eff` with effective sample sizes according to the Stan documentation, and `Rhat` the Gelman-Rubin statistic according to the Stan documentation.
#'
as.data.frame.FrameFit <- function(object, ...) {
  if(class(object)[1] != "FrameFit") stop("Object must be of type FrameFit.")

  if(object$model_name == "stationary") return(as.data.frame.FrameFit.stationary(object, ...))
  if(object$model_name %in% c("frame", "gp")) pattern <- c(f = "\\1\\3-\\2", r = "\\1\\2-\\3")
  if(object$model_name %in% c("spline", "dgp", "hdgp")) pattern <- c(f = "\\1\\2-\\3", r = "\\1\\2-\\3")

  res <- rstan::summary(object$stanfit)$summary %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "variable") %>%
    dplyr::filter(stringr::str_starts(variable, "f") | stringr::str_starts(variable, "r")) %>%
    dplyr::mutate(pattern = pattern[stringr::str_sub(variable, 1, 1)]) %>%
    dplyr::mutate(variable = stringr::str_replace(variable, "([a-zA-Z]+)\\[([0-9]+),([0-9]+)\\]", pattern)) %>%
    tidyr::separate(variable, c("variable", "index"), sep = "-") %>%
    dplyr::mutate(index = readr::parse_integer(index)) %>%
    dplyr::select(-pattern) %>%
    dplyr::arrange(index)

  if(!is.null(object$vnames)) res <- res %>% group_by(index) %>% mutate(source = factor(object$vnames, levels = object$vnames), .after = variable)
  if(!is.null(object$t)) res <- res %>% group_by(variable) %>% mutate(t = object$t, .after = index) %>% select(-index)

  return(as.data.frame(res))
}


as.data.frame.FrameFit.stationary <- function(object, ...) {
  res <- rstan::summary(object$stanfit)$summary %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "variable") %>%
    dplyr::filter(stringr::str_starts(variable, "f") | stringr::str_starts(variable, "r")) %>%
    dplyr::mutate(variable = stringr::str_replace(variable, "([a-zA-Z]+)\\[([0-9]+)\\]", "\\1\\2")) %>%
    tidyr::separate(variable, c("variable"), sep = "-")

  if(!is.null(object$vnames)) res <- res %>% mutate(source = factor(object$vnames, levels = object$vnames), .after = variable)

  return(res)
}


