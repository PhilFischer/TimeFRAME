#' Extract results from FrameFit object
#'
#' @export
#' @param object an object of type FrameFit.
#' @param point.est string, type of point estimate to compute. One of "mean|median" (Default "mean").
#' @param ... Other arguments.
coef.FrameFit <- function(object, point.est = "mean", ...) {
  if(class(object)[1] != "FrameFit") stop("Object must be of type FrameFit.")

  if(object$model_name == "stationary") return(coef.FrameFit.stationary(object, point.est, ...))
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
    dplyr::select(dplyr::all_of(c("variable", "index", ifelse(point.est == "mean", "mean", "50%")))) %>%
    dplyr::arrange(index)

  if(!is.null(object$vnames)) res <- res %>% group_by(index) %>% mutate(source = object$vnames, .after = variable)
  if(!is.null(object$t)) res <- res %>% group_by(variable) %>% mutate(t = object$t, .after = index) %>% select(-index)

  return(as.data.frame(res))
}


coef.FrameFit.stationary <- function(object, point.est = "mean", ...) {
  res <- rstan::summary(object$stanfit)$summary %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "variable") %>%
    dplyr::filter(stringr::str_starts(variable, "f") | stringr::str_starts(variable, "r")) %>%
    dplyr::mutate(variable = stringr::str_replace(variable, "([a-zA-Z]+)\\[([0-9]+)\\]", "\\1\\2")) %>%
    tidyr::separate(variable, c("variable"), sep = "-") %>%
    dplyr::select(dplyr::all_of(c("variable", ifelse(point.est == "mean", "mean", "50%")))) %>%
    tibble::column_to_rownames("variable")
  res <- as.matrix(t(res))
  if(!is.null(object$vnames)) colnames(res) <- object$vnames
  return(res)
}
