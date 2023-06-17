#' Plot `FrameFit` object using `ggplot2`
#'
#' @description
#' Plots the model results in objects of type `FrameFit` using `ggplot2`. Stationary models use a simple median estimate per source and line ranges for credible intervals, whereas time series models plot the mean per source over time with a ribbon for the credible interval.
#'
#' @export
#' @param object an object of type `FrameFit`.
#' @param ... other arguments.
autoplot.FrameFit <- function(object, ...) {
  if(class(object)[1] != "FrameFit") stop("Object must be of type FrameFit.")

  if(object$model_name == "stationary") return(autoplot.FrameFit.stationary(object, ...))

  data <- as.data.frame(object)
  xlabel <- "Source"
  if(is.null(object$vnames)) {
    data <- data %>% rename(source = variable)
    xlabel <- "Variable"
  }
  data %>%
    ggplot2::ggplot() +
    ggplot2::geom_ribbon(ggplot2::aes(x = t, ymin = `2.5%`, ymax = `97.5%`, fill = source), alpha = 0.4) +
    ggplot2::geom_line(ggplot2::aes(x = t, y = mean, col = source), linewidth = 1) +
    ggplot2::ylim(0, 1) +
    ggplot2::facet_grid(source ~ .) +
    ggplot2::xlab("Time") + ggplot2::ylab("Source Contribution and Fractionation Weight") +
    ggplot2::labs(col = xlabel, title = sprintf("Mean and 95%% Credible Interval (%s)", object$model_name)) +
    ggplot2::guides(fill = "none") +
    ggplot2::theme_bw()
}


autoplot.FrameFit.stationary <- function(object, ...) {
  data <- as.data.frame(object)
  xlabel <- "Source"
  if(is.null(object$vnames)) {
    data <- data %>% rename(source = variable)
    xlabel <- "Variable"
  }

  data %>%
    ggplot2::ggplot() +
    ggplot2::geom_linerange(ggplot2::aes(x = source, ymin = `2.5%`, ymax = `97.5%`, col = source), linewidth = 0.8) +
    ggplot2::geom_linerange(ggplot2::aes(x = source, ymin = `25%`, ymax = `75%`, col = source), linewidth = 2) +
    ggplot2::geom_point(ggplot2::aes(x = source, y = `50%`, col = source), size = 5) +
    ggplot2::ylim(0, 1) +
    ggplot2::xlab(xlabel) + ggplot2::ylab("Source Contribution and Fractionation Weight") +
    ggplot2::labs(title = sprintf("Median, 50%% and 95%% Credible Intervals (%s)", object$model_name)) +
    ggplot2::guides(col = "none") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none")
}


