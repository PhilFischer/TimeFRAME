test_that("Stationary model works", {
  m <- frame_model(matrix(c(5,0,0,5,1,1,1,1), ncol = 4))
  x <- data.frame(x = c(0.9, 1, 0.8), y = c(0.1, 0.3, 0.1))
  testthat::expect_no_error(fit <- fit_stationary(m, x))
  testthat::expect_vector(coef(fit, point.est = "mean"))
  testthat::expect_vector(coef(fit, point.est = "median"))
})

test_that("Stationary model works with fractionation", {
  m <- frame_model(matrix(c(5,0,0,5,1,1,1,1), ncol = 4), frac = matrix(c(-1,-1,1,1), ncol = 4))
  x <- data.frame(x = c(0.9, 1, 0.8), y = c(0.1, 0.3, 0.1))
  testthat::expect_no_error(fit <- fit_stationary(m, x))
  testthat::expect_vector(coef(fit, point.est = "mean"))
  testthat::expect_vector(coef(fit, point.est = "median"))
})


test_that("Independent model works", {
  d <- data.frame(X = c(5,0), Y = c(0, 5), DX = c(1,1), DY = c(1,1), row.names = c("A", "B"))
  m <- frame_model(d)
  x <- data.frame(x = c(0.9, 1, 0.8), y = c(0.1, 0.3, 0.1))
  testthat::expect_no_error(fit <- fit_frame(m, x))
  testthat::expect_vector(coef(fit, point.est = "mean"))
  testthat::expect_vector(coef(fit, point.est = "median"))
})

test_that("Indenendent model works with fractionation", {
  d <- data.frame(X = c(5,0), Y = c(0, 5), DX = c(1,1), DY = c(1,1), row.names = c("A", "B"))
  f <- data.frame(X = c(-1), Y = c(-1), sX = c(1), sY = c(1), row.names = c("F"))
  m <- frame_model(d, frac = f)
  x <- data.frame(x = c(0.9, 1, 0.8), y = c(0.1, 0.3, 0.1))
  testthat::expect_no_error(fit <- fit_frame(m, x))
  testthat::expect_vector(coef(fit, point.est = "mean"))
  testthat::expect_vector(coef(fit, point.est = "median"))
})


test_that("GP model works", {
  d <- data.frame(X = c(5,0), Y = c(0, 5), DX = c(1,1), DY = c(1,1), row.names = c("A", "B"))
  m <- frame_model(d)
  x <- data.frame(x = c(0.9, 1, 0.8), y = c(0.1, 0.3, 0.1))
  t <- seq(0, 1, length.out = 3)
  testthat::expect_no_error(fit <- fit_gp(m, x, t))
  testthat::expect_vector(coef(fit, point.est = "mean"))
  testthat::expect_vector(coef(fit, point.est = "median"))
})

test_that("GP model works with fractionation", {
  d <- data.frame(X = c(5,0), Y = c(0, 5), DX = c(1,1), DY = c(1,1), row.names = c("A", "B"))
  f <- data.frame(X = c(-1), Y = c(-1), sX = c(1), sY = c(1), row.names = c("F"))
  m <- frame_model(d, frac = f)
  x <- data.frame(x = c(0.9, 1, 0.8), y = c(0.1, 0.3, 0.1))
  t <- seq(0, 1, length.out = 3)
  testthat::expect_no_error(fit <- fit_gp(m, x, t))
  testthat::expect_vector(coef(fit, point.est = "mean"))
  testthat::expect_vector(coef(fit, point.est = "median"))
})


test_that("Spline model works", {
  d <- data.frame(X = c(5,0), Y = c(0, 5), DX = c(1,1), DY = c(1,1), row.names = c("A", "B"))
  m <- frame_model(d)
  x <- data.frame(x = c(0.9, 1, 0.8), y = c(0.1, 0.3, 0.1))
  t <- seq(0, 1, length.out = 3)
  testthat::expect_no_error(fit <- fit_glm(m, x, t))
  testthat::expect_vector(coef(fit, point.est = "mean"))
  testthat::expect_vector(coef(fit, point.est = "median"))
})

test_that("Spline model works with fractionation", {
  d <- data.frame(X = c(5,0), Y = c(0, 5), DX = c(1,1), DY = c(1,1), row.names = c("A", "B"))
  f <- data.frame(X = c(-1), Y = c(-1), sX = c(1), sY = c(1), row.names = c("F"))
  m <- frame_model(d, frac = f)
  x <- data.frame(x = c(0.9, 1, 0.8), y = c(0.1, 0.3, 0.1))
  t <- seq(0, 1, length.out = 3)
  testthat::expect_no_error(fit <- fit_glm(m, x, t))
  testthat::expect_vector(coef(fit, point.est = "mean"))
  testthat::expect_vector(coef(fit, point.est = "median"))
})


test_that("DGP model works", {
  d <- data.frame(X = c(5,0), Y = c(0, 5), DX = c(1,1), DY = c(1,1), row.names = c("A", "B"))
  m <- frame_model(d)
  x <- data.frame(x = c(0.9, 1, 0.8), y = c(0.1, 0.3, 0.1))
  t <- seq(0, 1, length.out = 3)
  testthat::expect_no_error(fit <- fit_dgp(m, x, t))
  testthat::expect_vector(coef(fit, point.est = "mean"))
  testthat::expect_vector(coef(fit, point.est = "median"))
})

test_that("DGP model works with fractionation", {
  d <- data.frame(X = c(5,0), Y = c(0, 5), DX = c(1,1), DY = c(1,1), row.names = c("A", "B"))
  f <- data.frame(X = c(-1), Y = c(-1), sX = c(1), sY = c(1), row.names = c("F"))
  m <- frame_model(d, frac = f)
  x <- data.frame(x = c(0.9, 1, 0.8), y = c(0.1, 0.3, 0.1))
  t <- seq(0, 1, length.out = 3)
  testthat::expect_no_error(fit <- fit_dgp(m, x, t))
  testthat::expect_vector(coef(fit, point.est = "mean"))
  testthat::expect_vector(coef(fit, point.est = "median"))
})

