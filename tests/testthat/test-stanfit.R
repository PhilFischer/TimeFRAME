test_that("Stationary model works", {
  m <- frame_model(matrix(c(5,0,0,5,1,1,1,1), ncol = 4))
  x <- data.frame(x = c(0.9, 1, 0.8), y = c(0.1, 0.3, 0.1))
  testthat::expect_no_error(fit_stationary(m, x))
})

test_that("Independent model works", {
  m <- frame_model(matrix(c(5,0,0,5,1,1,1,1), ncol = 4))
  x <- data.frame(x = c(0.9, 1, 0.8), y = c(0.1, 0.3, 0.1))
  testthat::expect_no_error(fit_frame(m, x))
})

test_that("Gaussian process model works", {
  m <- frame_model(matrix(c(5,0,0,5,1,1,1,1), ncol = 4))
  x <- data.frame(x = c(0.9, 1, 0.8), y = c(0.1, 0.3, 0.1))
  t <- seq(0, 2, length.out = 3)
  testthat::expect_no_error(fit_gp(m, x, t))
})

test_that("Spline model works", {
  m <- frame_model(matrix(c(5,0,0,5,1,1,1,1), ncol = 4))
  x <- data.frame(x = c(0.9, 1, 0.8), y = c(0.1, 0.3, 0.1))
  t <- seq(0, 2, length.out = 3)
  testthat::expect_no_error(fit_glm(m, x, t))
})

test_that("DGP model works", {
  m <- frame_model(matrix(c(5,0,0,5,1,1,1,1), ncol = 4))
  x <- data.frame(x = c(0.9, 1, 0.8), y = c(0.1, 0.3, 0.1))
  t <- seq(0, 2, length.out = 3)
  testthat::expect_no_error(fit_dgp(m, x, t))
})

