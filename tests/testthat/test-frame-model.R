test_that("FRAME model can be created", {
  testthat::expect_silent(frame_model(n2o_sources))
  testthat::expect_error(frame_model(cbind(n2o_sources, X = c(1,1,1))))
})

test_that("FRAME model with fractionation can be created", {
  testthat::expect_silent(frame_model(n2o_sources, frac = n2o_frac))
  testthat::expect_error(frame_model(cbind(n2o_sources, X = c(1,1,1)), frac = n2o_frac))
})

test_that("FRAME model warns for not enough dofs", {
  testthat::expect_warning(frame_model(matrix(c(1,1,1,1,2,3), ncol = 2)))
  testthat::expect_no_warning(frame_model(matrix(c(1,1,1,1,2,3,4,5), ncol = 4)))
  testthat::expect_warning(frame_model(matrix(c(1,1,1,1,2,3,4,5), ncol = 4), frac = matrix(c(1,1,1,1,2,2,2,2), ncol = 4)))
  testthat::expect_no_warning(frame_model(matrix(c(1,1,1,1,2,3,4,5), ncol = 4), frac = matrix(c(1,1,1,1), ncol = 4)))
})
