test_that("Param map works", {
  m <- frame_model(matrix(c(0,1,1,1), ncol = 2))
  testthat::expect_equal(param_map(m, c(0.5, 0.5)), c(0.5), ignore_attr = TRUE)
  testthat::expect_vector(param_map(m, matrix(c(0.5,0.3,0.9,0.5,0.7,0.1), ncol = 2)), size = 3)

  m <- frame_model(matrix(c(0,0,1,0,1,0,1,1,1,1,1,1), ncol = 4))
  testthat::expect_equal(param_map(m, c(0.5, 0.2, 0.3)), c(0.3, 0.2), ignore_attr = TRUE)
  testthat::expect_vector(param_map(m, matrix(c(0.5,0.4,0.9,0.3,0.4,0.05,0.2,0.2,0.05), ncol = 3)), size = 3)

  m <- frame_model(matrix(c(0,0,0,1,0,0,1,0,0,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1), ncol = 6))
  testthat::expect_equal(param_map(m, c(0.5, 0.2, 0.2, 0.1)), c(0.1, 0.2, 0.2), ignore_attr = TRUE)
  testthat::expect_vector(param_map(m, matrix(c(0.5,0.4,0.8,0.2,0.3,0.1,0.2,0.2,0.05,0.1,0.1,0.05), ncol = 4)), size = 3)

  m <- frame_model(matrix(c(0,1,1,0,1,1,1,1), ncol = 4))
  testthat::expect_equal(param_map(m, c(0.5, 0.5)), c(0.5, 0.5), ignore_attr = TRUE)
  testthat::expect_vector(param_map(m, matrix(c(0.5,0.3,0.9,0.5,0.7,0.1), ncol = 2)), size = 3)

  m <- frame_model(matrix(c(0,0,1,0,1,0,1,0,0,1,1,1,1,1,1,1,1,1), ncol = 6))
  testthat::expect_equal(param_map(m, c(0.5, 0.3, 0.2)), c(0.2, 0.3, 0.5), ignore_attr = TRUE)
  testthat::expect_vector(param_map(m, matrix(c(0.5,0.4,0.9,0.3,0.4,0.05,0.2,0.2,0.05), ncol = 3)), size = 3)

  m <- frame_model(matrix(c(0,0,0,1,0,0,1,0,0,1,0,0,1,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1), ncol = 8))
  testthat::expect_equal(param_map(m, c(0.5, 0.2, 0.2, 0.1)), c(0.1, 0.2, 0.2, 0.5), ignore_attr = TRUE)
  testthat::expect_vector(param_map(m, matrix(c(0.5,0.4,0.8,0.2,0.3,0.1,0.2,0.2,0.05,0.1,0.1,0.05), ncol = 4)), size = 3)

  m <- frame_model(matrix(c(0,1,1,0,1,1,1,1), ncol = 4),
                   frac = matrix(-c(1,1,1,1), ncol = 4))
  testthat::expect_equal(param_map(m, c(0.5, 0.5, 0.7)), c(0.8567, 0.8567), tolerance = 0.0001, ignore_attr = TRUE)
  testthat::expect_vector(param_map(m, matrix(c(0.5,0.3,0.9,0.5,0.7,0.1,0.5,0.6,0.7), ncol = 3)), size = 3)

  m <- frame_model(matrix(c(0,0,1,0,1,0,1,0,0,1,1,1,1,1,1,1,1,1), ncol = 6),
                   frac = matrix(-c(1,1,1,1,1,1), ncol = 6))
  testthat::expect_equal(param_map(m, c(0.5, 0.3, 0.2, 0.7)), c(0.5567, 0.6567, 0.8567), tolerance = 0.0001, ignore_attr = TRUE)
  testthat::expect_vector(param_map(m, matrix(c(0.5,0.4,0.9,0.3,0.4,0.05,0.2,0.2,0.05,0.5,0.6,0.7), ncol = 4)), size = 3)

  m <- frame_model(matrix(c(0,0,0,1,0,0,1,0,0,1,0,0,1,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1), ncol = 8),
                   frac = matrix(-c(1,1,1,1,1,1,1,1), ncol = 8))
  testthat::expect_equal(param_map(m, c(0.5, 0.2, 0.2, 0.1, 0.7)), c(0.4567, 0.5567, 0.5567, 0.8567), tolerance = 0.0001, ignore_attr = TRUE)
  testthat::expect_vector(param_map(m, matrix(c(0.5,0.4,0.8,0.2,0.3,0.1,0.2,0.2,0.05,0.1,0.1,0.05,0.5,0.6,0.7), ncol = 5)), size = 3)

  m <- frame_model(n2o_sources, n2o_frac)
  s <- param_map(m, c(0.5, 0.4, 0.1, 1))
  testthat::expect_equal(names(s), names(n2o_sources)[1:m$dims$K])
  s <- param_map(m, matrix(c(0.5, 0.5, 0.4, 0.4, 0.1, 0.1, 1, 1), ncol = 4))
  testthat::expect_equal(colnames(s), names(n2o_sources)[1:m$dims$K])

})


test_that("Sample measurements works", {
  m <- frame_model(matrix(c(0,1,1,1), ncol = 2))
  testthat::expect_vector(sample_measurements(m, c(0.5, 0.5), 1), size = 1)
  testthat::expect_vector(sample_measurements(m, matrix(c(0.5,0.3,0.9,0.5,0.7,0.1), ncol = 2), 1)) # edge case

  m <- frame_model(matrix(c(0,0,1,0,1,0,1,1,1,1,1,1), ncol = 4))
  testthat::expect_vector(sample_measurements(m, c(0.5, 0.2, 0.3), 1), size = 2)
  testthat::expect_vector(sample_measurements(m, matrix(c(0.5,0.4,0.9,0.3,0.4,0.05,0.2,0.2,0.05), ncol = 3), 1), size = 3)

  m <- frame_model(matrix(c(0,0,0,1,0,0,1,0,0,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1), ncol = 6))
  testthat::expect_vector(sample_measurements(m, c(0.5, 0.2, 0.2, 0.1), 1), size = 3)
  testthat::expect_vector(sample_measurements(m, matrix(c(0.5,0.4,0.8,0.2,0.3,0.1,0.2,0.2,0.05,0.1,0.1,0.05), ncol = 4), 1), size = 3)

  m <- frame_model(matrix(c(0,1,1,0,1,1,1,1), ncol = 4))
  testthat::expect_vector(sample_measurements(m, c(0.5, 0.5), 1), size = 2)
  testthat::expect_vector(sample_measurements(m, matrix(c(0.5,0.3,0.9,0.5,0.7,0.1), ncol = 2), 1), size = 3)

  m <- frame_model(matrix(c(0,0,1,0,1,0,1,0,0,1,1,1,1,1,1,1,1,1), ncol = 6))
  testthat::expect_vector(sample_measurements(m, c(0.5, 0.3, 0.2), 1), size = 3)
  testthat::expect_vector(sample_measurements(m, matrix(c(0.5,0.4,0.9,0.3,0.4,0.05,0.2,0.2,0.05), ncol = 3), 1), size = 3)

  m <- frame_model(matrix(c(0,0,0,1,0,0,1,0,0,1,0,0,1,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1), ncol = 8))
  testthat::expect_vector(sample_measurements(m, c(0.5, 0.2, 0.2, 0.1), 1), size = 4)
  testthat::expect_vector(sample_measurements(m, matrix(c(0.5,0.4,0.8,0.2,0.3,0.1,0.2,0.2,0.05,0.1,0.1,0.05), ncol = 4), 1), size = 3)

  m <- frame_model(matrix(c(0,1,1,0,1,1,1,1), ncol = 4),
                   frac = matrix(c(-1,-1,1,1), ncol = 4))
  testthat::expect_vector(sample_measurements(m, c(0.5, 0.5, 0.7), 1), size = 2)
  testthat::expect_vector(sample_measurements(m, matrix(c(0.5,0.3,0.9,0.5,0.7,0.1,0.5,0.6,0.7), ncol = 3), 1), size = 3)

  m <- frame_model(matrix(c(0,0,1,0,1,0,1,0,0,1,1,1,1,1,1,1,1,1), ncol = 6),
                   frac = matrix(c(-1,-1,-1,1,1,1), ncol = 6))
  testthat::expect_vector(sample_measurements(m, c(0.5, 0.3, 0.2, 0.7), 1), size = 3)
  testthat::expect_vector(sample_measurements(m, matrix(c(0.5,0.4,0.9,0.3,0.4,0.05,0.2,0.2,0.05,0.5,0.6,0.7), ncol = 4), 1), size = 3)

  m <- frame_model(matrix(c(0,0,0,1,0,0,1,0,0,1,0,0,1,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1), ncol = 8),
                   frac = matrix(c(-1,-1,-1,-1,1,1,1,1), ncol = 8))
  testthat::expect_vector(sample_measurements(m, c(0.5, 0.2, 0.2, 0.1, 0.7), 1), size = 4)
  testthat::expect_vector(sample_measurements(m, matrix(c(0.5,0.4,0.8,0.2,0.3,0.1,0.2,0.2,0.05,0.1,0.1,0.05,0.5,0.6,0.7), ncol = 5), 1), size = 3)

  m <- frame_model(n2o_sources, n2o_frac)
  s <- sample_measurements(m, c(0.5, 0.4, 0.1, 1), 1)
  testthat::expect_equal(names(s), names(n2o_sources)[1:m$dims$K])
  s <- sample_measurements(m, matrix(c(0.5, 0.5, 0.4, 0.4, 0.1, 0.1, 1, 1), ncol = 4), 1)
  testthat::expect_equal(colnames(s), names(n2o_sources)[1:m$dims$K])

})


test_that("Isotope map works with well-defined equations", {
  m <- frame_model(matrix(c(0,1,1,1), ncol = 2))
  testthat::expect_equal(isotope_map(m, c(0.5)), c(0.5,0.5), ignore_attr = TRUE)
  testthat::expect_vector(isotope_map(m, matrix(c(0.5,0.3,0.9), ncol = 1)), size = 3)

  m <- frame_model(matrix(c(0,0,1,0,1,0,1,1,1,1,1,1), ncol = 4))
  testthat::expect_equal(isotope_map(m, c(0.5, 0.2)), c(0.3, 0.2, 0.5), ignore_attr = TRUE)
  testthat::expect_vector(isotope_map(m, matrix(c(0.5,0.4,0.9,0.3,0.4,0.05), ncol = 2)), size = 3)

  m <- frame_model(matrix(c(0,0,0,1,0,0,1,0,0,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1), ncol = 6))
  testthat::expect_equal(isotope_map(m, c(0.5, 0.2, 0.2)), c(0.1, 0.2, 0.2, 0.5), ignore_attr = TRUE)
  testthat::expect_vector(isotope_map(m, matrix(c(0.1,0.1,0.05,0.2,0.2,0.05,0.2,0.3,0.1), ncol = 3)), size = 3)

  m <- frame_model(matrix(c(0,1,1,0,1,1,1,1), ncol = 4))
  testthat::expect_equal(isotope_map(m, c(0.5, 0.5)), c(0.5, 0.5), ignore_attr = TRUE)
  testthat::expect_vector(isotope_map(m, matrix(c(0.5,0.3,0.9,0.5,0.7,0.1), ncol = 2)), size = 3)

  m <- frame_model(matrix(c(0,0,1,0,1,0,1,0,0,1,1,1,1,1,1,1,1,1), ncol = 6))
  testthat::expect_equal(isotope_map(m, c(0.5, 0.3, 0.2)), c(0.2, 0.3, 0.5), ignore_attr = TRUE)
  testthat::expect_vector(isotope_map(m, matrix(c(0.5,0.4,0.9,0.3,0.4,0.05,0.2,0.2,0.05), ncol = 3)), size = 3)

  m <- frame_model(matrix(c(0,0,0,1,0,0,1,0,0,1,0,0,1,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1), ncol = 8))
  testthat::expect_equal(isotope_map(m, c(0.5, 0.2, 0.2, 0.1)), c(0.1, 0.2, 0.2, 0.5), ignore_attr = TRUE)
  testthat::expect_vector(isotope_map(m, matrix(c(0.5,0.4,0.8,0.2,0.3,0.1,0.2,0.2,0.05,0.1,0.1,0.05), ncol = 4)), size = 3)

  m <- frame_model(matrix(c(0,1,1,0,1,1,1,1), ncol = 4),
                   frac = matrix(-c(1,1,1,1), ncol = 4))
  testthat::expect_equal(isotope_map(m, c(0.8567, 0.8567)), c(0.5, 0.5, 0.7), tolerance = 0.0001, ignore_attr = TRUE)
  testthat::expect_vector(isotope_map(m, matrix(c(0.5,0.3,0.9,0.5,0.7,0.1), ncol = 2)), size = 3)

  m <- frame_model(matrix(c(0,0,1,0,1,0,1,0,0,1,1,1,1,1,1,1,1,1), ncol = 6),
                   frac = matrix(-c(1,1,1,1,1,1), ncol = 6))
  testthat::expect_equal(isotope_map(m, c(0.5567, 0.6567, 0.8567)), c(0.5, 0.3, 0.2, 0.7), tolerance = 0.0001, ignore_attr = TRUE)
  testthat::expect_vector(isotope_map(m, matrix(c(0.9,1.1,0.15,0.9,1.3,0.15,1.2,1.3,1.1), ncol = 3)), size = 3)

  m <- frame_model(matrix(c(0,0,0,1,0,0,1,0,0,1,0,0,1,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1), ncol = 8),
                   frac = matrix(-c(1,1,1,1,1,1,1,1), ncol = 8))
  testthat::expect_equal(isotope_map(m, c(0.4567, 0.5567, 0.5567, 0.8567)), c(0.5, 0.2, 0.2, 0.1, 0.7), tolerance = 0.0001, ignore_attr = TRUE)
  testthat::expect_vector(isotope_map(m, matrix(c(0.8,0.6,0.4,0.9,0.7,0.4,0.9,0.8,0.5,1.2,0.9,1.2), ncol = 4)), size = 3)

  m <- frame_model(n2o_sources, n2o_frac)
  s <- isotope_map(m, c(-41.7, 20.2, 23.6))
  testthat::expect_equal(names(s), c(row.names(n2o_sources), row.names(n2o_frac)))
  s <- isotope_map(m, matrix(c(-41.7, -41.7, 20.2, 20.2, 23.6, 23.6), ncol = 3))
  testthat::expect_equal(colnames(s), c(row.names(n2o_sources), row.names(n2o_frac)))

})

