test_that("CPM normalization works correctly", {
  test_data <- data.frame(
    sample1 = c(100, 200, 300, 400, 500),
    sample2 = c(50, 150, 250, 350, 450),
    sample3 = c(20, 120, 220, 320, 420),
    row.names = paste0("gene", 1:5)
  )
  ## Our primary calculation - returns a matrix
  expected_cpm <- as.data.frame(sweep(as.matrix(test_data), 2, colSums(test_data) / 1e6, "/"))
  rownames(expected_cpm) <- rownames(test_data)
  colnames(expected_cpm) <- colnames(test_data)
  expected_log_cpm <- log2(expected_cpm + 1)

  ### Using our function
  result <- cpm_normalization(test_data, log_trans = FALSE)
  result_log <- cpm_normalization(test_data, log_trans = TRUE)

  ## Check return type
  expect_true(is.data.frame(result))
  expect_true(is.data.frame(result_log))

  ## Same number of dimensions
  expect_equal(dim(result), dim(expected_cpm))
  expect_equal(dim(result_log), dim(expected_log_cpm))

  ## values should be equal
  expect_equal(as.matrix(result), as.matrix(expected_cpm), tolerance = 1e-6)
  expect_equal(as.matrix(result_log), as.matrix(expected_log_cpm), tolerance = 1e-6)
})
