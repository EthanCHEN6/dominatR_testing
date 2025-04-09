test_that("Quantile normalization works correctly", {
  mat <- matrix(c(10, 50, 30,
                  40, 60, 20),
                nrow = 2, byrow = TRUE)
  rownames(mat) <- c("gene1", "gene2")
  colnames(mat) <- c("sample1", "sample2", "sample3")

  sorted <- apply(mat, 2, sort)
  means <- rowMeans(sorted)
  ranks <- apply(mat, 2, rank, ties.method = "min")
  expected <- apply(ranks, 2, function(r) means[r])
  rownames(expected) <- rownames(mat)
  colnames(expected) <- colnames(mat)

  result <- quantile_normalization(mat, log_trans = FALSE)
  expect_equal(result, expected, tolerance = 1e-6)
})
