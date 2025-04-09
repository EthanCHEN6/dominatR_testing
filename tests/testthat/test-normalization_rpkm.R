test_that("regression_equivalence works", {
  mat <- matrix(c(100, 200, 300,
                  400, 500, 600),
                nrow = 2, byrow = TRUE)
  rownames(mat) <- c("gene1", "gene2")
  colnames(mat) <- c("sample1", "sample2", "sample3")
  gene_length <- c(1000, 2000)

  lib_size <- colSums(mat) / 1e6
  rpk <- sweep(mat, 2, lib_size, "/")
  gene_kb <- gene_length / 1000
  expected <- sweep(rpk, 1, gene_kb, "/")

  new_result <- rpkm_normalization(mat, gene_length = gene_length)
  expect_equal(new_result, expected, tolerance = 1e-6)
})
