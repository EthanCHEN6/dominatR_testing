#' Quantile Normalization
#'
#' @description
#' Normalizes each column to have the same distribution using quantile normalization.
#'   \enumerate{
#'     \item Sort values in each column
#'     \item Replace each rank with the mean across columns
#'     \item If \code{log_trans = TRUE}, applies \code{log2(... + 1)} transform
#'   }
#'
#' @param x A \code{matrix}, \code{data.frame}, or \code{SummarizedExperiment}.
#' @param log_trans Logical. If TRUE, applies log2 transformation to normalized data.
#' @param assay_name If \code{x} is a SummarizedExperiment, the assay to normalize.
#' @param new_assay_name If \code{x} is a SummarizedExperiment, name of a new assay
#'   to store the result. If NULL, overwrites the selected assay.
#'
#' @return A normalized matrix (or SummarizedExperiment with modified or new assay).
#'
#' @examples
#' # Example using a matrix
#' mat <- matrix(c(10, 50, 30,
#'                 40, 60, 20),
#'               nrow = 2, byrow = TRUE)
#' rownames(mat) <- c("gene1", "gene2")
#' colnames(mat) <- c("sample1", "sample2", "sample3")
#'
#' qn_mat <- quantile_normalization(mat, log_trans = TRUE)
#'
#' # Example using SummarizedExperiment
#' library(SummarizedExperiment)
#' se <- SummarizedExperiment(assays = list(counts = mat))
#' se_qn <- quantile_normalization(se, log_trans = TRUE)
#'
#' assay(se_qn, "counts")
#'
#' @export
quantile_normalization <- function(x,
                                   log_trans       = FALSE,
                                   assay_name      = NULL,
                                   new_assay_name  = NULL) {
  if (inherits(x, "SummarizedExperiment")) {
    if (is.null(assay_name)) {
      assay_name <- assayNames(x)[1]
    }

    mat <- assay(x, assay_name)
    if (!is.numeric(mat)) stop("Selected assay is not numeric.")

    norm_mat <- .quant_norm_calc(mat, log_trans)

    if (is.null(new_assay_name)) {
      assay(x, assay_name) <- norm_mat
    } else {
      assay(x, new_assay_name) <- norm_mat
    }

    return(x)

  } else if (is.matrix(x) || is.data.frame(x)) {
    if (is.data.frame(x)) {
      x <- as.matrix(x)
    }
    if (!is.numeric(x)) stop("Input must be numeric.")
    return(.quant_norm_calc(x, log_trans))

  } else {
    stop("Unsupported input type. Must be matrix, data.frame, or SummarizedExperiment.")
  }
}

# Internal quantile normalization logic
.quant_norm_calc <- function(mat, log_trans = FALSE) {
  ranks <- apply(mat, 2, rank, ties.method = "min")
  sorted <- apply(mat, 2, sort)
  means <- rowMeans(sorted)

  normalized <- apply(ranks, 2, function(r) means[r])

  rownames(normalized) <- rownames(mat)
  colnames(normalized) <- colnames(mat)

  if (log_trans) {
    normalized <- log2(normalized + 1)
  }

  return(normalized)
}
