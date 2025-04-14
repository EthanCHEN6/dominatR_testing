#' Quantile Normalization
#'
#' @description
#' Applies quantile normalization across samples (columns), enforcing the same
#' empirical distribution for each column. This helps to remove technical bias
#' across replicates in high-throughput datasets.
#'
#' @details
#' If \code{x} is a \code{SummarizedExperiment}, the function will extract the
#' assay using \code{assay_name}, apply quantile normalization, and return a
#' new or updated assay. If \code{x} is a matrix or data.frame, normalization is
#' applied directly to the input matrix.
#'
#' @param x A numeric \code{matrix} or \code{data.frame} of gene counts,
#'   or a \code{SummarizedExperiment} containing such counts.
#'   \describe{
#'     \item{If a \code{SummarizedExperiment},}{the function applies normalization
#'       to the specified assay (via \code{assay_name}).}
#'     \item{If a \code{data.frame}/\code{matrix},}{the normalization is applied directly.}
#'   }
#' @param log_trans Logical. If \code{TRUE}, apply \code{log2(... + 1)} transform
#'   to the quantile-normalized values.
#' @param assay_name If \code{x} is a SummarizedExperiment, name of the assay to
#'   normalize. Defaults to the first assay if not specified.
#' @param new_assay_name If \code{x} is a SummarizedExperiment, name of a new
#'   assay in which to store the quantile-normalized (or log2-transformed) values.
#'   If \code{NULL}, overwrites the original assay.
#'
#' @return A numeric \strong{matrix} of quantile-normalized (or log2-normalized)
#'   values if \code{x} is a data.frame or matrix. If \code{x} is a
#'   SummarizedExperiment, returns the modified SummarizedExperiment with the
#'   normalized data placed in the existing or new assay.
#'
#' @examples
#' # -------------------------------
#' # 1) Using a matrix
#' # -------------------------------
#' mat <- matrix(c(5, 4, 3,
#'                 2, 1, 6),
#'               nrow = 2, byrow = TRUE)
#' rownames(mat) <- c("gene1", "gene2")
#' colnames(mat) <- c("sample1", "sample2", "sample3")
#'
#' qn_mat <- quantile_normalization(mat)
#' qn_log <- quantile_normalization(mat, log_trans = TRUE)
#'
#' # -------------------------------
#' # 2) Using a SummarizedExperiment
#' # -------------------------------
#' if (requireNamespace("airway") && requireNamespace("SummarizedExperiment")) {
#'   library(SummarizedExperiment)
#'   data(airway, package = "airway")
#'   se <- airway
#'
#'   # Overwrite the 'counts' assay
#'   se_qn <- quantile_normalization(se)
#'   head(assay(se_qn))
#'
#'   # Store results in a new assay
#'   se_qn2 <- quantile_normalization(se, new_assay_name = "quant_counts")
#'   head(assay(se_qn2, "quant_counts"))
#'
#'   # Create a test assay and normalize
#'   mat <- assay(se)
#'   assay(se, "test_counts") <- mat
#'   se_qn3 <- quantile_normalization(se, assay_name = "test_counts", new_assay_name = "qn_test")
#'   head(assay(se_qn3, "qn_test"))
#' }
#'
#' @export
quantile_normalization <- function(x,
                                   log_trans = FALSE,
                                   assay_name = NULL,
                                   new_assay_name = NULL) {

  #---------------------------
  # SummarizedExperiment path
  #---------------------------
  if (inherits(x, "SummarizedExperiment")) {
    if (is.null(assay_name)) {
      all_assays <- SummarizedExperiment::assayNames(x)
      if (length(all_assays) < 1) {
        stop("No assays found in the SummarizedExperiment.")
      }
      assay_name <- all_assays[[1]]
    }

    mat <- SummarizedExperiment::assay(x, assay_name)
    if (is.null(mat)) {
      stop("No assay named '", assay_name, "' found in the SummarizedExperiment.")
    }
    if (!is.numeric(mat)) {
      stop("Selected assay is not numeric. Please provide numeric data for quantile normalization.")
    }

    # Perform quantile normalization
    norm_mat <- .quantile_normalize_core(mat)
    if (log_trans) {
      norm_mat <- log2(norm_mat + 1)
    }

    if (is.null(new_assay_name)) {
      assay(x, assay_name) <- norm_mat
    } else {
      assay(x, new_assay_name) <- norm_mat
    }

    return(x)

    #---------------------------
    # matrix / data.frame path
    #---------------------------
  } else if (is.data.frame(x) || is.matrix(x)) {
    if (is.data.frame(x)) {
      x <- as.matrix(x)
    }
    if (!is.numeric(x)) {
      stop("Input must be numeric matrix or data.frame.")
    }

    norm_mat <- .quantile_normalize_core(x)
    if (log_trans) {
      norm_mat <- log2(norm_mat + 1)
    }

    return(norm_mat)

  } else {
    stop("Input must be a matrix/data.frame or a SummarizedExperiment.")
  }
}

# Internal quantile normalization helper (kept local)
.quantile_normalize_core <- function(mat) {
  ranks <- apply(mat, 2, rank, ties.method = "min")
  sorted <- apply(mat, 2, sort)
  means <- rowMeans(sorted)
  norm <- apply(ranks, 2, function(r) means[r])
  rownames(norm) <- rownames(mat)
  colnames(norm) <- colnames(mat)
  return(norm)
}
