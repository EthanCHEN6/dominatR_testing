#' RPKM Normalization
#'
#' @description
#' Normalizes gene expression using Reads Per Kilobase per Million reads (RPKM):
#'   \enumerate{
#'     \item First normalize counts by library size (i.e., column totals),
#'     \item Then divide by gene length (in kilobases),
#'     \item If \code{log_trans = TRUE}, applies \code{log2(RPKM + 1)}.
#'   }
#'
#' @param x A numeric \code{matrix}, \code{data.frame}, or \code{SummarizedExperiment}.
#' @param gene_length A numeric vector of gene lengths (required if \code{x} is matrix or data.frame).
#' @param log_trans Logical. Whether to apply log2 transform to the output.
#' @param assay_name If \code{x} is a SummarizedExperiment, the assay to normalize.
#' @param new_assay_name If \code{x} is a SummarizedExperiment, name of a new assay
#'   to store the result. If NULL, overwrites the selected assay.
#'
#' @return A numeric matrix of RPKM values (or log2-RPKM), or a modified SummarizedExperiment object.
#'
#' @examples
#' # Example using a matrix
#' mat <- matrix(c(100, 200, 300,
#'                 400, 500, 600),
#'               nrow = 2, byrow = TRUE)
#' rownames(mat) <- c("gene1", "gene2")
#' colnames(mat) <- c("sample1", "sample2", "sample3")
#' gene_length <- c(1000, 2000)  # in base pairs
#'
#' rpkm_mat <- rpkm_normalization(mat, gene_length = gene_length, log_trans = TRUE)
#'
#' # Example using SummarizedExperiment
#' library(SummarizedExperiment)
#' se <- SummarizedExperiment(assays = list(counts = mat))
#' rowData(se)$gene_length <- gene_length
#' se_rpkm <- rpkm_normalization(se, log_trans = TRUE)
#'
#' assay(se_rpkm, "counts")
#'
#' @export
rpkm_normalization <- function(x,
                               gene_length     = NULL,
                               log_trans       = FALSE,
                               assay_name      = NULL,
                               new_assay_name  = NULL) {
  if (inherits(x, "SummarizedExperiment")) {
    if (is.null(assay_name)) {
      assay_name <- assayNames(x)[1]
    }

    mat <- assay(x, assay_name)
    if (!is.numeric(mat)) stop("Assay must be numeric.")

    rd <- rowData(x)
    if (!"gene_length" %in% colnames(rd)) {
      stop("No 'gene_length' column in rowData.")
    }

    gene_len <- rd[["gene_length"]]
    if (!is.numeric(gene_len)) stop("'gene_length' must be numeric.")

    norm_mat <- .rpkm_calc(mat, gene_len, log_trans)

    if (is.null(new_assay_name)) {
      assay(x, assay_name) <- norm_mat
    } else {
      assay(x, new_assay_name) <- norm_mat
    }

    return(x)

  } else if (is.matrix(x) || is.data.frame(x)) {
    if (is.data.frame(x)) x <- as.matrix(x)
    if (!is.numeric(x)) stop("Input must be numeric.")
    if (is.null(gene_length) || length(gene_length) != nrow(x)) {
      stop("Valid numeric 'gene_length' vector required, with length matching number of rows in input.")
    }

    return(.rpkm_calc(x, gene_length, log_trans))
  } else {
    stop("Unsupported input type. Must be matrix, data.frame, or SummarizedExperiment.")
  }
}

# Internal calculation helper
.rpkm_calc <- function(mat, gene_length, log_trans = FALSE) {
  lib_size <- colSums(mat, na.rm = TRUE) / 1e6
  lib_size[lib_size <= 0] <- 1

  rpk <- sweep(mat, 2, lib_size, "/")
  gene_length_kb <- gene_length / 1000
  rpkm <- sweep(rpk, 1, gene_length_kb, "/")

  if (log_trans) {
    rpkm <- log2(rpkm + 1)
  }

  return(rpkm)
}
