
#' Counts Per Million Normalization
#'
#' @param object A matrix, data.frame, list of data.frames, or a SummarizedExperiment object.
#' @param log_trans Logical, whether to apply log2 transformation.
#' @return Normalized object of the same structure.
#' @export
#'
#'

setGeneric("cpm_normalization", function(x, log_trans = FALSE) {
  standardGeneric("cpm_normalization")
})

normalize_cpm <- function(df, log_trans = FALSE) {
  if (!is.numeric(as.matrix(df))) {
    stop("Input matrix must be numeric.")
  }
  sum <- colSums(df) / 1e6
  norm_df <- sweep(df, 2, sum, "/")
  if (log_trans) {
    norm_df <- log2(norm_df + 1)
  }
  return(as.data.frame(norm_df))
}


# Method for SummarizedExperiment
#' @importFrom SummarizedExperiment assay assay<-
#' @export
setMethod("cpm_normalization", "SummarizedExperiment", function(x, log_trans = FALSE) {
  assay_data <- assay(x)
  normalized <- normalize_cpm(assay_data, log_trans)
  assay(x) <- as.matrix(normalized)
  return(x)
})

# Method for data.frame
#' @export
setMethod("cpm_normalization", "data.frame", function(x, log_trans = FALSE) {
  normalize_cpm(x, log_trans)
})

# Method for matrix
#' @export
setMethod("cpm_normalization", "matrix", function(x, log_trans = FALSE) {
  normalize_cpm(x, log_trans)
})


# Method for list
#' @export
setMethod("cpm_normalization", "list", function(x, log_trans = FALSE) {
  lapply(x, function(item) {
    if (!is.matrix(item) && !is.data.frame(item)) {
      stop("List elements must be matrices or data.frames.")
    }
    normalize_cpm(item, log_trans)
  })
})
