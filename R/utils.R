#' Get Normalized Counts
#'
#' Returns normalized counts using effective library sizes.
#'
#' @param x A `SeqCountData` object.
#' @return A normalized count matrix.
#'
#' @examples
#' counts <- matrix(1:12, nrow = 3)
#' colnames(counts) <- paste0("sample_", seq_len(ncol(counts)))
#' rownames(counts) <- paste0("gene_", seq_len(nrow(counts)))
#' exp_design <- data.frame(group = rep(c("A", "B"), each = 2))
#' x <- newSeqCountData(
#'     counts,
#'     exp_design = exp_design,
#'     norm_factors = rep(1, ncol(counts))
#' )
#' norm_counts(x)
#'
#' @export
norm_counts <- function(x) {
    if (is.null(x@meta$nf)) {
        stop("Normalization factors are not set. Run `calc_nf()` first.", call. = FALSE)
    }
    sweep(x@data, 2, .calc_sizefactors(x@data, x@meta$nf), "/")
}

.require_namespace <- function(package) {
    if (!requireNamespace(package, quietly = TRUE)) {
        stop("Package `", package, "` is required for this function.", call. = FALSE)
    }
    invisible(TRUE)
}

.calc_sizefactors <- function(x, nf) {
    effective_libsizes <- nf * colSums(x)
    effective_libsizes / mean(effective_libsizes)
}

.get_gene_names <- function(x) {
    if (!is.null(x@gene_names)) {
        return(x@gene_names)
    }
    if (!is.null(rownames(x@data))) {
        return(rownames(x@data))
    }
    paste0("gene_", seq_len(nrow(x@data)))
}
