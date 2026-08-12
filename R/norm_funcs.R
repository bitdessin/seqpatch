#' Estimate TMM Normalization Factors
#'
#' @param x A count matrix.
#' @param exp_design Sample metadata as a data frame. Must contain `group`.
#' @param ... Additional arguments passed to `edgeR::calcNormFactors()`.
#' @return Numeric normalization factors, one per sample.
#'
#' @examples
#' if (requireNamespace("edgeR", quietly = TRUE)) {
#'     set.seed(1)
#'     counts <- matrix(stats::rnbinom(600, mu = 50, size = 5), nrow = 100)
#'     colnames(counts) <- paste0("sample_", seq_len(ncol(counts)))
#'     rownames(counts) <- paste0("gene_", seq_len(nrow(counts)))
#'     exp_design <- data.frame(group = rep(c("A", "B"), each = 3))
#'     x <- newSeqCountData(counts, exp_design = exp_design)
#'     x <- calc_nf(x, norm_func = norm_tmm, test_func = test_edger, iter = 1)
#'     x@meta$nf
#' }
#'
#' @export
norm_tmm <- function(x, exp_design, ...) {
    if (!requireNamespace("edgeR", quietly = TRUE)) {
        stop("Package `edgeR` is required for this function.", call. = FALSE)
    }

    d <- edgeR::DGEList(counts = x, group = exp_design$group)
    d <- edgeR::calcNormFactors(d, ...)
    nf <- d$samples$norm.factors
    names(nf) <- colnames(x)
    nf
}


#' Estimate RLE Normalization Factors
#'
#' Uses DESeq2 size factor estimation, then converts size factors to
#' normalization factors by dividing by raw library sizes.
#'
#' @param x A count matrix.
#' @param exp_design Sample metadata as a data frame. Must contain `group`.
#' @param ... Reserved for future extension.
#'
#' @return Numeric normalization factors, one per sample.
#'
#' @examples
#' if (requireNamespace("DESeq2", quietly = TRUE)) {
#'     set.seed(1)
#'     counts <- matrix(stats::rnbinom(600, mu = 50, size = 5), nrow = 100)
#'     colnames(counts) <- paste0("sample_", seq_len(ncol(counts)))
#'     rownames(counts) <- paste0("gene_", seq_len(nrow(counts)))
#'     exp_design <- data.frame(group = rep(c("A", "B"), each = 3))
#'     x <- newSeqCountData(counts, exp_design = exp_design)
#'     x <- calc_nf(x, norm_func = norm_rle, test_func = test_deseq2, iter = 1)
#'     x@meta$nf
#' }
#'
#' @export
norm_rle <- function(x, exp_design, ...) {
    if (!requireNamespace("DESeq2", quietly = TRUE)) {
        stop("Package `DESeq2` is required for this function.", call. = FALSE)
    }
    exp_design <- as.data.frame(exp_design)
    rownames(exp_design) <- colnames(x)

    d <- DESeq2::DESeqDataSetFromMatrix(
        countData = x,
        colData = exp_design,
        design = ~ 1
    )
    d <- DESeq2::estimateSizeFactors(d)
    nf <- DESeq2::sizeFactors(d) / colSums(x)
    names(nf) <- colnames(x)
    nf
}
