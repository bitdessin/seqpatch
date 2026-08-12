.set_full_model <- function(exp_design = NULL) {
    if (is.null(exp_design)) stop("exp_design should not be NULL.", .call = FALSE)

    design_terms <- names(exp_design)
    if (length(design_terms) == 0L) {
        design <- stats::model.matrix(~ 1, data = exp_design)
    } else {
        design_terms <- vapply(
            design_terms,
            function(term) {
                if (grepl("^[A-Za-z0-9_.]+$", term)) {
                    term
                } else {
                    sprintf("`%s`", term)
                }
            },
            FUN.VALUE = character(1)
        )
        design_formula <- stats::as.formula(paste("~", paste(design_terms, collapse = "*")))
        design <- stats::model.matrix(design_formula, data = exp_design)
    }
    design
}

#' Run edgeR to Detect DEGs
#'
#' Runs quasi-likelihood F tests to detect DEGs
#' using edgeR. By default, this function tests for differences across the
#' groups defined by the `group` column of `exp_design`.
#'
#' @param x A count matrix.
#' @param exp_design A data frame of experimental design. Must contain `group`.
#' @param nf A vector of normalization factors.
#' @param design Optional model matrix. Defaults to `model.matrix(~ group,
#'   data = exp_design)`.
#' @param coef Coefficient or coefficients passed to `edgeR::glmQLFTest()`.
#' @param contrast Optional contrast passed to `edgeR::glmQLFTest()`.
#' @param ... Reserved for future extension.
#' @return A data frame with `p_value` and `q_value` columns.
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
#'     x@meta$DEGES$params
#' }
#'
#' @export
test_edger <- function(
    x,
    exp_design,
    nf = NULL,
    design = NULL,
    coef = NULL,
    contrast = NULL,
    ...
) {
    if (!requireNamespace("edgeR", quietly = TRUE)) {
        stop("Package `edgeR` is required for this function.", call. = FALSE)
    }

    # set default values
    if (is.null(nf)) nf <- rep(1, ncol(x))
    if (is.null(design)) design <- .set_full_model(exp_design)
    if (is.null(coef) && is.null(contrast)) coef <- if (ncol(design) > 1L) seq.int(2L, ncol(design)) else 1L

    x <- as.matrix(x)

    d <- edgeR::DGEList(counts = x, group = exp_design$group)
    d$samples$norm.factors <- nf
    d <- edgeR::estimateDisp(d, design)
    fit <- edgeR::glmQLFit(d, design)
    if (is.null(contrast)) {
        test_output <- edgeR::glmQLFTest(fit, coef = coef)
    } else {
        test_output <- edgeR::glmQLFTest(fit, contrast = contrast)
    }

    test_output <- as.data.frame(edgeR::topTags(test_output, n = nrow(x), sort.by = "none")$table)
    p <- test_output$PValue
    p[is.na(p)] <- 1
    q <- stats::p.adjust(p, method = "BH")

    data.frame(
        p_value = p,
        q_value = q,
        row.names = NULL,
        check.names = FALSE
    )
}

#' Run DESeq2 to Detect DEGs
#'
#' Runs a DESeq2 LRT to detect DEGs. By default, this function tests
#' for differences across the groups
#' defined by the `group` column of `exp_design`.
#'
#' @param x A count matrix.
#' @param exp_design A data frame of experimental design. Must contain `group`.
#' @param nf A vector of normalization factors.
#' @param design Optional DESeq2 design formula. Defaults to `~ group`.
#' @param full Optional full model for an LRT. Defaults to
#'   `design`.
#' @param reduced Optional reduced model for an LRT. Defaults
#'   to `~ 1`.
#' @param contrast Optional contrast passed to `DESeq2::results()`.
#' @param ... Reserved for future extension.
#'
#' @return A data frame with at least `p_value` and `q_value`
#'   columns.
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
#'     x@meta$DEGES$params
#' }
#'
#' @export
test_deseq2 <- function(
    x,
    exp_design,
    nf = NULL,
    design = NULL,
    full = NULL,
    reduced = NULL,
    contrast = NULL,
    ...) {
    if (!requireNamespace("DESeq2", quietly = TRUE)) {
        stop("Package `DESeq2` is required for this function.", call. = FALSE)
    }

    if (is.null(nf)) nf <- rep(1, ncol(x))
    if (is.null(design)) design <- ~ group
    if (is.null(full)) full <- design
    if (is.null(reduced)) reduced <- ~ 1

    x <- as.matrix(x)
    exp_design <- as.data.frame(exp_design)
    rownames(exp_design) <- colnames(x)
    effective_lib_sizes <- nf * colSums(x)

    d <- DESeq2::DESeqDataSetFromMatrix(
        countData = x,
        colData = exp_design,
        design = design
    )
    DESeq2::sizeFactors(d) <- effective_lib_sizes / mean(effective_lib_sizes)
    d <- DESeq2::estimateDispersions(d)
    d <- DESeq2::nbinomLRT(d, full = full, reduced = reduced)
    if (is.null(contrast)) {
        test_output <- DESeq2::results(d)
    } else {
        test_output <- DESeq2::results(d, contrast = contrast)
    }

    test_output <- as.data.frame(test_output)
    p <- test_output$pvalue
    p[is.na(p)] <- 1
    q <- stats::p.adjust(p, method = "BH")

    data.frame(
        p_value = p,
        q_value = q,
        row.names = NULL,
        check.names = FALSE
    )
}
