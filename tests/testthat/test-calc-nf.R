make_count_matrix <- function(n_genes = 100L, n_samples = 6L) {
    set.seed(1)
    counts <- matrix(
        stats::rnbinom(n_genes * n_samples, mu = 50, size = 10),
        nrow = n_genes,
        ncol = n_samples
    )
    rownames(counts) <- paste0("gene_", seq_len(n_genes))
    colnames(counts) <- paste0("sample_", seq_len(n_samples))
    counts
}

load_package <- function(package) {
    suppressWarnings(
        suppressPackageStartupMessages(
            library(package, character.only = TRUE)
        )
    )
}

test_deseq2_wald <- function(x, exp_design, nf = NULL, ...) {
    if (is.null(nf)) nf <- rep(1, ncol(x))

    x <- as.matrix(x)
    exp_design <- as.data.frame(exp_design)
    rownames(exp_design) <- colnames(x)
    effective_lib_sizes <- nf * colSums(x)

    dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = x,
        colData = exp_design,
        design = ~ group
    )
    DESeq2::sizeFactors(dds) <- effective_lib_sizes / mean(effective_lib_sizes)
    dds <- DESeq2::estimateDispersions(dds)
    dds <- DESeq2::nbinomWaldTest(dds)

    test_output <- as.data.frame(DESeq2::results(dds))
    p <- test_output$pvalue
    p[is.na(p)] <- 1

    data.frame(
        p_value = p,
        q_value = stats::p.adjust(p, method = "BH"),
        row.names = NULL,
        check.names = FALSE
    )
}

expect_valid_calc_nf_result <- function(x, iter) {
    expect_true(methods::is(x, "SeqCountData"))
    expect_type(x@meta$nf, "double")
    expect_length(x@meta$nf, ncol(x@data))
    expect_false(anyNA(x@meta$nf))
    expect_true(is.list(x@meta$DEGES))
    expect_equal(x@meta$DEGES$params$n_iters, iter)
    expect_true(is.matrix(x@meta$DEGES$potential_DEGs))
    expect_identical(dim(x@meta$DEGES$potential_DEGs), c(nrow(x@data), as.integer(iter)))
}

expect_norm_factors_equal <- function(observed, expected) {
    expect_equal(
        unname(observed),
        unname(expected),
        tolerance = 1e-8
    )
}

test_that("calc_nf runs with TMM/edgeR and RLE/DESeq2 pipelines", {
    skip_if_not_installed("edgeR")
    skip_if_not_installed("DESeq2")

    counts <- make_count_matrix()
    exp_design <- data.frame(group = factor(rep(c("A", "B"), each = 3)))

    x_tmm <- newSeqCountData(counts, exp_design = exp_design)
    x_tmm <- calc_nf(
        x_tmm,
        norm_func = norm_tmm,
        test_func = test_edger,
        iter = 1
    )
    expect_valid_calc_nf_result(x_tmm, iter = 1L)

    x_rle <- newSeqCountData(counts, exp_design = exp_design)
    x_rle <- suppressMessages(
        calc_nf(
            x_rle,
            norm_func = norm_rle,
            test_func = test_deseq2,
            iter = 1
        )
    )
    expect_valid_calc_nf_result(x_rle, iter = 1L)
})

test_that("two-group TMM/edgeR factors match TCC", {
    skip_if_not_installed("TCC")
    skip_if_not_installed("edgeR")
    load_package("edgeR")
    load_package("TCC")

    data("hypoData", package = "TCC")
    group <- factor(c(1, 1, 1, 2, 2, 2))
    exp_design <- data.frame(group = group)

    tcc <- methods::new("TCC", hypoData, group)
    tcc <- suppressMessages(
        TCC::calcNormFactors(
            tcc,
            norm.method = "tmm",
            test.method = "edger",
            iteration = 3,
            FDR = 0.1,
            floorPDEG = 0.05
        )
    )

    x <- newSeqCountData(hypoData, exp_design = exp_design)
    x <- calc_nf(
        x,
        norm_func = norm_tmm,
        test_func = test_edger,
        iter = 3,
        p_cutoff = 0,
        q_cutoff = 0.1,
        min_rm_genes = 0.05
    )

    expect_norm_factors_equal(x@meta$nf, tcc$norm.factors)
})

test_that("two-group RLE/DESeq2 factors match TCC with Wald screening", {
    skip_if_not_installed("TCC")
    skip_if_not_installed("DESeq2")
    load_package("DESeq2")
    load_package("TCC")

    data("hypoData", package = "TCC")
    group <- factor(c(1, 1, 1, 2, 2, 2))
    exp_design <- data.frame(group = group)

    tcc <- methods::new("TCC", hypoData, group)
    tcc <- suppressMessages(
        TCC::calcNormFactors(
            tcc,
            norm.method = "deseq2",
            test.method = "deseq2",
            iteration = 3,
            FDR = 0.1,
            floorPDEG = 0.05
        )
    )

    x <- newSeqCountData(hypoData, exp_design = exp_design)
    x <- suppressMessages(
        calc_nf(
            x,
            norm_func = norm_rle,
            test_func = test_deseq2_wald,
            iter = 3,
            p_cutoff = 0,
            q_cutoff = 0.1,
            min_rm_genes = 0.05
        )
    )

    expect_norm_factors_equal(x@meta$nf, tcc$norm.factors)
})

test_that("multi-group TMM/edgeR factors match TCC", {
    skip_if_not_installed("TCC")
    skip_if_not_installed("edgeR")
    load_package("edgeR")
    load_package("TCC")

    data("hypoData_mg", package = "TCC")
    group <- factor(c(1, 1, 1, 2, 2, 2, 3, 3, 3))
    exp_design <- data.frame(group = group)

    tcc <- methods::new("TCC", hypoData_mg, group)
    tcc <- suppressMessages(
        TCC::calcNormFactors(
            tcc,
            norm.method = "tmm",
            test.method = "edger",
            iteration = 3,
            FDR = 0.1,
            floorPDEG = 0.05
        )
    )

    x <- newSeqCountData(hypoData_mg, exp_design = exp_design)
    x <- calc_nf(
        x,
        norm_func = norm_tmm,
        test_func = test_edger,
        iter = 3,
        p_cutoff = 0,
        q_cutoff = 0.1,
        min_rm_genes = 0.05
    )

    expect_norm_factors_equal(x@meta$nf, tcc$norm.factors)
})

test_that("multi-group RLE/DESeq2 factors match TCC", {
    skip_if_not_installed("TCC")
    skip_if_not_installed("DESeq2")
    load_package("DESeq2")
    load_package("TCC")

    data("hypoData_mg", package = "TCC")
    group <- factor(c(1, 1, 1, 2, 2, 2, 3, 3, 3))
    exp_design <- data.frame(group = group)

    tcc <- methods::new("TCC", hypoData_mg, group)
    tcc <- suppressMessages(
        TCC::calcNormFactors(
            tcc,
            norm.method = "deseq2",
            test.method = "deseq2",
            iteration = 3,
            FDR = 0.1,
            floorPDEG = 0.05
        )
    )

    x <- newSeqCountData(hypoData_mg, exp_design = exp_design)
    x <- suppressMessages(
        calc_nf(
            x,
            norm_func = norm_rle,
            test_func = test_deseq2,
            iter = 3,
            p_cutoff = 0,
            q_cutoff = 0.1,
            min_rm_genes = 0.05
        )
    )

    expect_norm_factors_equal(x@meta$nf, tcc$norm.factors)
})
