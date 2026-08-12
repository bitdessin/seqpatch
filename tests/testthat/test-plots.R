make_profile_data <- function() {
    counts <- matrix(
        c(
            10, 12, 20, 22, 30, 32,
            8, 10, 16, 18, 24, 26,
            30, 28, 20, 18, 10, 8,
            24, 22, 16, 14, 8, 6
        ),
        nrow = 4,
        byrow = TRUE,
        dimnames = list(
            paste0("gene_", 1:4),
            paste0("sample_", 1:6)
        )
    )
    exp_design <- data.frame(
        group = factor(rep(c("A", "B", "C"), each = 2))
    )
    newSeqCountData(
        counts,
        exp_design = exp_design,
        norm_factors = rep(1, ncol(counts))
    )
}

test_that("plot_gene_profiles draws gene and cluster mean profiles", {
    x <- make_profile_data()
    p <- plot_gene_profiles(
        x,
        cl = c("up", "up", "down", "down"),
        by = "group"
    )

    expect_s3_class(p, "ggplot")
    expect_equal(nrow(p$data), 12)
    expect_equal(levels(p$data$sample), c("A", "B", "C"))
    expect_equal(levels(p$data$cluster), c("up (n = 2)", "down (n = 2)"))
    expect_equal(length(p$layers), 2)
    expect_equal(nrow(p$layers[[2]]$data), 6)
    expect_no_error(ggplot2::ggplot_build(p))

    gene_means <- tapply(p$data$expression, p$data$gene, mean)
    gene_sds <- tapply(p$data$expression, p$data$gene, stats::sd)
    expect_equal(as.vector(gene_means), rep(0, 4), tolerance = 1e-12)
    expect_equal(as.vector(gene_sds), rep(1, 4), tolerance = 1e-12)
})

test_that("plot_gene_profiles can plot unscaled sample profiles", {
    x <- make_profile_data()
    p <- plot_gene_profiles(
        x,
        cl = as.integer(c(1, 1, 2, 2)),
        by = NULL,
        zscale = FALSE,
        log_transform = FALSE
    )

    expect_equal(nrow(p$data), nrow(x@data) * ncol(x@data))
    expect_equal(levels(p$data$sample), colnames(x@data))
    expect_equal(p$labels$y, "Normalized expression")
    expect_no_error(ggplot2::ggplot_build(p))
})
