make_seed_counts <- function() {
    set.seed(1)
    matrix(
        stats::rnbinom(600, mu = 50, size = 5),
        nrow = 100,
        ncol = 6
    )
}

test_that("sim_gene_counts uses lognormal fold changes by default", {
    x <- sim_gene_counts(
        n_genes = 100,
        p_DEG = c(0.05, 0.15),
        fc_mean = c(4, 8),
        seed_counts = make_seed_counts()
    )

    fc_group_1 <- x@meta$sim_params@fc[x@meta$sim_params@fc[, 1] != 1, 1]
    fc_group_2 <- x@meta$sim_params@fc[x@meta$sim_params@fc[, 2] != 1, 2]

    expect_identical(x@meta$sim_params@fc_dist, "lognormal")
    expect_equal(x@meta$sim_params@fc_sd, c(0.3, 0.3))
    expect_true(all(fc_group_1 > 0))
    expect_true(all(fc_group_2 > 0))
    expect_gt(length(unique(fc_group_1)), 1)
    expect_gt(length(unique(fc_group_2)), 1)
})

test_that("sim_gene_counts can keep fixed fold changes explicitly", {
    x <- sim_gene_counts(
        n_genes = 100,
        p_DEG = c(0.05, 0.15),
        fc_mean = c(4, 8),
        fc_dist = "fixed",
        seed_counts = make_seed_counts()
    )

    expect_identical(x@meta$sim_params@fc_dist, "fixed")
    expect_equal(unique(x@meta$sim_params@fc[x@meta$sim_params@fc[, 1] != 1, 1]), 4)
    expect_equal(unique(x@meta$sim_params@fc[x@meta$sim_params@fc[, 2] != 1, 2]), 8)
})

test_that("sim_gene_counts can randomize fold changes", {
    seed <- make_seed_counts()

    set.seed(2)
    x_lognormal <- sim_gene_counts(
        n_genes = 100,
        p_DEG = c(0.05, 0.25),
        fc_mean = c(4, 4),
        fc_dist = "lognormal",
        fc_sd = 0.5,
        seed_counts = seed
    )
    fc_lognormal <- x_lognormal@meta$sim_params@fc[, 2]
    fc_lognormal <- fc_lognormal[fc_lognormal != 1]

    expect_identical(x_lognormal@meta$sim_params@fc_dist, "lognormal")
    expect_equal(x_lognormal@meta$sim_params@fc_sd, c(0.5, 0.5))
    expect_true(all(fc_lognormal > 0))
    expect_gt(length(unique(fc_lognormal)), 1)

    set.seed(2)
    x_gamma <- sim_gene_counts(
        n_genes = 100,
        p_DEG = c(0.05, 0.25),
        fc_mean = c(4, 4),
        fc_dist = "gamma",
        fc_sd = 0.5,
        seed_counts = seed
    )
    fc_gamma <- x_gamma@meta$sim_params@fc[, 2]
    fc_gamma <- fc_gamma[fc_gamma != 1]

    expect_identical(x_gamma@meta$sim_params@fc_dist, "gamma")
    expect_true(all(fc_gamma > 0))
    expect_gt(length(unique(fc_gamma)), 1)
})

test_that("sim_gene_ts_counts returns the expected data and design", {
    x <- sim_gene_ts_counts(
        n_genes = 40,
        ts = c(0, 2, 5),
        n_replicates = c(2, 1, 3),
        cv_libsize = 0,
        seed_counts = make_seed_counts()
    )

    expect_s4_class(x, "SeqCountData")
    expect_equal(dim(x@data), c(40, 6))
    expect_identical(rownames(x@data), x@gene_names)
    expect_identical(rownames(x@exp_design), colnames(x@data))
    expect_equal(x@meta$libsize_factors, stats::setNames(rep(1, 6), colnames(x@data)))
})

test_that("sim_gene_ts_counts records temporal truth and fixed fold changes", {
    x <- sim_gene_ts_counts(
        n_genes = 40,
        ts = c(0, 1, 3, 7),
        n_replicates = 1,
        p_patterns = c(increase = 0.5, null = 0.5),
        fc_mean = 4,
        fc_dist = "fixed",
        seed_counts = make_seed_counts()
    )

    truth <- x@meta$ts_truth
    fc <- x@meta$sim_params@fc

    expect_equal(truth$p_patterns[c("increase", "null")], c(increase = 0.5, null = 0.5))
    expect_equal(truth$realized_p_patterns[c("increase", "null")], c(increase = 0.5, null = 0.5))
    dynamic_range <- apply(fc, 1, max) / apply(fc, 1, min)

    expect_equal(truth$pattern_counts[c("increase", "null")], c(increase = 20L, null = 20L))
    expect_true(all(truth$pattern %in% c("increase", "null")))
    expect_equal(unname(truth$dynamic_fc[truth$pattern == "increase"]), rep(4, 20))
    expect_equal(unname(truth$dynamic_fc[truth$pattern == "null"]), rep(1, 20))
    expect_equal(unname(dynamic_range[truth$pattern == "increase"]), rep(4, 20))
    expect_equal(unname(dynamic_range[truth$pattern == "null"]), rep(1, 20))
})

test_that("sim_gene_ts_counts validates time-series arguments", {
    seed <- make_seed_counts()

    expect_error(
        sim_gene_ts_counts(ts = c(0, 2, 1), seed_counts = seed),
        "strictly increasing"
    )
    expect_error(
        sim_gene_ts_counts(
            ts = 0:2,
            n_replicates = c(2, 2),
            seed_counts = seed
        ),
        "same length as `ts`"
    )
    expect_error(
        sim_gene_ts_counts(
            p_patterns = c(increase = 1, unknown = 1),
            seed_counts = seed
        ),
        "Unknown temporal pattern"
    )
})
