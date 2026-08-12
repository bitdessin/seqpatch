#' Simulate RNA-seq Gene Counts
#'
#' Generates a synthetic RNA-seq count matrix from a seed count matrix.
#' If no seed matrix is supplied, the `arab` count matrix bundled with
#' seqpatch is used as the seed. The seed matrix is scaled to a common library
#' size and then used as the sampling population by computing row means and
#' variances. Counts are sampled from a negative binomial distribution after
#' applying group-level fold changes and optional sample-specific library-size
#' factors.
#'
#' @param n_genes Number of genes to simulate.
#' @param n_replicates Numeric vector giving the number of replicates per group.
#' @param p_DEG Numeric vector giving the proportion of genes assigned as DEGs
#'   for each group.
#' @param fc_mean Numeric vector giving the mean fold change for genes assigned
#'   to each group. When `fc_dist = "fixed"`, this exact value is used.
#' @param fc_dist Distribution used to generate fold changes for genes assigned
#'   as DEGs. One of `"lognormal"`, `"gamma"`, or `"fixed"`. The default is
#'   `"lognormal"`.
#' @param fc_sd Spread parameter for randomized fold changes. For
#'   `fc_dist = "lognormal"`, this is the standard deviation on the natural log
#'   fold change scale. For `fc_dist = "gamma"`, this is the coefficient of
#'   variation.
#' @param group_names Optional group names. If `NULL`, groups are named
#'   `"group_1"`, `"group_2"`, and so on.
#' @param cv_libsize Coefficient of variation of sample-specific library-size
#'   factors. Factors are drawn from a mean-one log-normal distribution and
#'   rescaled to have arithmetic mean one. Set to 0 to disable library-size
#'   variation. The default is 0, preserving the previous behavior.
#' @param seed_counts Optional seed count matrix. Rows are genes and columns
#'   are samples. If `NULL`, the package dataset `arab` is used.
#'
#' @return A \linkS4class{SeqCountData} object. Simulation parameters are stored
#'   in `x@meta$sim_params`, and realized sample-specific library-size factors
#'   are stored in `x@meta$libsize_factors`.
#'
#' @examples
#' x <- sim_gene_counts()
#' table(def_DEG(x, fc = 4))
#'
#' @export
sim_gene_counts <- function(n_genes = 10000,
                            n_replicates = c(3, 3),
                            p_DEG = c(0.04, 0.01),
                            fc_mean = c(4, 4),
                            fc_dist = c("lognormal", "gamma", "fixed"),
                            fc_sd = 0.3,
                            group_names = NULL,
                            cv_libsize = 0.2,
                            seed_counts = NULL) {
    n_genes <- .sim_validate_n_genes(n_genes)
    n_replicates <- .sim_validate_replicates(n_replicates)
    n_groups <- length(n_replicates)
    fc_dist <- match.arg(fc_dist)

    p_DEG <- .sim_recycle_numeric(p_DEG, n_groups, "p_DEG")
    fc_mean <- .sim_recycle_numeric(fc_mean, n_groups, "fc_mean")
    fc_sd <- .sim_recycle_numeric(fc_sd, n_groups, "fc_sd")

    if (any(p_DEG < 0 | p_DEG > 1)) {
        stop("`p_DEG` must contain values between 0 and 1.", call. = FALSE)
    }
    if (any(fc_mean <= 0)) {
        stop("`fc_mean` must contain positive values.", call. = FALSE)
    }

    if (is.null(group_names)) {
        group_names <- paste0("group_", seq_len(n_groups))
    } else {
        group_names <- as.character(group_names)
        if (length(group_names) != n_groups || anyNA(group_names) || any(group_names == "")) {
            stop("`group_names` must contain one non-empty name per group.", call. = FALSE)
        }
        if (anyDuplicated(group_names)) {
            stop("`group_names` must be unique.", call. = FALSE)
        }
    }

    population <- .sim_prepare_population(seed_counts)
    gene_names <- paste0("gene_", seq_len(n_genes))
    params <- .sim_sample_gene_params(population, gene_names)
    layout <- .sim_group_layout(n_replicates, group_names)

    fc <- .sim_fc_matrix(n_genes, p_DEG, fc_mean, fc_dist, fc_sd, group_names)
    rownames(fc) <- gene_names
    sample_fc <- .sim_sample_fc_matrix(fc, layout$group_idx, layout$sample_names)

    libsize_factors <- .sim_libsize_factors(length(layout$sample_names), cv_libsize)
    names(libsize_factors) <- layout$sample_names

    data <- .sim_nb_count_matrix(params, sample_fc, libsize_factors, layout$sample_names)
    exp_design <- data.frame(
        group = factor(group_names[layout$group_idx], levels = group_names),
        row.names = layout$sample_names
    )

    sim_params <- new(
        "SimParams",
        n_genes = as.numeric(n_genes),
        n_groups = as.numeric(n_groups),
        n_replicates = as.numeric(n_replicates),
        p_DEG = as.numeric(p_DEG),
        fc_dist = fc_dist,
        fc_sd = as.numeric(fc_sd),
        fc = fc,
        params_population = population,
        params = params
    )

    new(
        "SeqCountData",
        data = data,
        gene_names = gene_names,
        exp_design = exp_design,
        meta = list(
            sim_params = sim_params,
            libsize_factors = libsize_factors,
            cv_libsize = as.numeric(cv_libsize)
        )
    )
}


#' Simulate Time-series RNA-seq Gene Counts
#'
#' Generates synthetic time-series RNA-seq counts with known latent temporal
#' patterns. Gene-specific baseline abundance and negative-binomial dispersion
#' are sampled from the same seed-derived population used by
#' [sim_gene_counts()]. Each simulated gene is assigned to one of eight latent
#' temporal patterns, and the expected expression trajectory is generated as a
#' fold-change profile around the sampled baseline abundance. Sample-specific
#' library-size variation can be added independently of the temporal pattern.
#'
#' The eight built-in patterns are `"increase"`, `"decrease"`, `"early_peak"`,
#' `"late_peak"`, `"biphasic"`, `"oscillatory"`, `"delayed_increase"`, and
#' `"null"`. Pattern functions are evaluated on the supplied numeric sampling
#' times after scaling them to the interval from 0 to 1, so unevenly spaced
#' time points are represented explicitly.
#'
#' @param n_genes Number of genes to simulate.
#' @param ts Numeric vector of ordered sampling times.
#' @param n_replicates Number of biological replicates per time point. Either a
#'   single value or a vector with one value per time point.
#' @param p_patterns Relative proportions of the eight temporal patterns.
#'   If `NULL`, all patterns are equally represented. An unnamed numeric vector
#'   must have length eight and follows the pattern order given above. A named
#'   vector may specify any subset of the eight pattern names; omitted patterns
#'   receive probability zero. Values are normalized to sum to one.
#' @param fc_mean Mean dynamic-range fold change for non-null genes. The
#'   generated temporal profile is scaled so that, before count sampling, the
#'   ratio between its largest and smallest expected fold-change multipliers is
#'   the gene-specific value sampled from `fc_dist`.
#' @param fc_dist Distribution used to generate gene-specific dynamic-range fold
#'   changes. One of `"lognormal"`, `"gamma"`, or `"fixed"`.
#' @param fc_sd Spread parameter for randomized dynamic-range fold changes. Its
#'   interpretation is the same as in [sim_gene_counts()].
#' @param cv_libsize Coefficient of variation of sample-specific
#'   library-size factors. Factors are drawn from a mean-one log-normal
#'   distribution and rescaled to have arithmetic mean one. Set to 0 to disable
#'   library-size variation.
#' @param seed_counts Optional seed count matrix. Rows are genes and columns are
#'   samples. If `NULL`, the package dataset `arab` is used.
#'
#' @return A \linkS4class{SeqCountData} object. The `meta` slot contains
#'   `sim_params`, `libsize_factors`, and `ts_truth`. The `ts_truth` entry
#'   stores the true temporal pattern assigned to each gene, requested and
#'   realized pattern proportions, temporal pattern templates, gene-specific
#'   dynamic-range fold changes, and sampling times.
#'
#' @examples
#' x <- sim_gene_ts_counts()
#' table(x@meta$ts_truth$pattern)
#'
#' @export
sim_gene_ts_counts <- function(n_genes = 10000,
                               ts = 0:7,
                               n_replicates = 3,
                               p_patterns = NULL,
                               fc_mean = 4,
                               fc_dist = c("lognormal", "gamma", "fixed"),
                               fc_sd = 0.3,
                               cv_libsize = 0.2,
                               seed_counts = NULL) {
    n_genes <- .sim_validate_n_genes(n_genes)
    ts <- .sim_validate_ts(ts)
    n_time <- length(ts)
    n_replicates <- .sim_validate_replicates(n_replicates)
    if (length(n_replicates) == 1L) {
        n_replicates <- rep(n_replicates, n_time)
    } else if (length(n_replicates) != n_time) {
        stop(
            "`n_replicates` must have length 1 or the same length as `ts`.",
            call. = FALSE
        )
    }

    fc_dist <- match.arg(fc_dist)
    fc_mean <- as.numeric(fc_mean)
    fc_sd <- as.numeric(fc_sd)
    if (length(fc_mean) != 1L || !is.finite(fc_mean) || fc_mean < 1) {
        stop("`fc_mean` must be a single finite value greater than or equal to 1.", call. = FALSE)
    }
    if (length(fc_sd) != 1L || !is.finite(fc_sd) || fc_sd < 0) {
        stop("`fc_sd` must be a single non-negative finite value.", call. = FALSE)
    }

    probs <- .sim_ts_p_patterns(p_patterns)
    patterns <- .sim_ts_assign_patterns(n_genes, probs)
    templates <- .sim_ts_pattern_matrix(ts)

    population <- .sim_prepare_population(seed_counts)
    gene_names <- paste0("gene_", seq_len(n_genes))
    params <- .sim_sample_gene_params(population, gene_names)

    group_names <- .sim_time_group_names(ts)
    layout <- .sim_group_layout(n_replicates, group_names)

    dynamic_fc <- .sim_fc_values(n_genes, fc_mean, fc_dist, fc_sd)
    dynamic_fc <- pmax(dynamic_fc, 1)
    dynamic_fc[patterns == "null"] <- 1
    names(dynamic_fc) <- gene_names

    ts_fc <- .sim_ts_fc_matrix(patterns, templates, dynamic_fc, gene_names, group_names)
    sample_fc <- .sim_sample_fc_matrix(ts_fc, layout$group_idx, layout$sample_names)

    libsize_factors <- .sim_libsize_factors(length(layout$sample_names), cv_libsize)
    names(libsize_factors) <- layout$sample_names

    data <- .sim_nb_count_matrix(params, sample_fc, libsize_factors, layout$sample_names)

    exp_design <- data.frame(
        group = factor(group_names[layout$group_idx], levels = group_names),
        #time = ts[layout$group_idx],
        #replicate = layout$replicate_idx,
        row.names = layout$sample_names
    )

    p_dynamic <- colMeans(abs(log2(ts_fc)) > sqrt(.Machine$double.eps))
    sim_params <- new(
        "SimParams",
        n_genes = as.numeric(n_genes),
        n_groups = as.numeric(n_time),
        n_replicates = as.numeric(n_replicates),
        p_DEG = as.numeric(p_dynamic),
        fc_dist = fc_dist,
        fc_sd = rep(as.numeric(fc_sd), n_time),
        fc = ts_fc,
        params_population = population,
        params = params
    )

    realized_counts <- table(factor(patterns, levels = names(probs)))
    realized_probs <- as.numeric(realized_counts) / n_genes
    names(realized_probs) <- names(probs)

    new(
        "SeqCountData",
        data = data,
        gene_names = gene_names,
        exp_design = exp_design,
        meta = list(
            sim_params = sim_params,
            libsize_factors = libsize_factors,
            cv_libsize = as.numeric(cv_libsize),
            ts_truth = list(
                ts = ts,
                pattern = stats::setNames(patterns, gene_names),
                p_patterns = probs,
                pattern_counts = stats::setNames(as.integer(realized_counts), names(probs)),
                realized_p_patterns = realized_probs,
                pattern_templates = templates,
                dynamic_fc = dynamic_fc
            )
        )
    )
}


#' Define Simulated DEGs
#'
#' Uses the fold change matrix stored in `x@meta$sim_params@fc` to define
#' simulated DEG status. A gene is marked `TRUE` when the largest between-group fold
#' change is at least `fc`.
#'
#' @param x A \linkS4class{SeqCountData} object whose `meta` slot contains a
#'   `sim_params` entry of class \linkS4class{SimParams}.
#' @param fc Fold change cutoff. Defaults to 2.
#'
#' @return A logical vector with one value per gene.
#'
#' @examples
#' x <- sim_gene_counts()
#' head(def_DEG(x, fc = 4))
#'
#' @export
def_DEG <- function(x, fc = 2) {
    fc_mx <- x@meta$sim_params@fc
    max_fc <- apply(fc_mx, 1L, max)
    min_fc <- apply(fc_mx, 1L, min)
    is_deg <- (max_fc / min_fc) >= fc
    names(is_deg) <- rownames(fc_mx)
    is_deg
}


.sim_validate_n_genes <- function(n_genes) {
    n_genes <- as.integer(n_genes)
    if (length(n_genes) != 1L || is.na(n_genes) || n_genes < 1L) {
        stop("`n_genes` must be a positive integer.", call. = FALSE)
    }
    n_genes
}

.sim_validate_replicates <- function(n_replicates) {
    n_replicates <- as.integer(n_replicates)
    if (length(n_replicates) < 1L || anyNA(n_replicates) || any(n_replicates < 1L)) {
        stop("`n_replicates` must contain positive integers.", call. = FALSE)
    }
    n_replicates
}

.sim_recycle_numeric <- function(x, n, arg) {
    x <- as.numeric(x)
    if (length(x) == 1L) {
        x <- rep(x, n)
    }
    if (length(x) != n || any(!is.finite(x))) {
        stop(sprintf("`%s` must have length 1 or %d and contain finite values.", arg, n), call. = FALSE)
    }
    x
}

.sim_validate_ts <- function(ts) {
    ts <- as.numeric(ts)
    if (length(ts) < 2L || any(!is.finite(ts))) {
        stop("`ts` must contain at least two finite numeric values.", call. = FALSE)
    }
    if (is.unsorted(ts, strictly = TRUE)) {
        stop("`ts` must be strictly increasing.", call. = FALSE)
    }
    ts
}

.sim_init_seed_counts <- function(seed_counts) {
    if (is.null(seed_counts)) {
        data_env <- new.env(parent = emptyenv())
        utils::data("arab", package = "seqpatch", envir = data_env)
        seed_counts <- data_env$arab
    }

    seed_counts <- as.matrix(seed_counts)
    storage.mode(seed_counts) <- "numeric"
    if (nrow(seed_counts) < 1L || ncol(seed_counts) < 2L) {
        stop("`seed_counts` must contain at least one gene and two samples.", call. = FALSE)
    }
    if (any(!is.finite(seed_counts)) || any(seed_counts < 0)) {
        stop("`seed_counts` must contain finite, non-negative values.", call. = FALSE)
    }
    seed_counts
}

.sim_prepare_population <- function(seed_counts) {
    seed_counts <- .sim_init_seed_counts(seed_counts)
    population <- .sim_seed_population(seed_counts)
    if (nrow(population) < 1L) {
        stop(
            "No seed genes with finite positive mean and negative-binomial dispersion were available.",
            call. = FALSE
        )
    }
    population
}

.sim_seed_population <- function(seed_counts) {
    lib_size <- colSums(seed_counts)
    keep_samples <- is.finite(lib_size) & lib_size > 0
    seed_counts <- seed_counts[, keep_samples, drop = FALSE]
    lib_size <- lib_size[keep_samples]
    if (length(lib_size) < 2L) {
        stop("At least two seed samples with positive library sizes are required.", call. = FALSE)
    }

    seed_counts <- sweep(
        seed_counts,
        2L,
        stats::median(lib_size) / lib_size,
        "*"
    )

    mean_ab <- apply(seed_counts, 1L, mean)
    var_ab <- apply(seed_counts, 1L, stats::var)
    dispersion <- (var_ab - mean_ab) / (mean_ab * mean_ab)
    population <- data.frame(
        mean = mean_ab,
        var = var_ab,
        dispersion = dispersion,
        distribution = "NB"
    )
    keep <- is.finite(population$mean) &
        is.finite(population$var) &
        is.finite(population$dispersion) &
        population$mean > 0 &
        population$dispersion > 0
    population <- population[keep, , drop = FALSE]

    rownames(population) <- paste0("seed_", seq_len(nrow(population)))
    population
}

.sim_sample_gene_params <- function(population, gene_names) {
    sample_idx <- sample(seq_len(nrow(population)), length(gene_names), replace = TRUE)
    params <- population[sample_idx, , drop = FALSE]
    rownames(params) <- gene_names
    params
}

.sim_group_layout <- function(n_replicates, group_names) {
    group_idx <- rep(seq_along(group_names), times = n_replicates)
    replicate_idx <- sequence(n_replicates)
    sample_names <- paste0(group_names[group_idx], "_rep", replicate_idx)
    list(
        group_idx = group_idx,
        replicate_idx = replicate_idx,
        sample_names = sample_names
    )
}

.sim_libsize_factors <- function(n_samples, cv = 0) {
    cv <- as.numeric(cv)
    if (length(cv) != 1L || !is.finite(cv) || cv < 0) {
        stop("`cv_libsize` must be a single non-negative finite value.", call. = FALSE)
    }
    if (cv == 0 || n_samples == 1L) {
        return(rep(1, n_samples))
    }

    sigma <- sqrt(log1p(cv * cv))
    factors <- stats::rlnorm(
        n_samples,
        meanlog = -0.5 * sigma * sigma,
        sdlog = sigma
    )
    factors / mean(factors)
}

.sim_nb_count_matrix <- function(params,
                                 sample_multiplier,
                                 libsize_factors,
                                 sample_names = colnames(sample_multiplier)) {
    sample_multiplier <- as.matrix(sample_multiplier)
    n_genes <- nrow(params)
    n_samples <- ncol(sample_multiplier)

    if (nrow(sample_multiplier) != n_genes) {
        stop("`sample_multiplier` must have one row per simulated gene.", call. = FALSE)
    }
    if (length(libsize_factors) != n_samples) {
        stop("`libsize_factors` must have one value per simulated sample.", call. = FALSE)
    }

    mu <- sweep(sample_multiplier, 1L, params$mean, "*")
    mu <- sweep(mu, 2L, libsize_factors, "*")

    data <- vapply(
        seq_len(n_samples),
        function(j) {
            stats::rnbinom(
                n = n_genes,
                mu = mu[, j],
                size = 1 / params$dispersion
            )
        },
        numeric(n_genes)
    )

    matrix(
        data,
        nrow = n_genes,
        ncol = n_samples,
        dimnames = list(rownames(params), sample_names)
    )
}

.sim_fc_matrix <- function(n_genes, p_DEG, fc_mean, fc_dist, fc_sd, group_names) {
    n_groups <- length(group_names)
    fc <- matrix(
        1,
        nrow = n_genes,
        ncol = n_groups,
        dimnames = list(NULL, group_names)
    )
    n_deg <- as.integer(round(n_genes * p_DEG))
    while (sum(n_deg) > n_genes) {
        i <- which.max(n_deg)
        n_deg[i] <- n_deg[i] - 1L
    }

    idx_start <- 1L
    for (i in seq_len(n_groups)) {
        idx_end <- idx_start + n_deg[i] - 1L
        if (idx_end >= idx_start) {
            fc[idx_start:idx_end, i] <- .sim_fc_values(
                n = idx_end - idx_start + 1L,
                mean_fc = fc_mean[i],
                fc_dist = fc_dist,
                fc_sd = fc_sd[i]
            )
            idx_start <- idx_end + 1L
        }
    }
    fc
}

.sim_fc_values <- function(n, mean_fc, fc_dist, fc_sd) {
    if (n < 1L || fc_dist == "fixed" || !is.finite(fc_sd) || fc_sd <= 0) {
        return(rep(mean_fc, n))
    }

    if (fc_dist == "lognormal") {
        return(stats::rlnorm(
            n,
            meanlog = log(mean_fc) - 0.5 * fc_sd * fc_sd,
            sdlog = fc_sd
        ))
    }

    if (fc_dist == "gamma") {
        shape <- 1 / (fc_sd * fc_sd)
        scale <- mean_fc / shape
        return(stats::rgamma(n, shape = shape, scale = scale))
    }

    stop("Unsupported `fc_dist`: ", fc_dist, call. = FALSE)
}

.sim_sample_fc_matrix <- function(fc, group_idx, sample_names) {
    sample_fc <- fc[, group_idx, drop = FALSE]
    colnames(sample_fc) <- sample_names
    sample_fc
}

.sim_ts_pattern_names <- function() {
    c(
        "increase",
        "decrease",
        "early_peak",
        "late_peak",
        "biphasic",
        "oscillatory",
        "delayed_increase",
        "null"
    )
}

.sim_ts_p_patterns <- function(p_patterns) {
    pattern_names <- .sim_ts_pattern_names()

    if (is.null(p_patterns)) {
        return(stats::setNames(rep(1 / length(pattern_names), length(pattern_names)), pattern_names))
    }

    input_names <- names(p_patterns)
    p_patterns <- as.numeric(p_patterns)

    if (is.null(input_names)) {
        if (length(p_patterns) != length(pattern_names)) {
            stop(
                "An unnamed `p_patterns` vector must have length eight.",
                call. = FALSE
            )
        }
        names(p_patterns) <- pattern_names
        probs <- p_patterns
    } else {
        if (anyNA(input_names) || any(input_names == "") || anyDuplicated(input_names)) {
            stop("Named `p_patterns` must have unique, non-empty names.", call. = FALSE)
        }
        unknown <- setdiff(input_names, pattern_names)
        if (length(unknown) > 0L) {
            stop(
                "Unknown temporal pattern name(s): ",
                paste(unknown, collapse = ", "),
                call. = FALSE
            )
        }
        probs <- stats::setNames(rep(0, length(pattern_names)), pattern_names)
        probs[input_names] <- p_patterns
    }

    if (any(!is.finite(probs)) || any(probs < 0) || sum(probs) <= 0) {
        stop("`p_patterns` must contain non-negative finite values with a positive sum.", call. = FALSE)
    }
    probs / sum(probs)
}

.sim_ts_assign_patterns <- function(n_genes, probs) {
    expected <- n_genes * probs
    counts <- floor(expected)
    remainder <- n_genes - sum(counts)

    if (remainder > 0L) {
        frac <- expected - counts
        add_idx <- order(frac, decreasing = TRUE)[seq_len(remainder)]
        counts[add_idx] <- counts[add_idx] + 1L
    }

    patterns <- rep(names(probs), times = counts)
    sample(patterns, length(patterns), replace = FALSE)
}

.sim_ts_pattern_matrix <- function(ts) {
    u <- (ts - min(ts)) / diff(range(ts))

    raw <- rbind(
        increase = 2 * u - 1,
        decrease = 1 - 2 * u,
        early_peak = exp(-0.5 * ((u - 0.25) / 0.12)^2),
        late_peak = exp(-0.5 * ((u - 0.75) / 0.12)^2),
        biphasic = sin(2 * pi * u),
        oscillatory = sin(4 * pi * u),
        delayed_increase = stats::plogis((u - 0.65) / 0.08),
        null = rep(0, length(u))
    )

    templates <- t(apply(raw, 1L, .sim_ts_center_scale))
    colnames(templates) <- .sim_time_group_names(ts)
    templates
}

.sim_ts_center_scale <- function(x) {
    if (all(x == 0)) {
        return(x)
    }
    x_min <- min(x)
    x_max <- max(x)
    half_range <- (x_max - x_min) / 2
    if (!is.finite(half_range) || half_range == 0) {
        return(rep(0, length(x)))
    }
    midpoint <- (x_max + x_min) / 2
    (x - midpoint) / half_range
}

.sim_ts_fc_matrix <- function(patterns,
                              templates,
                              dynamic_fc,
                              gene_names,
                              group_names) {
    pattern_idx <- match(patterns, rownames(templates))
    scores <- templates[pattern_idx, , drop = FALSE]

    # Templates are range-scaled to [-1, 1]. Multiplying by half of
    # log2(dynamic_fc) makes the maximum-to-minimum fold-change ratio equal to
    # dynamic_fc for every non-null template.
    log2_amplitude <- 0.5 * log2(dynamic_fc)
    log2_fc <- sweep(scores, 1L, log2_amplitude, "*")
    fc <- 2^log2_fc
    dimnames(fc) <- list(gene_names, group_names)
    fc
}

.sim_time_group_names <- function(ts) {
    labels <- format(ts, trim = TRUE, scientific = FALSE)
    labels <- gsub("[^[:alnum:]_.-]", "_", labels)
    make.unique(paste0("time_", labels), sep = "_")
}
