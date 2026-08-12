#' Plot an MA Plot
#'
#' Creates an MA plot from normalized counts for two groups in `x@exp_design$group`.
#' Group averages are computed before log transformation for the A values.
#'
#' @param x A `SeqCountData` object.
#' @param comparison Optional vector of two group names to compare. If
#'   `NULL`, the first two groups in `x@exp_design$group` are used.
#' @param esp Value added before log transformation.
#' @param col Optional point colors. Can be a single color or one color per gene.
#'   When provided, points are drawn from the most frequent color category to
#'   the least frequent category.
#' @return A `ggplot` object.
#' @importFrom ggplot2 aes geom_point ggplot labs scale_color_identity
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'     set.seed(1)
#'     counts <- matrix(stats::rnbinom(600, mu = 50, size = 5), nrow = 100)
#'     colnames(counts) <- paste0("sample_", seq_len(ncol(counts)))
#'     rownames(counts) <- paste0("gene_", seq_len(nrow(counts)))
#'     exp_design <- data.frame(group = rep(c("A", "B"), each = 3))
#'     x <- newSeqCountData(
#'         counts,
#'         exp_design = exp_design,
#'         norm_factors = rep(1, ncol(counts))
#'     )
#'     plot_ma(x)
#' }
#'
#' @export
plot_ma <- function(x, comparison = NULL, esp = 1, col = NULL) {
    group_vec <- factor(x@exp_design$group)
    groups <- levels(group_vec)
    if (length(groups) < 2L) {
        stop("At least two groups are required for an MA plot.", call. = FALSE)
    }
    if (is.null(comparison)) {
        comparison <- groups[seq_len(2L)]
    } else {
        if (length(comparison) != 2L) {
            stop("`comparison` must be NULL or a length-2 vector.", call. = FALSE)
        }
        comparison <- as.character(comparison)
        missing_groups <- setdiff(comparison, groups)
        if (length(missing_groups) > 0L) {
            stop(
                "`comparison` contains groups not present in `x@exp_design$group`: ",
                paste(missing_groups, collapse = ", "),
                call. = FALSE
            )
        }
    }
    data <- norm_counts(x)

    y1 <- rowMeans(data[, group_vec == comparison[1], drop = FALSE])
    y2 <- rowMeans(data[, group_vec == comparison[2], drop = FALSE])
    log_y1 <- log2(y1 + esp)
    log_y2 <- log2(y2 + esp)

    data_df <- data.frame(
        A = unname(log2(((y1 + y2) / 2) + esp)),
        M = unname(log_y2 - log_y1),
        id = as.character(.get_gene_names(x)),
        check.names = FALSE
    )

    if (!is.null(col)) {
        data_df$col <- col
        col_counts <- sort(table(data_df$col), decreasing = TRUE)
        data_df <- data_df[
            order(match(data_df$col, names(col_counts))),
            ,
            drop = FALSE
        ]
    }

    p <- ggplot(data_df, aes(x = .data[["A"]], y = .data[["M"]], text = .data[["id"]]))
    if (is.null(col)) {
        p <- p + geom_point(alpha = 0.4, size = 1.2)
    } else {
        p <- p +
            geom_point(aes(color = .data[["col"]]), alpha = 0.4, size = 1.2) +
            scale_color_identity()
    }

    p + labs(x = paste0('log2((', comparison[1], ' + ', comparison[2], ') / 2)'),
             y = paste0('log2(', comparison[2], ') - log2(', comparison[1], ')'),
             title = paste(comparison[2], "vs", comparison[1]))
}

#' Plot a PCA Projection
#'
#' Creates a PCA plot from log-normalized counts.
#'
#' @param x A `SeqCountData` object.
#' @param pc_x Principal component to show on the horizontal axis.
#' @param pc_y Principal component to show on the vertical axis.
#' @param esp Value added before log transformation.
#' @return A `ggplot` object.
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'     set.seed(1)
#'     counts <- matrix(stats::rnbinom(600, mu = 50, size = 5), nrow = 100)
#'     colnames(counts) <- paste0("sample_", seq_len(ncol(counts)))
#'     rownames(counts) <- paste0("gene_", seq_len(nrow(counts)))
#'     exp_design <- data.frame(group = rep(c("A", "B"), each = 3))
#'     x <- newSeqCountData(
#'         counts,
#'         exp_design = exp_design,
#'         norm_factors = rep(1, ncol(counts))
#'     )
#'     plot_pca(x)
#' }
#'
#' @export
plot_pca <- function(x, pc_x = 1, pc_y = 2, esp = 1) {
    data <- log2(norm_counts(x) + esp)
    max_pc <- min(nrow(data), ncol(data))
    if (max_pc < 2L) {
        stop("At least two genes and two samples are required for a PCA plot.", call. = FALSE)
    }
    if (pc_x < 1L || pc_y < 1L || pc_x > max_pc || pc_y > max_pc) {
        stop("`pc_x` and `pc_y` must refer to available principal components.", call. = FALSE)
    }

    pca <- stats::prcomp(t(data), center = TRUE, scale. = FALSE)
    variance <- 100 * (pca$sdev^2 / sum(pca$sdev^2))

    data_df <- data.frame(
        sample_id = colnames(data),
        x = pca$x[, pc_x],
        y = pca$x[, pc_y],
        group = x@exp_design$group,
        check.names = FALSE
    )

    ggplot(data_df,
           aes(x = .data[["x"]], y = .data[["y"]], color = .data[["group"]])) +
        geom_point(size = 2.4) +
        labs(x = paste0("PC", pc_x, " (", round(variance[pc_x], 1), "%)"),
             y = paste0("PC", pc_y, " (", round(variance[pc_y], 1), "%)"),
             color = "group")
}


#' Plot Gene Expression Profiles
#'
#' Plots normalized expression profiles for genes grouped by a cluster or other
#' label. Individual gene profiles are shown as light lines, and the mean
#' profile for each gene group is overlaid as a darker line.
#'
#' Normalized counts are obtained with [norm_counts()]. When `by` is specified,
#' normalized counts are averaged across samples with the same value in the
#' corresponding column of `x@exp_design`. The values are then optionally log
#' transformed and scaled to z-scores within each gene.
#'
#' @param x A `SeqCountData` object with normalization factors stored in
#'   `x@meta$nf`.
#' @param cl Integer or character vector assigning each gene to a group or
#'   cluster. Its length must equal the number of genes in `x`.
#' @param by Column name in `x@exp_design` used to average biological
#'   replicates before plotting. Defaults to `"group"`. If `NULL`, every sample
#'   is plotted separately.
#' @param zscale Logical value indicating whether expression values should be
#'   converted to gene-wise z-scores.
#' @param log_transform Logical value indicating whether normalized counts
#'   should be transformed with `log2(value + pseudocount)`.
#' @param pseudocount Non-negative value added before log transformation.
#' @param ncol Number of columns used to arrange the gene-group facets.
#' @param line_color Color used for individual gene profiles.
#' @param line_alpha Opacity used for individual gene profiles.
#' @param line_width Width used for individual gene profiles.
#' @param mean_color Color used for gene-group mean profiles.
#' @param mean_width Width used for gene-group mean profiles.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' set.seed(1)
#' x <- sim_gene_ts_counts(2000)
#' cl <- x@meta$ts_truth$pattern
#'
#' x <- calc_nf(x)
#' plot_gene_profiles(x, cl)
#'
#' @importFrom methods is
#' @importFrom stats sd aggregate
#' @importFrom ggplot2 aes element_text facet_wrap geom_line ggplot labs theme vars
#' @export
plot_gene_profiles <- function(x,
                               cl,
                               by = "group",
                               zscale = TRUE,
                               log_transform = TRUE,
                               pseudocount = 1,
                               ncol = NULL,
                               line_color = "#7A7A7A",
                               line_alpha = 0.15,
                               line_width = 0.35,
                               mean_color = "#202020",
                               mean_width = 1.1) {

    data <- norm_counts(x)
    sample_labels <- colnames(data)
    if (is.null(sample_labels)) {
        sample_labels <- paste0("sample_", seq_len(ncol(data)))
    }

    if (!is.null(by)) {
        if (length(by) != 1L || !is.character(by) || !by %in% names(x@exp_design)) {
            stop("`by` must name one column in `x@exp_design`.", call. = FALSE)
        }
        by_values <- x@exp_design[[by]]
        if (anyNA(by_values)) {
            stop("The experimental design column selected by `by` contains missing values.", call. = FALSE)
        }
        sample_labels <- unique(as.character(by_values))
        data <- vapply(sample_labels, function(label) {
                rowMeans(data[, as.character(by_values) == label, drop = FALSE])
            }, numeric(nrow(data)))
        data <- matrix(data, nrow = nrow(x@data),
            dimnames = list(rownames(x@data), sample_labels))
    }

    if (log_transform) {
        data <- log2(data + pseudocount)
    }

    if (zscale) {
        gene_means <- rowMeans(data)
        gene_sds <- apply(data, 1L, stats::sd)
        gene_sds[!is.finite(gene_sds) | gene_sds == 0] <- 1
        data <- sweep(sweep(data, 1L, gene_means, "-"), 1L, gene_sds, "/")
    }

    gene_names <- as.character(.get_gene_names(x))
    cluster_levels <- unique(as.character(cl))
    cluster_sizes <- table(factor(as.character(cl), levels = cluster_levels))
    cluster_labels <- paste0(cluster_levels, " (n = ", as.integer(cluster_sizes), ")")
    names(cluster_labels) <- cluster_levels

    profile_df <- data.frame(
        gene = rep(gene_names, times = ncol(data)),
        sample = factor(rep(sample_labels, each = nrow(data)), levels = sample_labels),
        expression = as.vector(data),
        cluster = factor(rep(cluster_labels[as.character(cl)], times = ncol(data)),
                         levels = unname(cluster_labels)),
        check.names = FALSE)

    mean_df <- stats::aggregate(
        profile_df$expression,
        by = list(cluster = profile_df$cluster, sample = profile_df$sample),
        FUN = mean
    )
    names(mean_df)[3L] <- "expression"
    mean_df$cluster <- factor(mean_df$cluster, levels = levels(profile_df$cluster))
    mean_df$sample <- factor(mean_df$sample, levels = sample_labels)

    y_label <- if (zscale) {
        "Gene-wise z-score"
    } else if (log_transform) {
        "log2 normalized expression"
    } else {
        "Normalized expression"
    }

    ggplot(profile_df, aes(x = .data[["sample"]], y = .data[["expression"]], group = .data[["gene"]])) +
        ggplot2::geom_line(color = line_color, alpha = line_alpha, linewidth = line_width) +
        ggplot2::geom_line(data = mean_df,
            aes(x = .data[["sample"]], y = .data[["expression"]], group = .data[["cluster"]]),
            inherit.aes = FALSE, color = mean_color, linewidth = mean_width) +
        ggplot2::facet_wrap(ggplot2::vars(.data[["cluster"]]), ncol = ncol, scales = 'fixed') +
        labs(x = NULL, y = y_label) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}
