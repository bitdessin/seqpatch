#' Estimate Normalization Factors
#'
#' Performs a compact DEGES workflow:
#'
#' 1. Estimate normalization factors with `norm_func`.
#' 2. Screen candidate DEGs with `test_func`.
#' 3. Remove genes passing `p_cutoff` or `q_cutoff`. If fewer than
#'    `min_rm_genes` of genes would be removed, remove the top
#'    `min_rm_genes` proportion ranked by p-value.
#' 4. Estimate normalization factors again from the remaining genes.
#'
#' Steps 2 through 4 are repeated `iter` times.
#'
#' @param x A `SeqCountData` object.
#' @param norm_func Function used to estimate normalization factors. It should
#'   accept a count matrix `x`, sample metadata `exp_design`, and optional
#'   method-specific arguments.
#' @param test_func Function used for DEG testing. It should
#'   accept a count matrix `x`, sample metadata `exp_design`, normalization
#'   factors `nf`, and optional method-specific arguments.
#' @param iter Number of DEG elimination iterations after the initial
#'   normalization.
#' @param p_cutoff p-value cutoff for removing potential DEGs. Use `NULL`
#'   to ignore p-value filtering.
#' @param q_cutoff Adjusted p-value cutoff for removing potential DEGs. Use
#'   `NULL` to ignore adjusted p-value filtering.
#' @param min_rm_genes Minimum proportion of genes to remove at each iteration.
#' @param ... Additional arguments passed to `norm_func` and `test_func` when
#'   their signatures accept them.
#'
#' @return A `SeqCountData` object with updated normalization factors. DEGES
#'   settings and potential DEG indicators are stored in `x@meta$DEGES`.
#'
#' @examples
#' if (requireNamespace("edgeR", quietly = TRUE)) {
#'     set.seed(1)
#'     counts <- matrix(stats::rnbinom(600, mu = 50, size = 5), nrow = 100)
#'     colnames(counts) <- paste0("sample_", seq_len(ncol(counts)))
#'     rownames(counts) <- paste0("gene_", seq_len(nrow(counts)))
#'     exp_design <- data.frame(group = rep(c("A", "B"), each = 3))
#'     x <- newSeqCountData(counts, exp_design = exp_design)
#'     x <- calc_nf(x, iter = 1)
#'     x@meta$nf
#' }
#'
#' @export
calc_nf <- function(
    x,
    norm_func = norm_tmm,
    test_func = test_edger,
    iter = 3,
    p_cutoff = 0.05,
    q_cutoff = 0.10,
    min_rm_genes = 0.05,
    ...) {
    if (!is.function(norm_func)) {
        stop("`norm_func` must be a function.", call. = FALSE)
    }
    if (!is.function(test_func)) {
        stop("`test_func` must be a function.", call. = FALSE)
    }

    nf <- .call_norm_func(norm_func, x@data, x@exp_design, ...)
    x@meta$nf <- nf / mean(nf)
    x@meta$DEGES <- .init_DEGES_meta(
        x,
        iter,
        p_cutoff,
        q_cutoff,
        min_rm_genes
    )

    if (iter == 0L) {
        return(x)
    }

    for (i in seq_len(iter)) {
        # reset data
        dmat <- x@data
        exp_design <- x@exp_design
        nf <- x@meta$nf

        test_stats <- .call_test_func(test_func, dmat, exp_design, nf, ...)
        rm_targets <- .select_rm_targets(test_stats, p_cutoff, q_cutoff, min_rm_genes)
        x@meta$DEGES$potential_DEGs[, i] <- rm_targets
        if (!any(!rm_targets)) {
            warning(
                "All genes were selected for removal; returning factors from ",
                "the previous iteration.",
                call. = FALSE
            )
            break
        }

        dmat <- dmat[!rm_targets, , drop = FALSE]
        nf <- .call_norm_func(norm_func, dmat, exp_design, ...)
        nf <- nf * colSums(dmat) / colSums(x@data)
        x@meta$nf <- nf / mean(nf)
    }

    x
}

.init_DEGES_meta <- function(x, iter, p_cutoff, q_cutoff, min_rm_genes) {
    potential_DEGs <- matrix(FALSE, nrow = nrow(x@data), ncol = iter)
    colnames(potential_DEGs) <- if (iter > 0L) {
        paste0("iter_", seq_len(iter))
    } else {
        character()
    }

    list(
        params = list(
            n_iters = iter,
            p_cutoff = p_cutoff,
            q_cutoff = q_cutoff,
            min_rm_genes = min_rm_genes
        ),
        potential_DEGs = potential_DEGs
    )
}

.call_norm_func <- function(norm_func, x, exp_design, ...) {
    .call_with_supported_args(
        norm_func,
        c(list(x = x, exp_design = exp_design), list(...))
    )
}

.call_test_func <- function(test_func, x, exp_design, nf, ...) {
    .call_with_supported_args(
        test_func,
        c(list(x = x, exp_design = exp_design, nf = nf), list(...))
    )
}

.call_with_supported_args <- function(func, args) {
    formal_names <- names(formals(func))
    if (is.null(formal_names) || "..." %in% formal_names) {
        return(do.call(func, args))
    }
    do.call(func, args[intersect(names(args), formal_names)])
}


.select_rm_targets <- function(stat_df, p_cutoff, q_cutoff, min_rm_genes) {
    n_genes <- nrow(stat_df)
    rm_targets <- rep(FALSE, n_genes)

    if (!is.null(p_cutoff)) {
        rm_targets <- rm_targets | (stat_df$p_value <= p_cutoff)
    }

    if (!is.null(q_cutoff)) {
        rm_targets <- rm_targets | (stat_df$q_value <= q_cutoff)
    }
    rm_targets[is.na(rm_targets)] <- FALSE

    min_count <- ceiling(n_genes * min_rm_genes)
    if (sum(rm_targets) < min_rm_genes && min_count > 0L) {
        p_rank <- stat_df$p_value
        p_rank[is.na(p_rank)] <- Inf
        p_rank <- order(p_rank, decreasing = FALSE)
        rm_targets <- rep(FALSE, n_genes)
        rm_targets[p_rank[seq_len(min_count)]] <- TRUE
    }

    rm_targets
}
