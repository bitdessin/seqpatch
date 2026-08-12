#' Class to Store Simulation Parameters
#'
#' An S4 class for storing the parameters used to generate synthetic RNA-seq
#' count data with [sim_gene_counts()].
#'
#' @slot n_genes Number of genes.
#' @slot n_groups Number of groups.
#' @slot n_replicates Number of replicates per group.
#' @slot p_DEG Proportion of genes assigned as DEGs for each group.
#' @slot fc_dist Fold change distribution used for genes assigned as DEGs.
#' @slot fc_sd Spread parameter used when `fc_dist` is not `"fixed"`.
#' @slot fc Group-level fold change matrix. Rows are genes and columns are
#'   groups.
#' @slot params_population Data frame containing the seed-derived population of
#'   mean, variance, dispersion, and distribution values.
#' @slot params Data frame containing the sampled mean, variance, dispersion,
#'   and distribution values used for each simulated gene.
#'
#' @return An object of class \linkS4class{SimParams}.
#' @seealso [sim_gene_counts()]
#' @exportClass SimParams
setClass(
    "SimParams",
    slots = c(
        n_genes = "numeric",
        n_groups = "numeric",
        n_replicates = "numeric",
        p_DEG = "numeric",
        fc_dist = "character",
        fc_sd = "numeric",
        fc = "matrix",
        params_population = "data.frame",
        params = "data.frame"
    )
)

setValidity("SimParams", function(object) {
    errors <- character()

    if (length(object@n_genes) != 1L || object@n_genes < 1L) {
        errors <- c(errors, "`n_genes` must be a positive scalar.")
    }
    if (length(object@n_groups) != 1L || object@n_groups < 1L) {
        errors <- c(errors, "`n_groups` must be a positive scalar.")
    }
    if (length(object@n_replicates) != object@n_groups) {
        errors <- c(errors, "`n_replicates` must have one value per group.")
    }
    if (length(object@p_DEG) != object@n_groups) {
        errors <- c(errors, "`p_DEG` must have one value per group.")
    }
    if (length(object@fc_dist) != 1L) {
        errors <- c(errors, "`fc_dist` must be a scalar.")
    }
    if (length(object@fc_sd) != object@n_groups) {
        errors <- c(errors, "`fc_sd` must have one value per group.")
    }
    if (nrow(object@fc) != object@n_genes) {
        errors <- c(errors, "`fc` must have one row per gene.")
    }
    if (ncol(object@fc) != object@n_groups) {
        errors <- c(errors, "`fc` must have one column per group.")
    }
    if (nrow(object@params) != object@n_genes) {
        errors <- c(errors, "`params` must have one row per gene.")
    }

    if (length(errors)) {
        errors
    } else {
        TRUE
    }
})
