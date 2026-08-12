#' Class to Store RNA-seq Read Count Data
#'
#' An S4 class for storing RNA-seq read count data, gene names, experimental
#' design, and associated metadata.
#'
#' @slot data A read count matrix. Rows are genes and columns are samples.
#' @slot gene_names A character vector containing gene names.
#' @slot exp_design A data frame describing the experimental design.
#'      It must include a column named `group` indicating sample groupings.
#' @slot meta A list containing additional metadata related to the experiment
#'      or processing steps.
#'
#' @return An object of class \linkS4class{SeqCountData}.
#' @seealso [newSeqCountData()]
#' @exportClass SeqCountData
setClass(
    "SeqCountData",
    slots = c(
        data = "matrix",
        gene_names = "character",
        exp_design = "data.frame",
        meta = "list"
    )
)

setValidity("SeqCountData", function(object) {
    errors <- character()
    data <- object@data

    if (!is.numeric(data)) {
        errors <- c(errors, "`data` must be a numeric matrix.")
    }
    if (anyNA(data)) {
        errors <- c(errors, "`data` must not contain missing values.")
    }
    if (any(!is.finite(data))) {
        errors <- c(errors, "`data` must contain only finite values.")
    }
    if (any(data < 0)) {
        errors <- c(errors, "`data` must not contain negative values.")
    }
    if (nrow(object@exp_design) != ncol(data)) {
        errors <- c(errors, "`exp_design` must have one row per sample.")
    }
    if (!"group" %in% names(object@exp_design)) {
        errors <- c(errors, "`exp_design` must contain a `group` column.")
    } else if (anyNA(object@exp_design$group)) {
        errors <- c(errors, "`exp_design$group` must not contain missing values.")
    }
    if (!is.null(object@gene_names)) {
        if (length(object@gene_names) != nrow(data)) {
            errors <- c(errors, "`gene_names` must have one value per row of `data`.")
        }
        if (anyNA(object@gene_names)) {
            errors <- c(errors, "`gene_names` must not contain missing values.")
        }
    }

    if (length(errors)) {
        errors
    } else {
        TRUE
    }
})

#' Create a SeqCountData Object
#'
#' @param data Matrix object, edgeR `DGEList`, DESeq2 `DESeqDataSet`, or
#'   `SummarizedExperiment`. For matrix input, rows are genes and columns are
#'   samples.
#' @param exp_design Data frame describing sample groupings. For `DGEList`,
#'   `DESeqDataSet`, and `SummarizedExperiment` input, this defaults to the
#'   sample metadata stored in the input object. It must include a column named
#'   `group`.
#' @param gene_names Optional vector of gene names. Defaults to row names of
#'   `data` when present.
#' @param meta Optional list containing additional metadata.
#' @param norm_factors Optional numeric normalization factors. If provided, they
#'   are stored as `meta$nf`. For `DGEList` input, this defaults to
#'   `data$samples$norm.factors` when available. For `DESeqDataSet` input, this
#'   is derived from `DESeq2::sizeFactors(data)` when available.
#'
#' @return A `SeqCountData` object.
#'
#' @examples
#' counts <- matrix(1:12, nrow = 3)
#' colnames(counts) <- paste0("sample_", seq_len(ncol(counts)))
#' rownames(counts) <- paste0("gene_", seq_len(nrow(counts)))
#' exp_design <- data.frame(group = rep(c("A", "B"), each = 2))
#' x <- newSeqCountData(counts, exp_design = exp_design)
#' x
#'
#' @export
newSeqCountData <- function(data,
                gene_names = NULL,
                exp_design = NULL,
                meta = list(),
                norm_factors = NULL) {
    if (inherits(data, "DGEList")) {
        return(.newSeqCountData_from_DGEList(data, gene_names, exp_design, meta, norm_factors))
    }
    if (methods::is(data, "DESeqDataSet")) {
        return(.newSeqCountData_from_DESeqDataSet(data, gene_names, exp_design, meta, norm_factors))
    }
    if (methods::is(data, "SummarizedExperiment")) {
        return(.newSeqCountData_from_SummarizedExperiment(data, gene_names, exp_design, meta, norm_factors))
    }

    .newSeqCountData_from_matrix(data, gene_names, exp_design, meta, norm_factors)
}

.newSeqCountData_from_matrix <- function(data,
                                  gene_names = NULL,
                                  exp_design = NULL,
                                  meta = list(),
                                  norm_factors = NULL) {
    data <- as.matrix(data)
    if (is.null(meta)) meta <- list()
    if (is.vector(exp_design)) exp_design <- data.frame(group  = exp_design)
    exp_design <- as.data.frame(exp_design)
    if (is.null(gene_names)) {
        if (!is.null(rownames(data))) {
            gene_names <- rownames(data)
        }
    }
    if (!is.null(norm_factors)) {
        meta$nf <- .format_norm_factors(norm_factors, colnames(data))
    }

    new("SeqCountData", data = data, gene_names = gene_names, exp_design = exp_design, meta = meta)
}

.newSeqCountData_from_DGEList <- function(data,
                                   gene_names = NULL,
                                   exp_design = NULL,
                                   meta = list(),
                                   norm_factors = NULL) {
    counts <- as.matrix(data$counts)
    if (is.null(exp_design)) {
        exp_design <- as.data.frame(data$samples)
    }
    if (is.null(norm_factors) && "norm.factors" %in% names(data$samples)) {
        norm_factors <- data$samples$norm.factors
    }

    .newSeqCountData_from_matrix(counts, gene_names, exp_design, meta, norm_factors)
}

.newSeqCountData_from_DESeqDataSet <- function(data,
                                        gene_names = NULL,
                                        exp_design = NULL,
                                        meta = list(),
                                        norm_factors = NULL) {
    .require_namespace("DESeq2")
    .require_namespace("SummarizedExperiment")

    counts <- as.matrix(DESeq2::counts(data))
    if (is.null(exp_design)) {
        exp_design <- as.data.frame(SummarizedExperiment::colData(data))
    }
    if (is.null(norm_factors)) {
        size_factors <- DESeq2::sizeFactors(data)
        if (!is.null(size_factors)) {
            norm_factors <- size_factors / colSums(counts)
        }
    }

    .newSeqCountData_from_matrix(counts, gene_names, exp_design, meta, norm_factors)
}

.newSeqCountData_from_SummarizedExperiment <- function(data,
                                                gene_names = NULL,
                                                exp_design = NULL,
                                                meta = list(),
                                                norm_factors = NULL) {
    .require_namespace("SummarizedExperiment")

    assay_names <- SummarizedExperiment::assayNames(data)
    assay_name <- if ("counts" %in% assay_names) "counts" else 1L
    counts <- as.matrix(SummarizedExperiment::assay(data, assay_name))
    if (is.null(exp_design)) {
        exp_design <- as.data.frame(SummarizedExperiment::colData(data))
    }

    .newSeqCountData_from_matrix(counts, gene_names, exp_design, meta, norm_factors)
}

.format_norm_factors <- function(norm_factors, sample_names) {
    nf_names <- names(norm_factors)
    norm_factors <- as.numeric(norm_factors)
    if (!is.null(nf_names)) {
        names(norm_factors) <- nf_names
    } else {
        names(norm_factors) <- sample_names
    }
    norm_factors
}

#' Show a SeqCountData Object
#'
#' @param object A `SeqCountData` object.
#'
#' @examples
#' counts <- matrix(1:12, nrow = 3)
#' colnames(counts) <- paste0("sample_", seq_len(ncol(counts)))
#' rownames(counts) <- paste0("gene_", seq_len(nrow(counts)))
#' exp_design <- data.frame(group = rep(c("A", "B"), each = 2))
#' x <- newSeqCountData(counts, exp_design = exp_design)
#' show(x)
#'
#' @export
setMethod("show", "SeqCountData", function(object) {
    nf <- object@meta$nf
    cat("SeqCountData object\n")
    cat("  genes:   ", nrow(object@data), "\n", sep = "")
    cat("  samples: ", ncol(object@data), "\n", sep = "")
    cat("  groups:  ",
        paste(levels(factor(object@exp_design$group)), collapse = ", "),
        "\n", sep = "")
    if (is.null(nf)) {
        cat("  factors: not set\n", sep = "")
    } else {
        cat("  factors: ",
            paste(signif(nf, 4), collapse = ", "),
            "\n", sep = "")
    }
})

#' Convert to an edgeR DGEList
#'
#' @param x A `SeqCountData` object.
#' @return An edgeR `DGEList`.
#'
#' @examples
#' if (requireNamespace("edgeR", quietly = TRUE)) {
#'     counts <- matrix(1:12, nrow = 3)
#'     colnames(counts) <- paste0("sample_", seq_len(ncol(counts)))
#'     rownames(counts) <- paste0("gene_", seq_len(nrow(counts)))
#'     exp_design <- data.frame(group = rep(c("A", "B"), each = 2))
#'     x <- newSeqCountData(
#'         counts,
#'         exp_design = exp_design,
#'         norm_factors = rep(1, ncol(counts))
#'     )
#'     d <- as_DGEList(x)
#'     d$samples$norm.factors
#' }
#'
#' @export
as_DGEList <- function(x) {
    .require_namespace("edgeR")
    d <- edgeR::DGEList(counts = x@data, group = x@exp_design$group)
    d$samples$norm.factors <- x@meta$nf
    d
}

#' Convert to a SummarizedExperiment
#'
#' @param x A `SeqCountData` object.
#' @return A `SummarizedExperiment` object with counts stored in the `counts`
#'   assay, `x@exp_design` stored as `colData`, gene names stored as `rowData`,
#'   and `x@meta` stored as object metadata.
#'
#' @examples
#' if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
#'     counts <- matrix(1:12, nrow = 3)
#'     colnames(counts) <- paste0("sample_", seq_len(ncol(counts)))
#'     rownames(counts) <- paste0("gene_", seq_len(nrow(counts)))
#'     exp_design <- data.frame(group = rep(c("A", "B"), each = 2))
#'     x <- newSeqCountData(counts, exp_design = exp_design)
#'     se <- as_SummarizedExperiment(x)
#'     SummarizedExperiment::assayNames(se)
#' }
#'
#' @export
as_SummarizedExperiment <- function(x) {
    .require_namespace("SummarizedExperiment")

    row_data <- data.frame(
        gene_name = as.character(.get_gene_names(x)),
        stringsAsFactors = FALSE
    )
    col_data <- x@exp_design
    rownames(col_data) <- colnames(x@data)

    SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = x@data),
        rowData = row_data,
        colData = col_data,
        metadata = x@meta
    )
}

#' Convert to a DESeq2 DESeqDataSet
#'
#' @param x A `SeqCountData` object.
#' @param design A DESeq2 design formula. Defaults to `~ group`.
#' @return A DESeq2 `DESeqDataSet`.
#'
#' @examples
#' if (requireNamespace("DESeq2", quietly = TRUE)) {
#'     counts <- matrix(1:12, nrow = 3)
#'     colnames(counts) <- paste0("sample_", seq_len(ncol(counts)))
#'     rownames(counts) <- paste0("gene_", seq_len(nrow(counts)))
#'     exp_design <- data.frame(group = rep(c("A", "B"), each = 2))
#'     x <- newSeqCountData(
#'         counts,
#'         exp_design = exp_design,
#'         norm_factors = rep(1, ncol(counts))
#'     )
#'     dds <- as_DESeqDataSet(x)
#'     DESeq2::design(dds)
#' }
#'
#' @export
as_DESeqDataSet <- function(x, design = ~ group) {
    .require_namespace("DESeq2")
    col_data <- x@exp_design
    rownames(col_data) <- colnames(x@data)
    d <- DESeq2::DESeqDataSetFromMatrix(
        countData = x@data,
        colData = col_data,
        design = design)
    DESeq2::sizeFactors(d) <- .calc_sizefactors(x@data, x@meta$nf)
    d
}
