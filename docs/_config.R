rm(list = ls(all = TRUE))

library(seqpatch)
options(mc.cores = parallel::detectCores(logical = TRUE))
options(width = 10000)
options(ggplot2.discrete.colour = c("#999999", "#E69F00", "#56B4E9", "#009E73",
                                    "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))
docs_dir <- normalizePath(getwd(), mustWork = TRUE)
if (!file.exists(file.path(docs_dir, "_bookdown.yml"))) {
    stop("Render the book from the `docs` directory.", call. = FALSE)
}
knitr::opts_knit$set(root.dir = docs_dir)

knitr::opts_chunk$set(tidy = FALSE,
                      cache = FALSE,
                      error = FALSE,
                      warning = FALSE,
                      message = FALSE,
                      dev = "png")
set.seed(1)


pkg <- function(pkg_name, repos = NULL) {
    pkg_html <- ''
    if (is.null(repos)) {
        if (pkg_name == 'seqpatch') {
            pkg_html <- paste0('<span class="pkg-name">',
                               '<a href="https://bitdessin.github.io/seqpatch" target="_blank">',
                               pkg_name, '</a></span>')
        } else {
            pkg_html <- paste0('<span class="pkg-name">', pkg_name, '</span>')
        }
    } else if (tolower(repos) == 'cran') {
        pkg_html <- paste0('<span class="pkg-name">',
                           '<a href="https://CRAN.R-project.org/package=',
                           pkg_name, '" target="_blank">', pkg_name, '</a></span>')
    } else if (tolower(repos) == 'bioc') {
        pkg_html <- paste0('<span class="pkg-name">',
                           '<a href="https://bioconductor.org/packages/',
                           pkg_name, '" target="_blank">', pkg_name, '</a></span>')
    }
    pkg_html
}
