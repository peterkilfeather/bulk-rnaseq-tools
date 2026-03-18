#' Normalise a raw counts matrix using edgeR/limma workflows
#'
#' @param counts Raw counts matrix (genes × samples), or a file path to a TSV.
#' @param group  Optional factor/character vector of group labels per sample.
#'               When NULL, an intercept-only design is used.
#' @param filter Logical; apply edgeR::filterByExpr before normalization.
#'
#' @return A named list with elements:
#'   \item{logcpm}{logCPM matrix (edgeR::cpm, log = TRUE)}
#'   \item{voom}{limma-voom normalized expression matrix}
#'   \item{voom_sample_weights}{limma-voom with sample quality weights}
#'   \item{kept_genes}{Logical vector indicating genes that passed filtering}
#'   \item{dge}{DGEList after filtering and TMM normalization}
normalise_counts <- function(counts, group = NULL, filter = TRUE) {

    # -- Input handling --------------------------------------------------------
    if (is.character(counts) && length(counts) == 1L) {
        counts <- read.delim(counts, row.names = 1, check.names = FALSE)
    }
    counts <- as.matrix(counts)

    # -- DGEList ---------------------------------------------------------------
    dge <- edgeR::DGEList(counts = counts)

    # -- Design matrix ---------------------------------------------------------
    if (!is.null(group)) {
        group  <- factor(group)
        design <- model.matrix(~ group)
    } else {
        design <- model.matrix(~ 1,
            data = data.frame(row.names = colnames(counts)))
    }

    # -- Filtering -------------------------------------------------------------
    if (filter) {
        kept_genes <- edgeR::filterByExpr(dge, design = design)
        dge <- dge[kept_genes, , keep.lib.sizes = FALSE]
    } else {
        kept_genes <- rep(TRUE, nrow(dge))
        names(kept_genes) <- rownames(dge)
    }

    # -- TMM normalization -----------------------------------------------------
    dge <- edgeR::calcNormFactors(dge)

    # -- logCPM ----------------------------------------------------------------
    logcpm <- edgeR::cpm(dge, log = TRUE)

    # -- Voom ------------------------------------------------------------------
    voom_result <- limma::voom(dge, design = design)

    # -- Voom with sample weights ----------------------------------------------
    voom_sw_result <- limma::voomWithQualityWeights(dge, design = design)

    # -- Return ----------------------------------------------------------------
    list(
        logcpm              = logcpm,
        voom                = voom_result$E,
        voom_sample_weights = voom_sw_result$E,
        kept_genes          = kept_genes,
        dge                 = dge
    )
}
