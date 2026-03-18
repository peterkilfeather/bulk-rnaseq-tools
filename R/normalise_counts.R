#' Normalise a raw counts matrix using edgeR/limma workflows
#'
#' @param counts Raw counts matrix (genes × samples), or a file path to a TSV.
#' @param group  Optional factor/character vector of group labels per sample.
#'               When NULL, an intercept-only design is used.
#' @param block  Optional factor/character vector identifying a blocking variable
#'               (e.g., subject ID for repeated measures). When provided,
#'               \code{limma::duplicateCorrelation} is used to estimate
#'               intra-block correlation via a two-iteration procedure
#'               (voom → duplicateCorrelation → voom → duplicateCorrelation),
#'               following limma best practice.
#' @param filter Logical; apply edgeR::filterByExpr before normalization.
#'
#' @return A named list with elements:
#'   \item{logcpm}{logCPM matrix (edgeR::cpm, log = TRUE)}
#'   \item{voom}{limma-voom normalized expression matrix}
#'   \item{voom_sample_weights}{limma-voom with sample quality weights}
#'   \item{kept_genes}{Logical vector indicating genes that passed filtering}
#'   \item{dge}{DGEList after filtering and TMM normalization}
#'   \item{block_correlation}{Consensus intra-block correlation from the voom
#'         path, or NULL when \code{block} is NULL. Pass to
#'         \code{lmFit(correlation = ...)} for downstream analysis.}
#'   \item{block_correlation_sw}{Consensus intra-block correlation from the
#'         voomWithQualityWeights path, or NULL when \code{block} is NULL.}
#'   \item{sample_weights}{Sample quality weights from voomWithQualityWeights.}
normalise_counts <- function(counts, group = NULL, block = NULL, filter = TRUE) {

    # -- Input handling --------------------------------------------------------
    if (is.character(counts) && length(counts) == 1L) {
        message("Reading counts from file...")
        counts <- read.delim(counts, row.names = 1, check.names = FALSE)
    }
    counts <- as.matrix(counts)

    # -- DGEList ---------------------------------------------------------------
    message("Creating DGEList...")
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
        message("Filtering lowly expressed genes...")
        kept_genes <- edgeR::filterByExpr(dge, design = design)
        dge <- dge[kept_genes, , keep.lib.sizes = FALSE]
        message("Kept ", sum(kept_genes), " of ", length(kept_genes), " genes")
    } else {
        kept_genes <- rep(TRUE, nrow(dge))
        names(kept_genes) <- rownames(dge)
    }

    # -- TMM normalization -----------------------------------------------------
    message("Calculating TMM normalization factors...")
    dge <- edgeR::calcNormFactors(dge)

    # -- logCPM ----------------------------------------------------------------
    message("Computing logCPM...")
    logcpm <- edgeR::cpm(dge, log = TRUE)

    # -- Voom ------------------------------------------------------------------
    if (!is.null(block)) {
        block <- as.factor(block)

        # Voom with duplicateCorrelation (two iterations)
        message("Running voom...")
        v1 <- limma::voom(dge, design = design)
        message("Estimating duplicate correlation (pass 1)...")
        corfit <- limma::duplicateCorrelation(v1, design = design, block = block)
        message("Running voom with block correlation...")
        v2 <- limma::voom(dge, design = design,
                          block = block,
                          correlation = corfit$consensus.correlation)
        message("Estimating duplicate correlation (pass 2)...")
        corfit <- limma::duplicateCorrelation(v2, design = design, block = block)
        voom_result <- v2

        # Voom with sample weights + duplicateCorrelation (two iterations)
        message("Running voom with sample quality weights...")
        vsw1 <- limma::voomWithQualityWeights(dge, design = design)
        message("Estimating duplicate correlation for sample weights (pass 1)...")
        corfit_sw <- limma::duplicateCorrelation(vsw1, design = design,
                                                 block = block)
        message("Running voom with sample quality weights and block correlation...")
        vsw2 <- limma::voomWithQualityWeights(dge, design = design,
                          block = block,
                          correlation = corfit_sw$consensus.correlation)
        message("Estimating duplicate correlation for sample weights (pass 2)...")
        corfit_sw <- limma::duplicateCorrelation(vsw2, design = design,
                                                 block = block)
        voom_sw_result <- vsw2
    } else {
        message("Running voom...")
        voom_result <- limma::voom(dge, design = design)

        message("Running voom with sample quality weights...")
        voom_sw_result <- limma::voomWithQualityWeights(dge, design = design)

        corfit <- NULL
        corfit_sw <- NULL
    }

    # -- Return ----------------------------------------------------------------
    message("Done.")
    list(
        logcpm               = logcpm,
        voom                 = voom_result$E,
        voom_sample_weights  = voom_sw_result$E,
        kept_genes           = kept_genes,
        dge                  = dge,
        block_correlation    = if (!is.null(corfit)) corfit$consensus.correlation,
        block_correlation_sw = if (!is.null(corfit_sw)) corfit_sw$consensus.correlation,
        sample_weights       = voom_sw_result$targets$sample.weights
    )
}
