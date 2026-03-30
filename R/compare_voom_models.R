#' Compare voom model configurations for differential expression testing
#'
#' Fits voom (and optionally voomWithQualityWeights) pipelines with
#' \code{duplicateCorrelation} blocking and compares DE results across
#' model configurations. This helps decide whether sample quality weights
#' improve the analysis for a given contrast.
#'
#' @param dge       A filtered, TMM-normalized \code{DGEList}
#'                  (e.g., \code{normalise_counts()$dge}).
#' @param design    Design matrix from \code{model.matrix()}.
#' @param contrast  A character string (e.g.,
#'                  \code{"groupTreat.BL - groupCtrl.BL"}) passed to
#'                  \code{limma::makeContrasts()}, or a numeric contrast
#'                  vector of length \code{ncol(design)}.
#' @param block     Factor: blocking variable (e.g., subject ID). Required
#'                  for \code{duplicateCorrelation}.
#' @param include_unblocked Logical; also fit models without blocking for
#'                  comparison (default FALSE).
#' @param fdr_threshold Numeric; FDR cutoff for DE gene counting (default 0.05).
#' @param lfc_threshold Numeric; minimum absolute logFC threshold. When > 0,
#'                  uses \code{limma::treat()} instead of \code{eBayes()}.
#' @param plot      Logical; produce diagnostic comparison plots (default TRUE).
#' @param plot_file File path to save the plot. NULL = not saved to file.
#'
#' @return A named list:
#'   \item{summary}{Tibble with one row per model and key comparison metrics.}
#'   \item{models}{Named list of per-model details (fit objects, topTable, etc.).}
#'   \item{concordance}{Pairwise comparison metrics between models.}
#'   \item{plot}{ggplot object, or NULL if \code{plot = FALSE}.}
compare_voom_models <- function(dge, design, contrast, block,
                                include_unblocked = FALSE,
                                fdr_threshold = 0.05, lfc_threshold = 0,
                                plot = TRUE, plot_file = NULL) {

    # -- Input validation ------------------------------------------------------
    if (!inherits(dge, "DGEList")) {
        stop("`dge` must be a DGEList object")
    }
    if (!is.matrix(design)) {
        stop("`design` must be a matrix (from model.matrix())")
    }
    if (nrow(design) != ncol(dge)) {
        stop("`design` must have one row per sample (expected ", ncol(dge),
             " rows, got ", nrow(design), ")")
    }

    # Process contrast
    if (is.character(contrast) && length(contrast) == 1L) {
        contrast_mat <- tryCatch(
            limma::makeContrasts(contrasts = contrast, levels = design),
            error = function(e) {
                avail <- colnames(design)
                cmap  <- attr(design, "colname_map")
                hint  <- if (length(cmap) > 0L) {
                    paste0("\nNote: some column names were sanitized. ",
                           "Original names: ",
                           paste(cmap, collapse = ", "))
                } else {
                    ""
                }
                stop("Failed to build contrast from '", contrast,
                     "'. Available coefficients: ",
                     paste(avail, collapse = ", "),
                     hint,
                     "\nOriginal error: ", e$message, call. = FALSE)
            }
        )
    } else if (is.numeric(contrast)) {
        if (length(contrast) != ncol(design)) {
            stop("`contrast` vector must have length ", ncol(design),
                 " (one per coefficient), got ", length(contrast))
        }
        contrast_mat <- matrix(contrast, ncol = 1)
        rownames(contrast_mat) <- colnames(design)
    } else {
        stop("`contrast` must be a single character string or numeric vector")
    }

    # Block validation
    if (missing(block) || is.null(block)) {
        stop("`block` is required (experiment has repeated measures)")
    }
    block <- as.factor(block)
    if (length(block) != ncol(dge)) {
        stop("`block` must have length equal to number of samples (",
             ncol(dge), ")")
    }

    # Threshold validation
    if (fdr_threshold <= 0 || fdr_threshold >= 1) {
        stop("`fdr_threshold` must be between 0 and 1 (exclusive)")
    }
    if (lfc_threshold < 0) {
        stop("`lfc_threshold` must be >= 0")
    }

    # -- Internal helper: fit voom model ---------------------------------------
    fit_voom_model <- function(dge, design, block, use_weights) {
        if (!use_weights) {
            # Two-pass voom + duplicateCorrelation
            if (!is.null(block)) {
                v1 <- limma::voom(dge, design = design)
                corfit <- limma::duplicateCorrelation(v1, design = design,
                                                      block = block)
                v2 <- limma::voom(dge, design = design,
                                  block = block,
                                  correlation = corfit$consensus.correlation)
                corfit <- limma::duplicateCorrelation(v2, design = design,
                                                      block = block)
                block_cor <- corfit$consensus.correlation
                fit <- limma::lmFit(v2, design = design,
                                    block = block,
                                    correlation = block_cor)
            } else {
                v2 <- limma::voom(dge, design = design)
                block_cor <- NA_real_
                fit <- limma::lmFit(v2, design = design)
            }
        } else {
            # Two-pass voomWithQualityWeights + duplicateCorrelation
            if (!is.null(block)) {
                vsw1 <- limma::voomWithQualityWeights(dge, design = design)
                corfit <- limma::duplicateCorrelation(vsw1, design = design,
                                                      block = block)
                vsw2 <- limma::voomWithQualityWeights(
                    dge, design = design,
                    block = block,
                    correlation = corfit$consensus.correlation)
                corfit <- limma::duplicateCorrelation(vsw2, design = design,
                                                      block = block)
                block_cor <- corfit$consensus.correlation
                v2 <- vsw2
                fit <- limma::lmFit(v2, design = design,
                                    block = block,
                                    correlation = block_cor)
            } else {
                v2 <- limma::voomWithQualityWeights(dge, design = design)
                block_cor <- NA_real_
                fit <- limma::lmFit(v2, design = design)
            }
        }

        # Warn on unusual block correlations
        if (!is.na(block_cor)) {
            if (block_cor < 0) {
                warning("Negative block correlation (", round(block_cor, 3),
                        ") — this is unusual and may indicate a problem",
                        call. = FALSE)
            } else if (block_cor > 0.95) {
                warning("Very high block correlation (", round(block_cor, 3),
                        ") — near-1 values may indicate the blocking variable",
                        " explains almost all variance", call. = FALSE)
            }
        }

        list(voom_obj = v2, fit = fit, block_correlation = block_cor)
    }

    # -- Internal helper: weighted R-squared -----------------------------------
    compute_weighted_r2 <- function(voom_obj, fit, design) {
        y <- voom_obj$E
        w <- voom_obj$weights
        fitted_vals <- tcrossprod(fit$coefficients, design)
        residuals <- y - fitted_vals
        w_sum <- rowSums(w)
        y_wmean <- rowSums(w * y) / w_sum
        ss_tot <- rowSums(w * (y - y_wmean)^2)
        ss_res <- rowSums(w * residuals^2)
        1 - ss_res / ss_tot
    }

    # -- Internal helper: summarise a model ------------------------------------
    summarise_model <- function(model_result, contrast_mat, fdr, lfc) {
        voom_obj <- model_result$voom_obj
        fit <- model_result$fit
        block_cor <- model_result$block_correlation

        # Contrast fit → moderated stats
        cfit <- limma::contrasts.fit(fit, contrast_mat)
        if (lfc > 0) {
            efit <- limma::treat(cfit, lfc = lfc)
        } else {
            efit <- limma::eBayes(cfit)
        }

        # topTable
        top <- limma::topTable(efit, number = Inf, sort.by = "none")

        # DE decisions
        dt <- limma::decideTests(efit, p.value = fdr, lfc = lfc)
        n_up   <- sum(dt > 0)
        n_down <- sum(dt < 0)
        n_de   <- n_up + n_down
        de_genes <- rownames(dt)[dt != 0]

        # Per-gene fit quality (from pre-contrast fit)
        rmse <- fit$sigma
        r_squared <- compute_weighted_r2(voom_obj, fit, fit$design)

        list(
            voom_obj         = voom_obj,
            fit              = fit,
            cfit             = cfit,
            efit             = efit,
            top_table        = top,
            rmse             = rmse,
            r_squared        = r_squared,
            de_genes         = de_genes,
            n_de             = n_de,
            n_up             = n_up,
            n_down           = n_down,
            median_rmse      = stats::median(rmse),
            median_r_squared = stats::median(r_squared),
            prior_df         = efit$df.prior,
            prior_var        = efit$s2.prior,
            block_correlation = block_cor
        )
    }

    # -- Fit all models --------------------------------------------------------
    models <- list()

    message("Fitting voom + block model...")
    models[["voom_block"]] <- summarise_model(
        fit_voom_model(dge, design, block, use_weights = FALSE),
        contrast_mat, fdr_threshold, lfc_threshold
    )

    message("Fitting voom + sample weights + block model...")
    models[["voom_sw_block"]] <- summarise_model(
        fit_voom_model(dge, design, block, use_weights = TRUE),
        contrast_mat, fdr_threshold, lfc_threshold
    )

    if (include_unblocked) {
        message("Fitting voom model (unblocked)...")
        models[["voom"]] <- summarise_model(
            fit_voom_model(dge, design, block = NULL, use_weights = FALSE),
            contrast_mat, fdr_threshold, lfc_threshold
        )

        message("Fitting voom + sample weights model (unblocked)...")
        models[["voom_sw"]] <- summarise_model(
            fit_voom_model(dge, design, block = NULL, use_weights = TRUE),
            contrast_mat, fdr_threshold, lfc_threshold
        )
    }

    model_names <- names(models)

    # -- Summary tibble --------------------------------------------------------
    summary_df <- tibble::tibble(
        model             = model_names,
        n_de              = vapply(models, `[[`, integer(1), "n_de"),
        n_up              = vapply(models, `[[`, integer(1), "n_up"),
        n_down            = vapply(models, `[[`, integer(1), "n_down"),
        median_rmse       = vapply(models, `[[`, numeric(1), "median_rmse"),
        median_r_squared  = vapply(models, `[[`, numeric(1), "median_r_squared"),
        prior_df          = vapply(models, `[[`, numeric(1), "prior_df"),
        prior_var         = vapply(models, `[[`, numeric(1), "prior_var"),
        block_correlation = vapply(models, `[[`, numeric(1), "block_correlation")
    )

    # -- Cross-model comparisons -----------------------------------------------
    n_models <- length(model_names)

    # Extract logFC and p-value vectors per model
    logfc_list <- lapply(models, function(m) m$top_table$logFC)
    pval_list  <- lapply(models, function(m) m$top_table$P.Value)

    # Pairwise logFC Pearson correlation
    logfc_mat <- do.call(cbind, logfc_list)
    logfc_cor <- stats::cor(logfc_mat, method = "pearson")

    # Pairwise p-value rank Spearman correlation
    rank_mat <- apply(do.call(cbind, pval_list), 2, rank)
    rank_cor <- stats::cor(rank_mat, method = "spearman")

    # Pairwise Jaccard index of DE gene sets
    jaccard_mat <- matrix(NA_real_, nrow = n_models, ncol = n_models,
                          dimnames = list(model_names, model_names))
    for (i in seq_len(n_models)) {
        for (j in seq_len(n_models)) {
            a <- models[[i]]$de_genes
            b <- models[[j]]$de_genes
            union_size <- length(union(a, b))
            if (union_size == 0L) {
                jaccard_mat[i, j] <- NA_real_
            } else {
                jaccard_mat[i, j] <- length(intersect(a, b)) / union_size
            }
        }
    }

    # DE overlap tibble
    all_de <- unique(unlist(lapply(models, `[[`, "de_genes")))
    if (length(all_de) > 0L) {
        overlap_cols <- lapply(models, function(m) all_de %in% m$de_genes)
        de_overlap <- tibble::tibble(gene = all_de)
        for (nm in model_names) {
            de_overlap[[nm]] <- overlap_cols[[nm]]
        }
    } else {
        de_overlap <- tibble::tibble(gene = character(0))
        for (nm in model_names) {
            de_overlap[[nm]] <- logical(0)
        }
    }

    concordance <- list(
        logfc_cor  = logfc_cor,
        rank_cor   = rank_cor,
        jaccard    = jaccard_mat,
        de_overlap = de_overlap
    )

    # -- Plot ------------------------------------------------------------------
    p <- NULL
    if (plot && n_models >= 2L) {
        # Panel A: RMSE density overlay
        rmse_df <- do.call(rbind, lapply(model_names, function(nm) {
            data.frame(model = nm, rmse = models[[nm]]$rmse,
                       stringsAsFactors = FALSE)
        }))
        pa <- ggplot2::ggplot(rmse_df,
                              ggplot2::aes(x = rmse, color = model)) +
            ggplot2::geom_density() +
            ggplot2::theme_minimal() +
            ggplot2::labs(title = "Per-gene RMSE", x = "RMSE", y = "Density")

        # Panel B: R-squared density overlay
        r2_df <- do.call(rbind, lapply(model_names, function(nm) {
            data.frame(model = nm, r_squared = models[[nm]]$r_squared,
                       stringsAsFactors = FALSE)
        }))
        pb <- ggplot2::ggplot(r2_df,
                              ggplot2::aes(x = r_squared, color = model)) +
            ggplot2::geom_density() +
            ggplot2::theme_minimal() +
            ggplot2::labs(title = expression(Weighted ~ R^2),
                         x = expression(R^2), y = "Density")

        # Panel C: DE gene count bar chart
        de_df <- do.call(rbind, lapply(model_names, function(nm) {
            data.frame(
                model     = nm,
                direction = c("Up", "Down"),
                count     = c(models[[nm]]$n_up, -models[[nm]]$n_down),
                stringsAsFactors = FALSE
            )
        }))
        de_df$model <- factor(de_df$model, levels = model_names)
        pc <- ggplot2::ggplot(de_df,
                              ggplot2::aes(x = model, y = count,
                                           fill = direction)) +
            ggplot2::geom_col(position = "identity") +
            ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
            ggplot2::scale_fill_manual(values = c(Up = "#D73027",
                                                   Down = "#4575B4")) +
            ggplot2::theme_minimal() +
            ggplot2::labs(title = "DE genes", x = NULL, y = "Count") +
            ggplot2::theme(axis.text.x = ggplot2::element_text(
                angle = 45, hjust = 1))

        # Panel D: logFC scatter (model 1 vs model 2)
        m1 <- model_names[1]
        m2 <- model_names[2]
        scatter_df <- data.frame(
            logfc_1 = models[[m1]]$top_table$logFC,
            logfc_2 = models[[m2]]$top_table$logFC,
            de_status = ifelse(
                rownames(models[[m1]]$top_table) %in% models[[m1]]$de_genes &
                rownames(models[[m1]]$top_table) %in% models[[m2]]$de_genes,
                "Both",
                ifelse(
                    rownames(models[[m1]]$top_table) %in% models[[m1]]$de_genes,
                    m1,
                    ifelse(
                        rownames(models[[m1]]$top_table) %in% models[[m2]]$de_genes,
                        m2,
                        "None"
                    )
                )
            ),
            stringsAsFactors = FALSE
        )
        pd <- ggplot2::ggplot(scatter_df,
                              ggplot2::aes(x = logfc_1, y = logfc_2,
                                           color = de_status)) +
            ggplot2::geom_point(size = 0.5, alpha = 0.4) +
            ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                                 color = "grey40") +
            ggplot2::scale_color_manual(
                values = c(Both = "#D73027", None = "grey70",
                           stats::setNames("#4575B4", m1),
                           stats::setNames("#FDB863", m2))
            ) +
            ggplot2::theme_minimal() +
            ggplot2::labs(title = "LogFC comparison",
                         x = paste("logFC —", m1),
                         y = paste("logFC —", m2),
                         color = "DE in")

        p <- patchwork::wrap_plots(pa, pb, pc, pd, ncol = 2)

        if (!is.null(plot_file)) {
            ggplot2::ggsave(plot_file, plot = p, width = 12, height = 10)
            message("Saved plot to ", plot_file)
        }
    }

    # -- Return ----------------------------------------------------------------
    message("Done.")
    list(
        summary     = summary_df,
        models      = models,
        concordance = concordance,
        plot        = p
    )
}
