#' Select the best model from compare_voom_models() output
#'
#' Automates model selection by evaluating fit quality metrics (median RMSE,
#' median R-squared) across voom strategies and, optionally, across multiple
#' designs. Selects the best voom strategy within each design, then the best
#' design overall, preferring simpler models when fit improvement is negligible.
#'
#' @param comparisons Either a single \code{compare_voom_models()} output, or a
#'   named list of them (one per design, e.g., from running
#'   \code{compare_voom_models()} on each design from \code{suggest_models()}).
#' @param verbose Logical; print a recommendation report via \code{message()}
#'   (default TRUE).
#'
#' @return A named list:
#'   \item{recommended}{Character: winning model name (e.g., "voom_sw_block").}
#'   \item{design}{Character: winning design name.}
#'   \item{reason}{Character: explanation of the recommendation.}
#'   \item{model}{The winning model's detail list (voom_obj, fit, efit,
#'     top_table, de_genes, etc.).}
#'   \item{design_summary}{Tibble: cross-design comparison with one row per
#'     design, showing the winning voom strategy and its metrics.}
select_model <- function(comparisons, verbose = TRUE) {

    # -- Input handling --------------------------------------------------------
    comparisons <- validate_comparisons(comparisons)

    # -- Stage 1: Within each design, pick the best voom strategy --------------
    design_winners <- lapply(names(comparisons), function(design_name) {
        comp <- comparisons[[design_name]]
        pick_voom_strategy(comp, design_name, verbose = verbose)
    })
    names(design_winners) <- names(comparisons)

    # -- Stage 2: Across designs, pick the best design -------------------------
    best <- pick_best_design(design_winners, comparisons = comparisons,
                              verbose = verbose)

    # -- Build design_summary tibble -------------------------------------------
    design_summary <- tibble::tibble(
        design       = names(design_winners),
        winner       = vapply(design_winners, `[[`, character(1), "winner"),
        median_rmse  = vapply(design_winners, function(w) {
            w$winner_row$median_rmse
        }, numeric(1)),
        median_r_squared = vapply(design_winners, function(w) {
            w$winner_row$median_r_squared
        }, numeric(1)),
        n_de         = vapply(design_winners, function(w) {
            w$winner_row$n_de
        }, integer(1))
    )

    # -- Build reason ----------------------------------------------------------
    reason <- best$reason

    # -- Verbose: cross-design comparison --------------------------------------
    if (verbose && length(comparisons) > 1L) {
        message("\n-- Design comparison --")

        # Header
        message(format_table_row(
            c("Design", "Winner", "median_rmse", "median_r2", "n_de"),
            widths = c(20, 18, 14, 14, 8)
        ))

        for (i in seq_len(nrow(design_summary))) {
            message(format_table_row(
                c(design_summary$design[i],
                  design_summary$winner[i],
                  format_num(design_summary$median_rmse[i]),
                  format_num(design_summary$median_r_squared[i]),
                  as.character(design_summary$n_de[i])),
                widths = c(20, 18, 14, 14, 8)
            ))
        }

        message("\nRecommendation: ", best$design, " / ", best$winner)
        message("Reason: ", reason)
    }

    if (verbose && length(comparisons) == 1L) {
        message("\nRecommendation: ", best$winner)
        message("Reason: ", reason)
    }

    # -- Return ----------------------------------------------------------------
    # Prefer cross-design concordance plot; fall back to within-design
    conc_plot <- best$concordance_plot
    if (is.null(conc_plot)) {
        conc_plot <- design_winners[[best$design]]$concordance_plot
    }

    list(
        recommended      = best$winner,
        design           = best$design,
        reason           = reason,
        model            = comparisons[[best$design]]$models[[best$winner]],
        design_summary   = design_summary,
        concordance_plot = conc_plot
    )
}


# -- Internal helpers ----------------------------------------------------------

#' Validate and normalise the comparisons input
#' @noRd
validate_comparisons <- function(comparisons) {
    is_single_comp <- function(x) {
        is.list(x) &&
            all(c("summary", "models", "concordance") %in% names(x))
    }

    if (is_single_comp(comparisons)) {
        return(list(design_1 = comparisons))
    }

    if (is.list(comparisons) && length(comparisons) >= 1L &&
        !is.null(names(comparisons)) &&
        all(vapply(comparisons, is_single_comp, logical(1)))) {
        return(comparisons)
    }

    stop("Expected a compare_voom_models() output or a named list of them.\n",
         "A single output should have $summary, $models, and $concordance.\n",
         "A multi-design input should be a named list of such outputs.",
         call. = FALSE)
}


#' Pick the best voom strategy within one design
#' @noRd
pick_voom_strategy <- function(comp, design_name, verbose = TRUE) {
    summ <- comp$summary

    # Only consider blocked models for recommendation
    blocked <- summ[!is.na(summ$block_correlation), , drop = FALSE]
    if (nrow(blocked) == 0L) {
        # Fall back to all models if none are blocked (unusual)
        blocked <- summ
    }

    # If only one model, it wins by default
    if (nrow(blocked) == 1L) {
        winner_name <- blocked$model[1L]
        winner_row  <- blocked[1L, ]

        if (verbose) {
            message("\n-- Design: ", design_name, " --")
            message("  Only one model (", winner_name,
                    "); selected by default")
        }

        return(list(
            winner           = winner_name,
            winner_row       = winner_row,
            summary          = summ,
            concordance_plot = NULL
        ))
    }

    # Compare the two blocked models: voom_block vs voom_sw_block
    m1_idx <- which(blocked$model == "voom_block")
    m2_idx <- which(blocked$model == "voom_sw_block")

    if (length(m1_idx) == 0L || length(m2_idx) == 0L) {
        # If expected names not found, compare by position
        m1_idx <- 1L
        m2_idx <- 2L
    }

    m1 <- blocked[m1_idx, ]
    m2 <- blocked[m2_idx, ]

    # Direction winners
    rmse_winner <- if (m1$median_rmse < m2$median_rmse) m1$model
                   else if (m2$median_rmse < m1$median_rmse) m2$model
                   else "tied"
    r2_winner   <- if (m1$median_r_squared > m2$median_r_squared) m1$model
                   else if (m2$median_r_squared > m1$median_r_squared) m2$model
                   else "tied"

    # Relative difference magnitudes (% of better value)
    rmse_min <- min(m1$median_rmse, m2$median_rmse)
    if (rmse_min > 0) {
        rmse_rel_pct <- abs(m1$median_rmse - m2$median_rmse) / rmse_min * 100
    } else {
        rmse_rel_pct <- 0
        rmse_winner  <- "tied"
    }

    r2_max <- max(m1$median_r_squared, m2$median_r_squared)
    if (r2_max > 0) {
        r2_rel_pct <- abs(m1$median_r_squared - m2$median_r_squared) / r2_max * 100
    } else {
        r2_rel_pct <- 0
        r2_winner  <- "tied"
    }

    negligible_pct <- 1.0  # differences below 1% are essentially tied
    concordance_evaluated <- FALSE
    concordance_plot <- NULL

    # Decision logic — magnitude-aware
    if (rmse_winner == r2_winner && rmse_winner != "tied") {
        winner_name <- rmse_winner
        decision_reason <- "both metrics agree"
    } else if (rmse_rel_pct < negligible_pct && r2_rel_pct < negligible_pct) {
        # Fit metrics are comparable — check if DE gene counts diverge
        concordance_evaluated <- TRUE
        m1_name <- m1$model
        m2_name <- m2$model

        if (m2$n_de > m1$n_de) {
            more_de_model  <- m2_name
            fewer_de_model <- m1_name
            more_de_n      <- m2$n_de
            fewer_de_n     <- m1$n_de
        } else {
            more_de_model  <- m1_name
            fewer_de_model <- m2_name
            more_de_n      <- m1$n_de
            fewer_de_n     <- m2$n_de
        }

        conc <- check_concordance_power(
            fewer_n_de     = fewer_de_n,
            more_n_de      = more_de_n,
            fewer_prior_df = blocked[blocked$model == fewer_de_model, ]$prior_df,
            more_prior_df  = blocked[blocked$model == more_de_model, ]$prior_df,
            fewer_de_genes = comp$models[[fewer_de_model]]$de_genes,
            more_de_genes  = comp$models[[more_de_model]]$de_genes,
            logfc_cor_val  = comp$concordance$logfc_cor[fewer_de_model,
                                                        more_de_model],
            fewer_label    = fewer_de_model,
            more_label     = more_de_model
        )

        if (conc$trustworthy) {
            winner_name <- more_de_model
            decision_reason <- sprintf(
                paste0("fit comparable (<%.0f%%); %s finds more DE genes ",
                       "(%d vs %d) with good concordance; %s"),
                negligible_pct, more_de_model, more_de_n, fewer_de_n,
                conc$detail)
        } else {
            winner_name <- m1$model
            decision_reason <- sprintf(
                paste0("fit comparable (<%.0f%%); %s"),
                negligible_pct, conc$detail)
        }

        if (verbose && conc$substantial) {
            message("\n  Concordance diagnostics (comparable fit, DE gene check):")
            message(sprintf("    DE genes: %s=%d, %s=%d (ratio: %.1fx)",
                            fewer_de_model, fewer_de_n,
                            more_de_model, more_de_n,
                            more_de_n / max(fewer_de_n, 1L)))
            message(sprintf("    Subset overlap: %.0f%% of %s DE genes in %s",
                            conc$subset_frac * 100,
                            fewer_de_model, more_de_model))
            message(sprintf("    LogFC correlation: %.4f", conc$logfc_cor))
            message(sprintf("    Prior df ratio (%s/%s): %.2f",
                            more_de_model, fewer_de_model,
                            conc$prior_df_ratio))
            message(sprintf("    Verdict: %s",
                            if (conc$trustworthy) "power gain supported"
                            else "insufficient concordance"))
        }

        # Build concordance diagnostic plot when DE counts diverge substantially
        if (conc$substantial) {
            concordance_plot <- build_concordance_plot(
                comp, fewer_de_model, more_de_model, conc)
        }
    } else if (rmse_winner != "tied" && r2_winner != "tied" &&
               rmse_winner != r2_winner) {
        # Disagree — whichever metric has the larger relative difference wins
        if (r2_rel_pct > rmse_rel_pct) {
            winner_name <- r2_winner
            decision_reason <- sprintf(
                "metrics disagree: R2 difference (%.1f%%) outweighs RMSE difference (%.1f%%)",
                r2_rel_pct, rmse_rel_pct)
        } else if (rmse_rel_pct > r2_rel_pct) {
            winner_name <- rmse_winner
            decision_reason <- sprintf(
                "metrics disagree: RMSE difference (%.1f%%) outweighs R2 difference (%.1f%%)",
                rmse_rel_pct, r2_rel_pct)
        } else {
            winner_name <- m1$model
            decision_reason <- sprintf(
                "metrics disagree with equal relative differences (%.1f%%); preferring simpler model",
                rmse_rel_pct)
        }
    } else if (rmse_winner != "tied") {
        winner_name <- rmse_winner
        decision_reason <- "R2 tied; decided by RMSE"
    } else if (r2_winner != "tied") {
        winner_name <- r2_winner
        decision_reason <- "RMSE tied; decided by R2"
    } else {
        winner_name <- m1$model
        decision_reason <- "both metrics tied; preferring simpler model"
    }

    winner_row <- blocked[blocked$model == winner_name, ]

    # Suspicious pattern: worse fit but more DE genes
    # (skip when concordance was already evaluated for this scenario)
    if (!concordance_evaluated) {
        loser_name <- setdiff(c(m1$model, m2$model), winner_name)
        loser_row  <- blocked[blocked$model == loser_name, ]
        if (nrow(loser_row) == 1L && loser_row$n_de > winner_row$n_de &&
            loser_row$median_rmse > winner_row$median_rmse) {
            warning("Design '", design_name, "': ", loser_name,
                    " has worse fit (higher RMSE) but more DE genes (",
                    loser_row$n_de, " vs ", winner_row$n_de,
                    "). This may indicate spurious discoveries.",
                    call. = FALSE)
        }
    }

    # Verbose output
    if (verbose) {
        print_voom_comparison(design_name, m1, m2, winner_name, decision_reason)
    }

    list(
        winner           = winner_name,
        winner_row       = winner_row,
        summary          = summ,
        concordance_plot = concordance_plot
    )
}


#' Check whether a model with more DE genes represents a trustworthy power gain
#'
#' Works for both within-design (pre-computed concordance via \code{comp}) and
#' cross-design comparisons (raw model data via \code{fewer_model_data} /
#' \code{more_model_data}).
#' @noRd
check_concordance_power <- function(fewer_n_de, more_n_de,
                                     fewer_prior_df, more_prior_df,
                                     fewer_de_genes, more_de_genes,
                                     logfc_cor_val,
                                     fewer_label = "fewer",
                                     more_label = "more") {

    de_ratio <- more_n_de / max(fewer_n_de, 1L)

    # Gate: DE gene difference not substantial enough to evaluate
    if (de_ratio < 1.5) {
        return(list(
            trustworthy    = FALSE,
            substantial    = FALSE,
            subset_frac    = NA_real_,
            logfc_cor      = NA_real_,
            prior_df_ratio = NA_real_,
            detail         = sprintf(
                "DE ratio %.1fx (< 1.5x threshold); difference not substantial",
                de_ratio)
        ))
    }

    # Subset fraction: what fraction of the fewer model's DE genes are in the
    # more model's DE set?
    if (fewer_n_de > 0L) {
        n_both <- length(intersect(fewer_de_genes, more_de_genes))
        subset_frac <- n_both / fewer_n_de
    } else {
        subset_frac <- 1.0  # vacuously true when fewer model finds nothing
    }

    # Prior df ratio
    prior_df_ratio <- more_prior_df / fewer_prior_df

    # Thresholds
    subset_ok   <- subset_frac >= 0.80
    logfc_ok    <- !is.na(logfc_cor_val) && logfc_cor_val >= 0.95
    prior_df_ok <- prior_df_ratio >= 0.5

    trustworthy <- subset_ok && logfc_ok && prior_df_ok

    # Build detail string
    checks <- c(
        sprintf("subset=%.0f%%%s", subset_frac * 100,
                if (subset_ok) "" else " [FAIL]"),
        sprintf("logFC_cor=%.3f%s", logfc_cor_val,
                if (logfc_ok) "" else " [FAIL]"),
        sprintf("prior_df_ratio=%.2f%s", prior_df_ratio,
                if (prior_df_ok) "" else " [FAIL]")
    )
    verdict <- if (trustworthy) "power gain supported" else "insufficient concordance"
    detail <- paste0(paste(checks, collapse = ", "), " -> ", verdict)

    list(
        trustworthy    = trustworthy,
        substantial    = TRUE,
        subset_frac    = subset_frac,
        logfc_cor      = logfc_cor_val,
        prior_df_ratio = prior_df_ratio,
        detail         = detail
    )
}


#' Build concordance diagnostic plot (logFC scatter + DE overlap bar)
#' @noRd
build_concordance_plot <- function(comp, fewer_model, more_model, conc) {
    fewer_top <- comp$models[[fewer_model]]$top_table
    more_top  <- comp$models[[more_model]]$top_table
    fewer_de  <- comp$models[[fewer_model]]$de_genes
    more_de   <- comp$models[[more_model]]$de_genes

    genes <- rownames(fewer_top)

    # -- Panel A: LogFC scatter with DE membership colouring -------------------
    de_status <- ifelse(
        genes %in% fewer_de & genes %in% more_de, "Both",
        ifelse(genes %in% fewer_de, fewer_model,
               ifelse(genes %in% more_de, more_model, "Neither")))

    scatter_df <- data.frame(
        logfc_fewer = fewer_top$logFC,
        logfc_more  = more_top$logFC,
        de_status   = de_status,
        stringsAsFactors = FALSE
    )
    # Plot non-DE genes first, DE genes on top
    scatter_df$de_status <- factor(scatter_df$de_status,
        levels = c("Neither", fewer_model, more_model, "Both"))
    scatter_df <- scatter_df[order(scatter_df$de_status), ]

    pa <- ggplot2::ggplot(scatter_df,
                          ggplot2::aes(x = logfc_fewer, y = logfc_more,
                                       color = de_status)) +
        ggplot2::geom_point(size = 0.5, alpha = 0.4) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                             color = "grey40") +
        ggplot2::scale_color_manual(
            values = c(Neither = "grey80",
                       stats::setNames("#4575B4", fewer_model),
                       stats::setNames("#FDB863", more_model),
                       Both = "#D73027"),
            drop = FALSE
        ) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
            title = "LogFC concordance",
            x = paste("logFC \u2014", fewer_model),
            y = paste("logFC \u2014", more_model),
            color = "DE in"
        )

    # -- Panel B: DE overlap bar chart -----------------------------------------
    n_both       <- sum(genes %in% fewer_de & genes %in% more_de)
    n_only_fewer <- sum(genes %in% fewer_de & !(genes %in% more_de))
    n_only_more  <- sum(!(genes %in% fewer_de) & genes %in% more_de)

    bar_df <- data.frame(
        category = c("Shared", paste("Only", fewer_model),
                     paste("Only", more_model)),
        count    = c(n_both, n_only_fewer, n_only_more),
        stringsAsFactors = FALSE
    )
    bar_df$category <- factor(bar_df$category, levels = bar_df$category)

    subset_label <- if (!is.na(conc$subset_frac)) {
        sprintf("Subset overlap: %.0f%%", conc$subset_frac * 100)
    } else {
        ""
    }

    pb <- ggplot2::ggplot(bar_df, ggplot2::aes(x = category, y = count,
                                                fill = category)) +
        ggplot2::geom_col(show.legend = FALSE) +
        ggplot2::scale_fill_manual(
            values = c(Shared = "#D73027",
                       stats::setNames("#4575B4",
                                       paste("Only", fewer_model)),
                       stats::setNames("#FDB863",
                                       paste("Only", more_model)))
        ) +
        ggplot2::geom_text(ggplot2::aes(label = count), vjust = -0.3,
                           size = 3.5) +
        ggplot2::annotate("text", x = 2, y = max(bar_df$count) * 0.9,
                          label = subset_label, size = 3.5,
                          fontface = "italic") +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = "DE gene overlap", x = NULL, y = "Genes") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30,
                                                            hjust = 1))

    patchwork::wrap_plots(pa, pb, ncol = 2, widths = c(1.5, 1))
}


#' Pick the best design across multiple designs
#' @noRd
pick_best_design <- function(design_winners, comparisons = NULL,
                              verbose = TRUE) {
    n_designs <- length(design_winners)

    if (n_designs == 1L) {
        dw <- design_winners[[1L]]
        return(list(
            design          = names(design_winners)[1L],
            winner          = dw$winner,
            reason          = paste0(dw$winner, " selected as best voom strategy"),
            concordance_plot = NULL
        ))
    }

    # Collect metrics for each design's winner
    design_names <- names(design_winners)
    rmses <- vapply(design_winners, function(w) {
        w$winner_row$median_rmse
    }, numeric(1))
    n_des <- vapply(design_winners, function(w) {
        w$winner_row$n_de
    }, integer(1))

    # Rank by RMSE (lower = better)
    best_idx    <- which.min(rmses)
    best_rmse   <- rmses[best_idx]
    best_design <- design_names[best_idx]

    # Check parsimony: does a simpler design have nearly equal RMSE?
    # Design ordering from suggest_models(): "base", "base_plus_1", ...
    # Earlier in the list = simpler
    threshold <- 0.01  # 1% relative difference

    selected_idx <- best_idx
    for (i in seq_along(design_names)) {
        if (i >= best_idx) break
        rel_diff <- (rmses[i] - best_rmse) / best_rmse
        if (rel_diff <= threshold) {
            selected_idx <- i
            break
        }
    }

    # Concordance check: if parsimony selected a simpler design, check whether
    # any more complex design with comparable RMSE has substantially more DE
    # genes with good concordance.
    cross_conc_plot <- NULL
    if (selected_idx < n_designs && !is.null(comparisons)) {
        sel_winner <- design_winners[[selected_idx]]$winner
        sel_n_de   <- n_des[selected_idx]
        sel_prior  <- design_winners[[selected_idx]]$winner_row$prior_df
        sel_model  <- comparisons[[design_names[selected_idx]]]$models[[sel_winner]]

        for (j in (selected_idx + 1L):n_designs) {
            # Only consider designs within the RMSE parsimony window
            rel_diff_j <- (rmses[j] - best_rmse) / best_rmse
            if (rel_diff_j > threshold) next

            cand_winner <- design_winners[[j]]$winner
            cand_n_de   <- n_des[j]
            cand_prior  <- design_winners[[j]]$winner_row$prior_df
            cand_model  <- comparisons[[design_names[j]]]$models[[cand_winner]]

            # Compute cross-design logFC correlation on the fly
            logfc_cor_val <- stats::cor(sel_model$top_table$logFC,
                                        cand_model$top_table$logFC,
                                        method = "pearson")

            # Determine which has more DE genes
            if (cand_n_de > sel_n_de) {
                fewer_de_n  <- sel_n_de
                more_de_n   <- cand_n_de
                fewer_genes <- sel_model$de_genes
                more_genes  <- cand_model$de_genes
                fewer_prior <- sel_prior
                more_prior  <- cand_prior
                fewer_label <- paste0(design_names[selected_idx], "/", sel_winner)
                more_label  <- paste0(design_names[j], "/", cand_winner)
                override_idx <- j
            } else {
                next  # candidate has fewer or equal DE genes — no reason to override
            }

            conc <- check_concordance_power(
                fewer_n_de     = fewer_de_n,
                more_n_de      = more_de_n,
                fewer_prior_df = fewer_prior,
                more_prior_df  = more_prior,
                fewer_de_genes = fewer_genes,
                more_de_genes  = more_genes,
                logfc_cor_val  = logfc_cor_val,
                fewer_label    = fewer_label,
                more_label     = more_label
            )

            if (verbose && conc$substantial) {
                message("\n  Cross-design concordance (",
                        fewer_label, " vs ", more_label, "):")
                message(sprintf("    DE genes: %s=%d, %s=%d (ratio: %.1fx)",
                                fewer_label, fewer_de_n,
                                more_label, more_de_n,
                                more_de_n / max(fewer_de_n, 1L)))
                message(sprintf("    Subset overlap: %.0f%% of %s DE genes in %s",
                                conc$subset_frac * 100, fewer_label, more_label))
                message(sprintf("    LogFC correlation: %.4f", conc$logfc_cor))
                message(sprintf("    Prior df ratio (%s/%s): %.2f",
                                more_label, fewer_label, conc$prior_df_ratio))
                message(sprintf("    Verdict: %s",
                                if (conc$trustworthy) "power gain supported"
                                else "insufficient concordance"))
            }

            if (conc$trustworthy) {
                selected_idx <- override_idx

                # Build cross-design concordance plot
                cross_conc_plot <- build_cross_design_concordance_plot(
                    sel_model, cand_model, fewer_label, more_label, conc)
                break
            }
        }
    }

    selected_design <- design_names[selected_idx]
    selected_winner <- design_winners[[selected_idx]]$winner
    selected_rmse   <- rmses[selected_idx]

    # Build reason
    if (selected_idx == best_idx) {
        if (n_designs > 1L) {
            # The best RMSE design was selected
            second_best_rmse <- sort(rmses)[2L]
            pct_improvement  <- (second_best_rmse - best_rmse) /
                second_best_rmse * 100
            reason <- paste0(
                selected_design, " has lower median RMSE (",
                format_num(best_rmse), " vs ",
                format_num(second_best_rmse), ", -",
                format_num(pct_improvement), "% improvement)")
        } else {
            reason <- paste0(selected_winner,
                             " selected as best voom strategy")
        }
    } else if (!is.null(cross_conc_plot)) {
        # Concordance override: more complex design selected for DE power
        reason <- paste0(
            selected_design, " selected: comparable RMSE to simpler design ",
            "but more DE genes with good concordance")
    } else {
        # Simpler design preferred via parsimony
        rel_diff_pct <- (selected_rmse - best_rmse) / best_rmse * 100
        reason <- paste0(
            "Simpler design '", selected_design,
            "' preferred (median RMSE within ",
            format_num(rel_diff_pct), "% of best design '",
            best_design, "')")
    }

    list(
        design           = selected_design,
        winner           = selected_winner,
        reason           = reason,
        concordance_plot = cross_conc_plot
    )
}


#' Build cross-design concordance plot from two winning models' data
#' @noRd
build_cross_design_concordance_plot <- function(fewer_model_data,
                                                 more_model_data,
                                                 fewer_label, more_label,
                                                 conc) {
    fewer_top <- fewer_model_data$top_table
    more_top  <- more_model_data$top_table
    fewer_de  <- fewer_model_data$de_genes
    more_de   <- more_model_data$de_genes

    genes <- rownames(fewer_top)

    # -- Panel A: LogFC scatter ------------------------------------------------
    de_status <- ifelse(
        genes %in% fewer_de & genes %in% more_de, "Both",
        ifelse(genes %in% fewer_de, fewer_label,
               ifelse(genes %in% more_de, more_label, "Neither")))

    scatter_df <- data.frame(
        logfc_fewer = fewer_top$logFC,
        logfc_more  = more_top$logFC,
        de_status   = de_status,
        stringsAsFactors = FALSE
    )
    scatter_df$de_status <- factor(scatter_df$de_status,
        levels = c("Neither", fewer_label, more_label, "Both"))
    scatter_df <- scatter_df[order(scatter_df$de_status), ]

    pa <- ggplot2::ggplot(scatter_df,
                          ggplot2::aes(x = logfc_fewer, y = logfc_more,
                                       color = de_status)) +
        ggplot2::geom_point(size = 0.5, alpha = 0.4) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                             color = "grey40") +
        ggplot2::scale_color_manual(
            values = c(Neither = "grey80",
                       stats::setNames("#4575B4", fewer_label),
                       stats::setNames("#FDB863", more_label),
                       Both = "#D73027"),
            drop = FALSE
        ) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
            title = "Cross-design logFC concordance",
            x = paste("logFC \u2014", fewer_label),
            y = paste("logFC \u2014", more_label),
            color = "DE in"
        )

    # -- Panel B: DE overlap bar -----------------------------------------------
    n_both       <- sum(genes %in% fewer_de & genes %in% more_de)
    n_only_fewer <- sum(genes %in% fewer_de & !(genes %in% more_de))
    n_only_more  <- sum(!(genes %in% fewer_de) & genes %in% more_de)

    bar_df <- data.frame(
        category = c("Shared", paste("Only", fewer_label),
                     paste("Only", more_label)),
        count    = c(n_both, n_only_fewer, n_only_more),
        stringsAsFactors = FALSE
    )
    bar_df$category <- factor(bar_df$category, levels = bar_df$category)

    subset_label <- if (!is.na(conc$subset_frac)) {
        sprintf("Subset overlap: %.0f%%", conc$subset_frac * 100)
    } else {
        ""
    }

    pb <- ggplot2::ggplot(bar_df, ggplot2::aes(x = category, y = count,
                                                fill = category)) +
        ggplot2::geom_col(show.legend = FALSE) +
        ggplot2::scale_fill_manual(
            values = c(Shared = "#D73027",
                       stats::setNames("#4575B4",
                                       paste("Only", fewer_label)),
                       stats::setNames("#FDB863",
                                       paste("Only", more_label)))
        ) +
        ggplot2::geom_text(ggplot2::aes(label = count), vjust = -0.3,
                           size = 3.5) +
        ggplot2::annotate("text", x = 2, y = max(bar_df$count) * 0.9,
                          label = subset_label, size = 3.5,
                          fontface = "italic") +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = "DE gene overlap (cross-design)",
                      x = NULL, y = "Genes") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30,
                                                            hjust = 1))

    patchwork::wrap_plots(pa, pb, ncol = 2, widths = c(1.5, 1))
}


#' Print the voom strategy comparison table for one design
#' @noRd
print_voom_comparison <- function(design_name, m1, m2, winner_name,
                                  decision_reason = NULL) {
    message("\n-- Design: ", design_name, " --")

    metrics <- c("median_rmse", "median_r_squared", "n_de", "n_up", "n_down",
                 "prior_df", "block_correlation")

    # "lower is better" metrics vs "higher is better"
    lower_better <- c("median_rmse")
    higher_better <- c("median_r_squared", "prior_df")

    # Header
    message(format_table_row(
        c("Metric", m1$model, m2$model, "rel_diff"),
        widths = c(20, 18, 18, 12)
    ))

    for (metric in metrics) {
        v1 <- m1[[metric]]
        v2 <- m2[[metric]]

        # Determine which is better
        star1 <- ""
        star2 <- ""
        if (metric %in% lower_better) {
            if (v1 < v2) star1 <- " *" else if (v2 < v1) star2 <- " *"
        } else if (metric %in% higher_better) {
            if (v1 > v2) star1 <- " *" else if (v2 > v1) star2 <- " *"
        }

        if (is.integer(v1)) {
            s1 <- as.character(v1)
            s2 <- as.character(v2)
        } else {
            s1 <- format_num(v1)
            s2 <- format_num(v2)
        }

        # Relative difference for the two key metrics
        rel_str <- ""
        if (metric == "median_rmse") {
            denom <- min(v1, v2)
            if (denom > 0) {
                rel_str <- paste0(format_num(abs(v1 - v2) / denom * 100), "%")
            }
        } else if (metric == "median_r_squared") {
            denom <- max(v1, v2)
            if (denom > 0) {
                rel_str <- paste0(format_num(abs(v1 - v2) / denom * 100), "%")
            }
        }

        message(format_table_row(
            c(metric, paste0(s1, star1), paste0(s2, star2), rel_str),
            widths = c(20, 18, 18, 12)
        ))
    }

    message("  * = better on this metric")
    message("  Winner: ", winner_name)
    if (!is.null(decision_reason)) {
        message("  Reason: ", decision_reason)
    }
}


#' Format a number for display
#' @noRd
format_num <- function(x) {
    if (is.na(x)) return("NA")
    if (abs(x) >= 100) return(formatC(x, format = "f", digits = 1))
    formatC(x, format = "f", digits = 4)
}


#' Format a table row with fixed-width columns
#' @noRd
format_table_row <- function(values, widths) {
    parts <- vapply(seq_along(values), function(i) {
        formatC(values[i], width = widths[i], flag = "-")
    }, character(1))
    paste0("  ", paste(parts, collapse = "  "))
}
