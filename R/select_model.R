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
    best <- pick_best_design(design_winners, verbose = verbose)

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
    list(
        recommended    = best$winner,
        design         = best$design,
        reason         = reason,
        model          = comparisons[[best$design]]$models[[best$winner]],
        design_summary = design_summary
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
            winner     = winner_name,
            winner_row = winner_row,
            summary    = summ
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

    # Decision logic — magnitude-aware
    if (rmse_winner == r2_winner && rmse_winner != "tied") {
        winner_name <- rmse_winner
        decision_reason <- "both metrics agree"
    } else if (rmse_rel_pct < negligible_pct && r2_rel_pct < negligible_pct) {
        winner_name <- m1$model
        decision_reason <- sprintf(
            "metrics disagree but both differences negligible (<%.0f%%); preferring simpler model",
            negligible_pct)
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

    # Verbose output
    if (verbose) {
        print_voom_comparison(design_name, m1, m2, winner_name, decision_reason)
    }

    list(
        winner     = winner_name,
        winner_row = winner_row,
        summary    = summ
    )
}


#' Pick the best design across multiple designs
#' @noRd
pick_best_design <- function(design_winners, verbose = TRUE) {
    n_designs <- length(design_winners)

    if (n_designs == 1L) {
        dw <- design_winners[[1L]]
        return(list(
            design = names(design_winners)[1L],
            winner = dw$winner,
            reason = paste0(dw$winner, " selected as best voom strategy")
        ))
    }

    # Collect metrics for each design's winner
    design_names <- names(design_winners)
    rmses <- vapply(design_winners, function(w) {
        w$winner_row$median_rmse
    }, numeric(1))

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
        design = selected_design,
        winner = selected_winner,
        reason = reason
    )
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
