#' Suggest stepwise design matrices from PC association results
#'
#' Bridges \code{test_pc_associations()} and \code{compare_voom_models()} by
#' ranking metadata variables by their weighted association with principal
#' components and building design matrices with incrementally added covariates.
#'
#' @param pc_assoc   Output list from \code{test_pc_associations()} (uses
#'                   \code{$results} tibble).
#' @param base_formula Formula with the variable of interest, e.g.,
#'                   \code{~ 0 + group}.
#' @param metadata   Same data.frame used in \code{test_pc_associations()}.
#' @param block      Character; blocking variable column name.
#' @param pca        Optional \code{prcomp} object for variance-explained
#'                   weighting. If NULL, uniform \code{1/k} weights are used.
#' @param p_threshold Adjusted p-value threshold for significance (default 0.05).
#' @param max_covariates Maximum number of covariates to add stepwise (default 3).
#' @param combinations  Logical; if TRUE, test all covariate subsets up to
#'                      \code{max_covariates} instead of building a strict
#'                      sequential ladder.  Default FALSE preserves backward
#'                      compatibility.
#'
#' @return A named list:
#'   \item{covariate_ranking}{Tibble of candidate covariates ranked by weighted
#'         association score.}
#'   \item{models}{List of model specifications, each containing name, formula,
#'         design matrix, and covariates added.}
#'   \item{skipped_covariates}{Tibble of covariates skipped due to rank
#'         deficiency.}
#'   \item{sample_weights_flag}{Logical; TRUE if sample weights are
#'         significantly associated with any PC.}
#'   \item{sample_weights_vars}{Character vector of variables significantly
#'         associated with sample weights.}
suggest_models <- function(pc_assoc, base_formula, metadata, block,
                           pca = NULL, p_threshold = 0.05,
                           max_covariates = 3, combinations = FALSE) {

    # -- Input validation ------------------------------------------------------
    if (!is.list(pc_assoc) || is.null(pc_assoc$results)) {
        stop("`pc_assoc` must be the output of test_pc_associations() ",
             "(a list with a $results tibble)")
    }
    results <- pc_assoc$results
    required_cols <- c("variable", "response", "type", "p_adj", "r_squared")
    missing_cols <- setdiff(required_cols, colnames(results))
    if (length(missing_cols) > 0L) {
        stop("`pc_assoc$results` is missing required columns: ",
             paste(missing_cols, collapse = ", "))
    }

    if (!inherits(base_formula, "formula")) {
        stop("`base_formula` must be a formula (e.g., ~ 0 + group)")
    }

    if (!is.data.frame(metadata)) {
        stop("`metadata` must be a data.frame or tibble")
    }

    if (missing(block) || is.null(block)) {
        stop("`block` is required (experiment has repeated measures)")
    }
    if (length(block) != 1L || !is.character(block)) {
        stop("`block` must be a single character string")
    }
    if (!block %in% colnames(metadata)) {
        stop("`block` column '", block, "' not found in metadata")
    }

    if (!is.logical(combinations) || length(combinations) != 1L) {
        stop("`combinations` must be TRUE or FALSE")
    }

    # Check formula variables exist in metadata
    formula_vars <- all.vars(base_formula)
    missing_vars <- setdiff(formula_vars, colnames(metadata))
    if (length(missing_vars) > 0L) {
        stop("Formula variables not found in metadata: ",
             paste(missing_vars, collapse = ", "))
    }

    if (!is.null(pca)) {
        if (!inherits(pca, "prcomp")) {
            stop("`pca` must be a prcomp object")
        }
        if (nrow(pca$x) != nrow(metadata)) {
            stop("`pca` has ", nrow(pca$x), " samples but `metadata` has ",
                 nrow(metadata), " rows")
        }
    }

    # -- Identify candidate covariates -----------------------------------------
    formula_terms <- formula_vars
    candidates <- setdiff(
        unique(results$variable),
        c(formula_terms, block)
    )

    if (length(candidates) == 0L) {
        message("No candidate covariates found beyond formula terms and block.")
        base_design <- model.matrix(base_formula, data = metadata)
        return(list(
            covariate_ranking  = tibble::tibble(
                variable = character(0), score = numeric(0),
                type = character(0), significant_pcs = character(0),
                n_significant_pcs = integer(0), max_r_squared = numeric(0)
            ),
            models = list(
                list(name = "base", formula = base_formula,
                     design = base_design, covariates = character(0))
            ),
            skipped_covariates  = tibble::tibble(
                variable = character(0), score = numeric(0),
                reason = character(0)
            ),
            sample_weights_flag = FALSE,
            sample_weights_vars = character(0)
        ))
    }

    # -- Compute PC variance weights -------------------------------------------
    pc_responses <- grep("^PC\\d+", unique(results$response), value = TRUE)
    n_pcs <- length(pc_responses)

    if (!is.null(pca)) {
        all_var <- pca$sdev^2 / sum(pca$sdev^2)
        # Match to the PCs present in results
        pc_indices <- as.integer(sub("^PC", "", pc_responses))
        var_weights <- all_var[pc_indices]
    } else {
        var_weights <- rep(1 / n_pcs, n_pcs)
    }
    # Normalise to sum to 1
    var_weights <- var_weights / sum(var_weights)
    names(var_weights) <- pc_responses

    # -- Score each candidate --------------------------------------------------
    pc_results <- results[grepl("^PC\\d+", results$response), , drop = FALSE]
    sig_pc <- pc_results[!is.na(pc_results$p_adj) &
                         pc_results$p_adj < p_threshold, , drop = FALSE]

    ranking_rows <- lapply(candidates, function(v) {
        v_sig <- sig_pc[sig_pc$variable == v, , drop = FALSE]
        if (nrow(v_sig) == 0L) {
            return(data.frame(
                variable = v,
                score = 0,
                type = unique(results$type[results$variable == v])[1],
                significant_pcs = "",
                n_significant_pcs = 0L,
                max_r_squared = 0,
                stringsAsFactors = FALSE
            ))
        }
        pcs_hit <- v_sig$response
        r2_vals <- v_sig$r_squared
        weights_hit <- var_weights[pcs_hit]
        score <- sum(r2_vals * weights_hit, na.rm = TRUE)

        data.frame(
            variable = v,
            score = score,
            type = unique(v_sig$type)[1],
            significant_pcs = paste(pcs_hit, collapse = ", "),
            n_significant_pcs = length(pcs_hit),
            max_r_squared = max(r2_vals, na.rm = TRUE),
            stringsAsFactors = FALSE
        )
    })
    ranking <- do.call(rbind, ranking_rows)
    ranking <- ranking[order(-ranking$score), , drop = FALSE]
    covariate_ranking <- tibble::as_tibble(ranking)

    # Keep only candidates with score > 0 for stepwise addition
    sig_candidates <- covariate_ranking$variable[covariate_ranking$score > 0]

    # -- Check sample weights flag ---------------------------------------------
    sw_results <- results[results$response == "sample_weights", , drop = FALSE]
    sw_sig <- sw_results[!is.na(sw_results$p_adj) &
                         sw_results$p_adj < p_threshold, , drop = FALSE]
    sample_weights_flag <- nrow(sw_sig) > 0L
    sample_weights_vars <- sw_sig$variable

    if (sample_weights_flag) {
        message("Sample weights significantly associated with: ",
                paste(sample_weights_vars, collapse = ", "),
                " — consider using voomWithQualityWeights")
    }

    # -- Stepwise design matrix construction -----------------------------------
    models <- list()
    skipped <- list()

    # Base model
    base_design <- model.matrix(base_formula, data = metadata)
    n_samples <- nrow(metadata)

    if (nrow(base_design) != n_samples) {
        warning("Base design matrix has ", nrow(base_design), " rows but ",
                "metadata has ", n_samples,
                " — NAs in formula variables may cause row drops")
    }

    models[[1L]] <- list(
        name = "base",
        formula = base_formula,
        design = base_design,
        covariates = character(0)
    )

    if (length(sig_candidates) == 0L) {
        message("No significantly associated covariates — returning base model only.")
        return(list(
            covariate_ranking   = covariate_ranking,
            models              = models,
            skipped_covariates  = tibble::tibble(
                variable = character(0), score = numeric(0),
                reason = character(0)
            ),
            sample_weights_flag = sample_weights_flag,
            sample_weights_vars = sample_weights_vars
        ))
    }

    if (combinations) {
        # -- Combinatorial branch --------------------------------------------------
        # Pre-filter: keep covariates that are individually valid
        valid_covariates <- character(0)
        for (cov in sig_candidates) {
            if (length(valid_covariates) >= max_covariates) break
            result <- try_design(base_formula, cov, metadata, n_samples)
            if (result$valid) {
                valid_covariates <- c(valid_covariates, cov)
            } else {
                message("Covariate '", cov, "' invalid individually: ",
                        result$reason, " — excluding from combinations")
                skipped[[length(skipped) + 1L]] <- data.frame(
                    variable = cov,
                    score = covariate_ranking$score[
                        covariate_ranking$variable == cov],
                    reason = result$reason,
                    stringsAsFactors = FALSE
                )
            }
        }

        # Generate all subsets, ordered by size (1, 2, ..., k)
        if (length(valid_covariates) > 0L) {
            for (size in seq_len(length(valid_covariates))) {
                combos <- utils::combn(valid_covariates, size,
                                       simplify = FALSE)
                for (covset in combos) {
                    result <- try_design(base_formula, covset,
                                         metadata, n_samples)
                    if (result$valid) {
                        model_name <- paste(covset, collapse = "_")
                        models[[length(models) + 1L]] <- list(
                            name = model_name,
                            formula = result$formula,
                            design = result$design,
                            covariates = covset
                        )
                        message("Valid combination: '", model_name,
                                "' (", ncol(result$design), " columns)")
                    } else {
                        skipped[[length(skipped) + 1L]] <- data.frame(
                            variable = paste(covset, collapse = " + "),
                            score = NA_real_,
                            reason = result$reason,
                            stringsAsFactors = FALSE
                        )
                        message("Combination '",
                                paste(covset, collapse = " + "),
                                "' skipped: ", result$reason)
                    }
                }
            }
        }

        if (length(models) == 1L) {
            warning("All covariate combinations caused rank deficiency — ",
                    "returning base model only")
        }
    } else {
        # -- Sequential ladder (original behaviour) --------------------------------
        current_covariates <- character(0)
        n_added <- 0L

        for (cov in sig_candidates) {
            if (n_added >= max_covariates) break

            trial_covariates <- c(current_covariates, cov)
            result <- try_design(base_formula, trial_covariates,
                                 metadata, n_samples)

            if (!result$valid) {
                message("Adding '", cov, "' failed: ", result$reason,
                        " — skipping")
                skipped[[length(skipped) + 1L]] <- data.frame(
                    variable = cov,
                    score = covariate_ranking$score[
                        covariate_ranking$variable == cov],
                    reason = result$reason,
                    stringsAsFactors = FALSE
                )
                next
            }

            # Valid — add this covariate
            current_covariates <- trial_covariates
            n_added <- n_added + 1L
            model_name <- paste0("base_plus_", n_added)

            models[[length(models) + 1L]] <- list(
                name = model_name,
                formula = result$formula,
                design = result$design,
                covariates = current_covariates
            )

            message("Added '", cov, "' -> model '", model_name,
                    "' (", ncol(result$design), " columns, rank ",
                    qr(result$design)$rank, ")")
        }

        if (n_added == 0L) {
            warning("All candidate covariates caused rank deficiency — ",
                    "returning base model only")
        }
    }

    # -- Assemble skipped tibble -----------------------------------------------
    if (length(skipped) > 0L) {
        skipped_df <- tibble::as_tibble(do.call(rbind, skipped))
    } else {
        skipped_df <- tibble::tibble(
            variable = character(0), score = numeric(0),
            reason = character(0)
        )
    }

    # -- Return ----------------------------------------------------------------
    message("Done. ", length(models), " model(s) suggested.")
    list(
        covariate_ranking   = covariate_ranking,
        models              = models,
        skipped_covariates  = skipped_df,
        sample_weights_flag = sample_weights_flag,
        sample_weights_vars = sample_weights_vars
    )
}


# -- Internal helper: check whether a covariate set produces a valid design ----
#' @noRd
try_design <- function(base_formula, covariates, metadata, n_samples) {
    trial_formula <- extend_formula(base_formula, covariates)
    trial_design <- tryCatch(
        model.matrix(trial_formula, data = metadata),
        error = function(e) NULL
    )
    if (is.null(trial_design))
        return(list(valid = FALSE, reason = "model.matrix error"))
    if (nrow(trial_design) != n_samples)
        return(list(valid = FALSE, reason = paste0("drops ",
            n_samples - nrow(trial_design), " samples due to NAs")))
    design_rank <- qr(trial_design)$rank
    if (design_rank < ncol(trial_design))
        return(list(valid = FALSE, reason = paste0("rank deficient (rank ",
            design_rank, " < ", ncol(trial_design), " columns)")))
    list(valid = TRUE, design = trial_design, formula = trial_formula)
}


# -- Internal helper: extend a formula with additional covariates --------------
extend_formula <- function(base_formula, covariates) {
    if (length(covariates) == 0L) return(base_formula)
    base_str <- deparse(base_formula, width.cutoff = 500L)
    cov_str <- paste(covariates, collapse = " + ")
    as.formula(paste(base_str, "+", cov_str))
}
