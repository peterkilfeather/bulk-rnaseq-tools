#' Test associations between metadata variables and principal components
#'
#' For each metadata variable, fits a linear model (or mixed model when
#' \code{block} is provided) against each of the top PCs and, optionally,
#' sample quality weights.
#'
#' @param pca        A \code{prcomp} object (pre-computed PCA).
#' @param metadata   A data.frame or tibble with one row per sample. Numeric
#'                   columns are tested as continuous (t-test); character, factor
#'                   and logical columns as categorical (F-test).
#' @param n_pcs      Number of top PCs to test (default 10, capped at available).
#' @param block      Column name in \code{metadata} identifying the blocking /
#'                   subject variable. When provided, \code{lmerTest::lmer()}
#'                   with a random intercept \code{(1|block)} is used.
#' @param adjust     Character vector of column names in \code{metadata} to
#'                   include as covariates. These columns are NOT tested.
#' @param sample_weights Numeric vector of sample quality weights (e.g., from
#'                   \code{normalise_counts()$sample_weights}). Tested as an
#'                   additional response alongside PCs.
#' @param plot       Logical; produce a heatmap of -log10(p.adj) (default TRUE).
#' @param plot_file  File path to save the heatmap. When NULL the plot is not
#'                   saved to file.
#'
#' @return A named list:
#'   \item{results}{Tidy tibble with one row per variable x response pair.}
#'   \item{plot}{ggplot object, or NULL if \code{plot = FALSE}.}
test_pc_associations <- function(pca, metadata, n_pcs = 10, block = NULL,
                                  adjust = NULL, sample_weights = NULL,
                                  plot = TRUE, plot_file = NULL) {

    # -- Input validation ------------------------------------------------------
    if (!inherits(pca, "prcomp")) {
        stop("`pca` must be a prcomp object")
    }
    if (!is.data.frame(metadata)) {
        stop("`metadata` must be a data.frame or tibble")
    }
    n_samples <- nrow(pca$x)
    if (nrow(metadata) != n_samples) {
        stop("`metadata` must have the same number of rows as samples in `pca` (",
             n_samples, " expected, got ", nrow(metadata), ")")
    }
    n_pcs <- min(as.integer(n_pcs), ncol(pca$x))
    if (!is.null(block)) {
        if (length(block) != 1L || !block %in% colnames(metadata)) {
            stop("`block` must be a single column name present in `metadata`")
        }
    }
    if (!is.null(adjust)) {
        missing_adj <- setdiff(adjust, colnames(metadata))
        if (length(missing_adj) > 0L) {
            stop("adjust columns not found in metadata: ",
                 paste(missing_adj, collapse = ", "))
        }
    }
    if (!is.null(sample_weights)) {
        if (length(sample_weights) != n_samples) {
            stop("`sample_weights` must have length equal to number of samples (",
                 n_samples, ")")
        }
    }

    # -- Variance explained per PC ---------------------------------------------
    var_explained <- pca$sdev^2 / sum(pca$sdev^2) * 100
    pc_labels <- sprintf("PC%d (%.1f%%)", seq_len(n_pcs), var_explained[seq_len(n_pcs)])
    pc_names  <- paste0("PC", seq_len(n_pcs))

    # -- Build response matrix -------------------------------------------------
    responses <- as.data.frame(pca$x[, seq_len(n_pcs), drop = FALSE])
    colnames(responses) <- pc_names
    response_labels <- stats::setNames(pc_labels, pc_names)

    if (!is.null(sample_weights)) {
        responses$sample_weights <- sample_weights
        response_labels <- c(response_labels, sample_weights = "Sample weights")
    }
    response_names <- colnames(responses)

    # -- Determine test variables and types ------------------------------------
    exclude_cols <- c(block, adjust)
    test_cols <- setdiff(colnames(metadata), exclude_cols)

    var_types <- vapply(test_cols, function(col) {
        x <- metadata[[col]]
        if (is.numeric(x)) "continuous" else "categorical"
    }, character(1))

    # Drop variables with < 2 unique non-NA values
    keep <- vapply(test_cols, function(col) {
        length(unique(stats::na.omit(metadata[[col]]))) >= 2L
    }, logical(1))
    if (any(!keep)) {
        warning("Skipping variables with < 2 unique non-NA values: ",
                paste(test_cols[!keep], collapse = ", "))
    }
    test_cols <- test_cols[keep]
    var_types <- var_types[keep]

    if (length(test_cols) == 0L) {
        stop("No testable variables remaining in `metadata`")
    }

    # -- Helper: marginal R-squared for lmer -----------------------------------
    marginal_r2 <- function(fit) {
        var_fixed  <- stats::var(stats::predict(fit, re.form = NA))
        var_random <- sum(vapply(lme4::VarCorr(fit),
                                 function(v) v[1], numeric(1)))
        var_resid  <- stats::sigma(fit)^2
        var_fixed / (var_fixed + var_random + var_resid)
    }

    # -- Helper: fit one (variable, response) pair -----------------------------
    fit_one <- function(var_name, var_type, resp_name) {
        df <- data.frame(y = responses[[resp_name]])
        test_var <- metadata[[var_name]]
        if (var_type == "categorical") {
            test_var <- factor(test_var)
        }
        df$test_var <- test_var

        # Add adjustment covariates
        if (!is.null(adjust)) {
            for (a in adjust) {
                val <- metadata[[a]]
                if (is.character(val) || is.logical(val)) val <- factor(val)
                df[[a]] <- val
            }
        }

        # Add block
        if (!is.null(block)) {
            df$.block <- metadata[[block]]
        }

        # Drop rows with NA in test variable or response
        complete <- stats::complete.cases(df)
        if (sum(complete) < 3L) {
            return(data.frame(
                variable  = var_name,
                response  = resp_name,
                type      = var_type,
                estimate  = NA_real_,
                eta_sq    = NA_real_,
                r_squared = NA_real_,
                statistic = NA_real_,
                df1       = NA_real_,
                df2       = NA_real_,
                p_value   = NA_real_,
                stringsAsFactors = FALSE
            ))
        }
        df <- df[complete, , drop = FALSE]

        # Check sufficient levels after NA removal
        if (var_type == "categorical" && length(unique(df$test_var)) < 2L) {
            return(data.frame(
                variable  = var_name,
                response  = resp_name,
                type      = var_type,
                estimate  = NA_real_,
                eta_sq    = NA_real_,
                r_squared = NA_real_,
                statistic = NA_real_,
                df1       = NA_real_,
                df2       = NA_real_,
                p_value   = NA_real_,
                stringsAsFactors = FALSE
            ))
        }

        # Build formula strings
        adj_terms <- if (!is.null(adjust)) paste(adjust, collapse = " + ") else NULL
        full_rhs <- paste(
            c("test_var", adj_terms),
            collapse = " + "
        )
        reduced_rhs <- if (!is.null(adj_terms)) adj_terms else "1"

        na_row <- data.frame(
            variable  = var_name,
            response  = resp_name,
            type      = var_type,
            estimate  = NA_real_,
            eta_sq    = NA_real_,
            r_squared = NA_real_,
            statistic = NA_real_,
            df1       = NA_real_,
            df2       = NA_real_,
            p_value   = NA_real_,
            stringsAsFactors = FALSE
        )

        singular_hit <- FALSE
        result <- tryCatch(
            withCallingHandlers({
                if (!is.null(block)) {
                    # -- Mixed model -------------------------------------------
                    full_formula    <- stats::as.formula(
                        paste("y ~", full_rhs, "+ (1 | .block)"))
                    reduced_formula <- stats::as.formula(
                        paste("y ~", reduced_rhs, "+ (1 | .block)"))

                    full_fit    <- lmerTest::lmer(full_formula, data = df)
                    reduced_fit <- lme4::lmer(reduced_formula, data = df)

                    r2_inc <- marginal_r2(full_fit) - marginal_r2(reduced_fit)

                    if (var_type == "continuous") {
                        s    <- summary(full_fit)$coefficients
                        est  <- s["test_var", "Estimate"]
                        stat <- s["test_var", "t value"]
                        pval <- s["test_var", "Pr(>|t|)"]
                        d1   <- 1
                        d2   <- s["test_var", "df"]
                        eta  <- NA_real_
                    } else {
                        a    <- anova(full_fit)
                        stat <- a["test_var", "F value"]
                        pval <- a["test_var", "Pr(>F)"]
                        d1   <- a["test_var", "NumDF"]
                        d2   <- a["test_var", "DenDF"]
                        est  <- NA_real_
                        eta  <- (stat * d1) / (stat * d1 + d2)
                    }
                } else {
                    # -- Standard linear model ---------------------------------
                    full_formula    <- stats::as.formula(paste("y ~", full_rhs))
                    reduced_formula <- stats::as.formula(paste("y ~", reduced_rhs))

                    full_fit    <- stats::lm(full_formula, data = df)
                    reduced_fit <- stats::lm(reduced_formula, data = df)

                    r2_inc <- summary(full_fit)$r.squared -
                              summary(reduced_fit)$r.squared

                    if (var_type == "continuous") {
                        s    <- summary(full_fit)$coefficients
                        est  <- s["test_var", "Estimate"]
                        stat <- s["test_var", "t value"]
                        pval <- s["test_var", "Pr(>|t|)"]
                        d1   <- 1
                        d2   <- as.numeric(full_fit$df.residual)
                        eta  <- NA_real_
                    } else {
                        a    <- stats::anova(full_fit)
                        stat <- a["test_var", "F value"]
                        pval <- a["test_var", "Pr(>F)"]
                        d1   <- as.numeric(a["test_var", "Df"])
                        d2   <- as.numeric(a["Residuals", "Df"])
                        est  <- NA_real_
                        eta  <- a["test_var", "Sum Sq"] /
                                sum(a[["Sum Sq"]])
                    }
                }

                data.frame(
                    variable  = var_name,
                    response  = resp_name,
                    type      = var_type,
                    estimate  = est,
                    eta_sq    = eta,
                    r_squared = r2_inc,
                    statistic = stat,
                    df1       = d1,
                    df2       = d2,
                    p_value   = pval,
                    stringsAsFactors = FALSE
                )
            }, warning = function(w) {
                if (grepl("singular|converge", w$message, ignore.case = TRUE)) {
                    warning("Model warning for ", var_name, " vs ", resp_name,
                            ": ", w$message, call. = FALSE)
                    singular_hit <<- TRUE
                }
                invokeRestart("muffleWarning")
            }),
            error = function(e) {
                warning("Model failed for ", var_name, " vs ", resp_name,
                        ": ", e$message, call. = FALSE)
                na_row
            }
        )

        if (singular_hit) result <- na_row

        result
    }

    # -- Fit all models --------------------------------------------------------
    n_vars <- length(test_cols)
    n_resp <- length(response_names)
    message("Testing associations for ", n_vars, " variables across ",
            n_resp, " responses...")

    rows <- vector("list", n_vars * n_resp)
    idx  <- 0L
    for (i in seq_along(test_cols)) {
        message("Testing variable ", i, "/", n_vars, ": ", test_cols[i], "...")
        for (j in seq_along(response_names)) {
            idx <- idx + 1L
            rows[[idx]] <- fit_one(test_cols[i], var_types[i], response_names[j])
        }
    }

    # -- Assemble results ------------------------------------------------------
    message("Applying BH correction...")
    results <- do.call(rbind, rows)
    results$p_adj <- stats::p.adjust(results$p_value, method = "BH")
    results <- tibble::as_tibble(results)

    # -- Heatmap ---------------------------------------------------------------
    p <- NULL
    if (plot) {
        plot_df <- results
        plot_df$neg_log10_p <- pmin(-log10(plot_df$p_adj), 10)
        plot_df$neg_log10_p[is.na(plot_df$neg_log10_p)] <- 0

        # Significance stars
        plot_df$stars <- ifelse(is.na(plot_df$p_adj), "",
                         ifelse(plot_df$p_adj < 0.001, "***",
                         ifelse(plot_df$p_adj < 0.01,  "**",
                         ifelse(plot_df$p_adj < 0.05,  "*", ""))))

        # Pretty response labels
        all_labels <- response_labels
        plot_df$response_label <- all_labels[plot_df$response]
        plot_df$response_label <- factor(plot_df$response_label,
                                         levels = all_labels)

        p <- ggplot2::ggplot(
            plot_df,
            ggplot2::aes(x = response_label, y = variable,
                         fill = neg_log10_p)
        ) +
            ggplot2::geom_tile(color = "grey80") +
            ggplot2::geom_text(ggplot2::aes(label = stars), size = 4) +
            ggplot2::scale_fill_gradient(
                low = "white", high = "red",
                name = expression(-log[10](p[adj]))
            ) +
            ggplot2::theme_minimal() +
            ggplot2::labs(x = NULL, y = NULL) +
            ggplot2::theme(
                axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
            )

        if (!is.null(plot_file)) {
            w <- 2 + n_resp * 0.8
            h <- 1.5 + n_vars * 0.5
            ggplot2::ggsave(plot_file, plot = p, width = w, height = h)
            message("Saved heatmap to ", plot_file)
        }
    }

    message("Done.")
    list(results = results, plot = p)
}
