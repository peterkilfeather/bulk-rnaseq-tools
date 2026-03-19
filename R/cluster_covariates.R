#' Cluster correlated metadata variables and select representatives
#'
#' Computes pairwise associations between metadata variables (Spearman
#' correlation, sqrt(eta-squared), or Cramer's V), clusters them
#' hierarchically, and selects one representative per cluster. Designed
#' to be called between \code{test_pc_associations()} and
#' \code{suggest_models()} to collapse redundant covariates.
#'
#' @param metadata   data.frame — same metadata passed to
#'                   \code{test_pc_associations()}.
#' @param pc_assoc   Output of \code{test_pc_associations()}, or NULL.
#'                   Used to score representatives by PC association strength.
#' @param pca        \code{prcomp} object or NULL. Used for PC variance
#'                   weighting when scoring. Uniform 1/k weights if NULL.
#' @param exclude    Character vector of column names to exclude (e.g.,
#'                   \code{"group"}, sample ID columns).
#' @param block      Character: blocking variable name. Always excluded from
#'                   clustering.
#' @param threshold  Numeric in [0, 1]: association threshold for clustering.
#'                   Variables with pairwise association >= threshold are
#'                   grouped. Default 0.6.
#' @param p_threshold Numeric: FDR threshold for "significant" PC associations
#'                   when scoring representatives. Default 0.05.
#' @param method     Character: \code{hclust} linkage method. Default
#'                   \code{"average"} (UPGMA).
#' @param plot       Logical: produce diagnostic plots. Default TRUE.
#' @param plot_file  Character or NULL: path to save the plot.
#'
#' @return A named list:
#'   \item{clusters}{Tibble with columns: variable, type, cluster,
#'         is_representative, score.}
#'   \item{representatives}{Character vector of representative variable names.}
#'   \item{assoc_matrix}{Numeric matrix of pairwise associations.}
#'   \item{hclust}{hclust object.}
#'   \item{threshold}{Threshold used.}
#'   \item{pc_assoc_filtered}{pc_assoc with results filtered to
#'         representatives only, or NULL.}
#'   \item{plot}{ggplot object, or NULL.}
cluster_covariates <- function(metadata, pc_assoc = NULL, pca = NULL,
                               exclude = NULL, block = NULL,
                               threshold = 0.6, p_threshold = 0.05,
                               method = "average",
                               plot = TRUE, plot_file = NULL) {

    # -- Input validation ------------------------------------------------------
    if (!is.data.frame(metadata)) {
        stop("`metadata` must be a data.frame or tibble")
    }
    if (threshold < 0 || threshold > 1) {
        stop("`threshold` must be in [0, 1]")
    }
    if (!is.null(pc_assoc)) {
        if (!is.list(pc_assoc) || is.null(pc_assoc$results)) {
            stop("`pc_assoc` must be the output of test_pc_associations() ",
                 "(a list with a $results tibble)")
        }
    }
    if (!is.null(pca) && !inherits(pca, "prcomp")) {
        stop("`pca` must be a prcomp object")
    }

    # -- Step 1: Select and classify variables ---------------------------------
    exclude_cols <- c(block, exclude)
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
        message("Dropping variables with < 2 unique non-NA values: ",
                paste(test_cols[!keep], collapse = ", "))
    }
    test_cols <- test_cols[keep]
    var_types <- var_types[keep]

    # Warn about high-cardinality categoricals
    n_samples <- nrow(metadata)
    for (col in test_cols[var_types == "categorical"]) {
        n_levels <- length(unique(stats::na.omit(metadata[[col]])))
        if (n_levels > n_samples / 2) {
            warning("'", col, "' has ", n_levels, " levels (> n/2 = ",
                    n_samples / 2, ") — likely a sample ID that should ",
                    "be excluded")
        }
    }

    # Early return if < 2 variables
    if (length(test_cols) < 2L) {
        message("< 2 variables remaining — no clustering performed.")
        clusters_tbl <- tibble::tibble(
            variable          = test_cols,
            type              = unname(var_types),
            cluster           = seq_along(test_cols),
            is_representative = TRUE,
            score             = NA_real_
        )
        pc_assoc_filt <- NULL
        if (!is.null(pc_assoc)) {
            pc_assoc_filt <- pc_assoc
            pc_assoc_filt$results <- pc_assoc$results[
                pc_assoc$results$variable %in% test_cols, , drop = FALSE
            ]
        }
        return(list(
            clusters          = clusters_tbl,
            representatives   = test_cols,
            assoc_matrix      = NULL,
            hclust            = NULL,
            threshold         = threshold,
            pc_assoc_filtered = pc_assoc_filt,
            plot              = NULL
        ))
    }

    # If pc_assoc provided, use intersection of variables
    if (!is.null(pc_assoc)) {
        pc_vars <- unique(pc_assoc$results$variable)
        in_both <- intersect(test_cols, pc_vars)
        only_meta <- setdiff(test_cols, pc_vars)
        only_pc <- setdiff(pc_vars, c(test_cols, exclude_cols))
        if (length(only_meta) > 0L) {
            message("Variables in metadata but not pc_assoc (will cluster ",
                    "but cannot score): ", paste(only_meta, collapse = ", "))
        }
        if (length(only_pc) > 0L) {
            message("Variables in pc_assoc but not metadata (ignored): ",
                    paste(only_pc, collapse = ", "))
        }
    }

    # -- Step 2: Compute pairwise association matrix ---------------------------
    assoc_mat <- build_assoc_matrix(metadata, test_cols, var_types)

    # -- Step 3: Hierarchical clustering ---------------------------------------
    dist_mat <- stats::as.dist(1 - assoc_mat)
    hc <- stats::hclust(dist_mat, method = method)
    clusters <- stats::cutree(hc, h = 1 - threshold)

    n_clusters <- max(clusters)
    message("Analysed ", length(test_cols), " variables -> ",
            n_clusters, " cluster(s) at threshold ", threshold)

    # Warn if all in one cluster
    if (n_clusters == 1L && length(test_cols) > 1L) {
        warning("All variables in a single cluster — threshold may be too low")
    }

    # -- Step 4: Select representatives ----------------------------------------
    scores <- rep(NA_real_, length(test_cols))
    names(scores) <- test_cols

    if (!is.null(pc_assoc)) {
        scores <- score_variables(pc_assoc, pca, test_cols, p_threshold)
    }

    representatives <- character(n_clusters)
    for (k in seq_len(n_clusters)) {
        members <- test_cols[clusters == k]

        if (length(members) == 1L) {
            representatives[k] <- members
            next
        }

        if (!is.null(pc_assoc)) {
            member_scores <- scores[members]
            if (any(member_scores > 0, na.rm = TRUE)) {
                best <- members[which.max(member_scores)]
                representatives[k] <- best
                message("  Cluster ", k, " (", paste(members, collapse = ", "),
                        "): representative = '", best,
                        "' (highest PC association score: ",
                        round(scores[best], 4), ")")
                next
            }
        }

        # Fallback: most central variable
        sub_mat <- assoc_mat[members, members, drop = FALSE]
        centrality <- rowMeans(sub_mat) # includes self (=1), but consistent
        best <- members[which.max(centrality)]
        representatives[k] <- best
        message("  Cluster ", k, " (", paste(members, collapse = ", "),
                "): representative = '", best, "' (highest centrality: ",
                round(centrality[best], 3), ")")
    }

    # -- Build output tibble ---------------------------------------------------
    clusters_tbl <- tibble::tibble(
        variable          = test_cols,
        type              = unname(var_types),
        cluster           = unname(clusters),
        is_representative = test_cols %in% representatives,
        score             = unname(scores[test_cols])
    )

    # -- Step 5: Build pc_assoc_filtered ---------------------------------------
    pc_assoc_filt <- NULL
    if (!is.null(pc_assoc)) {
        pc_assoc_filt <- pc_assoc
        pc_assoc_filt$results <- pc_assoc$results[
            pc_assoc$results$variable %in% representatives, , drop = FALSE
        ]
    }

    # -- Step 6: Plot ----------------------------------------------------------
    p <- NULL
    if (plot) {
        p <- build_cluster_plot(assoc_mat, hc, clusters, representatives,
                                threshold)

        if (!is.null(plot_file)) {
            n_vars <- length(test_cols)
            w <- 4 + n_vars * 0.5
            h <- 2 + n_vars * 0.4
            ggplot2::ggsave(plot_file, plot = p, width = w, height = h)
            message("Saved plot to ", plot_file)
        }
    }

    # -- Return ----------------------------------------------------------------
    list(
        clusters          = clusters_tbl,
        representatives   = representatives,
        assoc_matrix      = assoc_mat,
        hclust            = hc,
        threshold         = threshold,
        pc_assoc_filtered = pc_assoc_filt,
        plot              = p
    )
}


# -- Internal helpers ----------------------------------------------------------

#' Compute pairwise association between two variables
#' @noRd
compute_association <- function(x, y, type_x, type_y) {
    # Pairwise complete observations
    complete <- !is.na(x) & !is.na(y)
    if (sum(complete) < 3L) return(0)
    x <- x[complete]
    y <- y[complete]

    # Both numeric: absolute Spearman correlation
    if (type_x == "continuous" && type_y == "continuous") {
        if (stats::sd(x) == 0 || stats::sd(y) == 0) return(0)
        r <- stats::cor(x, y, method = "spearman", use = "pairwise.complete.obs")
        result <- abs(r)
        return(if (is.nan(result)) 0 else result)
    }

    # Numeric x Categorical: sqrt(eta-squared)
    if (type_x == "continuous" && type_y == "categorical") {
        return(compute_eta(x, y))
    }
    if (type_x == "categorical" && type_y == "continuous") {
        return(compute_eta(y, x))
    }

    # Both categorical: Cramer's V
    compute_cramers_v(x, y)
}

#' Compute sqrt(eta-squared) from one-way ANOVA
#' @noRd
compute_eta <- function(num, cat) {
    cat <- factor(cat)
    if (nlevels(cat) < 2L) return(0)
    if (stats::sd(num) == 0) return(0)
    fit <- stats::lm(num ~ cat)
    a <- stats::anova(fit)
    ss_between <- a[1L, "Sum Sq"]
    ss_total <- sum(a[["Sum Sq"]])
    if (ss_total == 0) return(0)
    result <- sqrt(ss_between / ss_total)
    if (is.nan(result)) 0 else result
}

#' Compute Cramer's V
#' @noRd
compute_cramers_v <- function(x, y) {
    x <- factor(x)
    y <- factor(y)
    if (nlevels(x) < 2L || nlevels(y) < 2L) return(0)
    tbl <- table(x, y)
    chi2 <- suppressWarnings(stats::chisq.test(tbl)$statistic)
    n <- sum(tbl)
    k <- min(nrow(tbl), ncol(tbl))
    if (k <= 1L || n == 0L) return(0)
    result <- sqrt(as.numeric(chi2) / (n * (k - 1)))
    if (is.nan(result) || is.na(result)) 0 else min(result, 1)
}

#' Build pairwise association matrix
#' @noRd
build_assoc_matrix <- function(metadata, var_names, var_types) {
    n <- length(var_names)
    mat <- matrix(0, nrow = n, ncol = n,
                  dimnames = list(var_names, var_names))
    diag(mat) <- 1

    for (i in seq_len(n - 1L)) {
        for (j in (i + 1L):n) {
            val <- compute_association(
                metadata[[var_names[i]]], metadata[[var_names[j]]],
                var_types[i], var_types[j]
            )
            mat[i, j] <- val
            mat[j, i] <- val
        }
    }
    mat
}

#' Score variables by PC association strength
#'
#' Replicates the scoring logic from suggest_models: for each variable,
#' score = sum(r_squared * variance_weight) across significant PCs.
#' @noRd
score_variables <- function(pc_assoc, pca, var_names, p_threshold) {
    results <- pc_assoc$results
    pc_responses <- grep("^PC\\d+", unique(results$response), value = TRUE)
    n_pcs <- length(pc_responses)

    if (n_pcs == 0L) {
        return(stats::setNames(rep(0, length(var_names)), var_names))
    }

    # Variance weights
    if (!is.null(pca)) {
        all_var <- pca$sdev^2 / sum(pca$sdev^2)
        pc_indices <- as.integer(sub("^PC", "", pc_responses))
        var_weights <- all_var[pc_indices]
    } else {
        var_weights <- rep(1 / n_pcs, n_pcs)
    }
    var_weights <- var_weights / sum(var_weights)
    names(var_weights) <- pc_responses

    # Significant PC results only
    pc_results <- results[grepl("^PC\\d+", results$response), , drop = FALSE]
    sig_pc <- pc_results[!is.na(pc_results$p_adj) &
                         pc_results$p_adj < p_threshold, , drop = FALSE]

    scores <- vapply(var_names, function(v) {
        v_sig <- sig_pc[sig_pc$variable == v, , drop = FALSE]
        if (nrow(v_sig) == 0L) return(0)
        sum(v_sig$r_squared * var_weights[v_sig$response], na.rm = TRUE)
    }, numeric(1))

    scores
}

#' Extract dendrogram segment data from hclust for ggplot2
#' @noRd
dendro_data <- function(hc) {
    n <- length(hc$labels)
    # Leaf x positions according to hc$order
    leaf_x <- numeric(n)
    leaf_x[hc$order] <- seq_len(n)

    segments <- list()
    # node_pos tracks the x position of each merged node
    node_pos <- numeric(n - 1L)

    for (i in seq_len(n - 1L)) {
        left  <- hc$merge[i, 1]
        right <- hc$merge[i, 2]
        h <- hc$height[i]

        # Get x positions and heights of children
        if (left < 0) {
            x_left <- leaf_x[-left]
            h_left <- 0
        } else {
            x_left <- node_pos[left]
            h_left <- hc$height[left]
        }
        if (right < 0) {
            x_right <- leaf_x[-right]
            h_right <- 0
        } else {
            x_right <- node_pos[right]
            h_right <- hc$height[right]
        }

        node_pos[i] <- (x_left + x_right) / 2

        # Vertical segments from children up to merge height
        segments[[length(segments) + 1L]] <- data.frame(
            x = x_left, y = h_left, xend = x_left, yend = h)
        segments[[length(segments) + 1L]] <- data.frame(
            x = x_right, y = h_right, xend = x_right, yend = h)
        # Horizontal segment connecting children
        segments[[length(segments) + 1L]] <- data.frame(
            x = x_left, y = h, xend = x_right, yend = h)
    }

    seg_df <- do.call(rbind, segments)

    # Leaf data
    leaf_df <- data.frame(
        x     = seq_len(n),
        label = hc$labels[hc$order],
        stringsAsFactors = FALSE
    )

    list(segments = seg_df, leaves = leaf_df)
}

#' Build two-panel cluster diagnostic plot
#' @noRd
build_cluster_plot <- function(assoc_mat, hc, clusters, representatives,
                               threshold) {
    var_names <- colnames(assoc_mat)
    n_vars <- length(var_names)

    # -- Panel A: Dendrogram ---------------------------------------------------
    dd <- dendro_data(hc)

    # Colour leaves by cluster
    ordered_labels <- hc$labels[hc$order]
    leaf_clusters <- clusters[ordered_labels]
    dd$leaves$cluster <- factor(leaf_clusters)
    dd$leaves$is_rep <- ordered_labels %in% representatives

    p_dendro <- ggplot2::ggplot() +
        ggplot2::geom_segment(
            data = dd$segments,
            ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
            color = "grey40"
        ) +
        ggplot2::geom_hline(yintercept = 1 - threshold,
                            linetype = "dashed", color = "red") +
        ggplot2::geom_text(
            data = dd$leaves,
            ggplot2::aes(x = x, y = -0.02, label = label, color = cluster),
            angle = 90, hjust = 1, size = 3, show.legend = FALSE
        ) +
        ggplot2::scale_y_continuous(
            expand = ggplot2::expansion(mult = c(0.3, 0.05))
        ) +
        ggplot2::labs(y = "Distance (1 - association)", x = NULL) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.x  = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            panel.grid.major.x = ggplot2::element_blank()
        )

    # -- Panel B: Heatmap ------------------------------------------------------
    ordered_vars <- hc$labels[hc$order]

    # Build long-form data
    heat_df <- expand.grid(
        var1 = ordered_vars,
        var2 = ordered_vars,
        stringsAsFactors = FALSE
    )
    heat_df$value <- vapply(seq_len(nrow(heat_df)), function(i) {
        assoc_mat[heat_df$var1[i], heat_df$var2[i]]
    }, numeric(1))
    heat_df$var1 <- factor(heat_df$var1, levels = ordered_vars)
    heat_df$var2 <- factor(heat_df$var2, levels = rev(ordered_vars))

    # Mark representatives in axis labels
    labels_x <- ordered_vars
    labels_y <- rev(ordered_vars)
    labels_x[labels_x %in% representatives] <- paste0(
        labels_x[labels_x %in% representatives], " *")
    labels_y[labels_y %in% representatives] <- paste0(
        labels_y[labels_y %in% representatives], " *")

    p_heat <- ggplot2::ggplot(
        heat_df,
        ggplot2::aes(x = var1, y = var2, fill = value)
    ) +
        ggplot2::geom_tile(color = "white")

    # Add cell values for small matrices
    if (n_vars <= 15) {
        p_heat <- p_heat +
            ggplot2::geom_text(
                ggplot2::aes(label = sprintf("%.2f", value)),
                size = max(2, 4 - n_vars * 0.15)
            )
    }

    # Cluster boundary rectangles
    cluster_ranges <- list()
    for (k in unique(clusters[ordered_vars])) {
        idx <- which(clusters[ordered_vars] == k)
        cluster_ranges[[length(cluster_ranges) + 1L]] <- data.frame(
            xmin = min(idx) - 0.5,
            xmax = max(idx) + 0.5,
            ymin = n_vars - max(idx) + 0.5,
            ymax = n_vars - min(idx) + 1.5
        )
    }
    rect_df <- do.call(rbind, cluster_ranges)

    p_heat <- p_heat +
        ggplot2::geom_rect(
            data = rect_df,
            ggplot2::aes(xmin = xmin, xmax = xmax,
                         ymin = ymin, ymax = ymax),
            fill = NA, color = "red", linewidth = 0.8,
            inherit.aes = FALSE
        ) +
        ggplot2::scale_fill_gradient(
            low = "white", high = "darkblue",
            name = "Association", limits = c(0, 1)
        ) +
        ggplot2::scale_x_discrete(labels = labels_x) +
        ggplot2::scale_y_discrete(labels = labels_y) +
        ggplot2::labs(x = NULL, y = NULL) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
        )

    # Combine with patchwork
    patchwork::wrap_plots(p_dendro, p_heat, ncol = 2, widths = c(1, 1.5))
}
