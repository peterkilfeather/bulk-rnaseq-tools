#' PCA score plots colored by metadata variables
#'
#' Produces a multi-panel scatter plot of PCA scores, one panel per metadata
#' variable tested in \code{test_pc_associations()}. Each panel colours points
#' by that variable's values. Serves as a visual companion to the association
#' heatmap from \code{test_pc_associations()}.
#'
#' @param pca       A \code{prcomp} object.
#' @param metadata  A data.frame with one row per sample (same as passed to
#'                  \code{test_pc_associations()}).
#' @param pc_assoc  Output of \code{test_pc_associations()} (a list with a
#'                  \code{$results} tibble).
#' @param pcs       Integer vector of length 2: default PC axes (default
#'                  \code{c(1L, 2L)}).
#' @param per_variable_pcs Logical; if TRUE, each panel plots the top two PCs
#'                  most significantly associated with that variable instead of
#'                  using the fixed \code{pcs} pair. Default FALSE.
#' @param p_threshold Adjusted p-value threshold for significance stars and
#'                  for selecting per-variable PCs. Default 0.05.
#' @param point_size Numeric; point size passed to \code{geom_point()}.
#'                  Default 2.
#' @param plot_file File path to save the combined plot. When NULL the plot
#'                  is not saved to file.
#'
#' @return A named list:
#'   \item{combined_plot}{patchwork object (all panels arranged).}
#'   \item{individual_plots}{Named list of ggplot objects, one per variable.}
#'   \item{pc_data}{Named list of data.frames, one per variable, each with
#'         columns \code{pc_x}, \code{pc_y}, and the colouring variable.}
plot_pca <- function(pca, metadata, pc_assoc,
                     pcs = c(1L, 2L),
                     per_variable_pcs = FALSE,
                     p_threshold = 0.05,
                     point_size = 2,
                     plot_file = NULL) {

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
    if (!is.list(pc_assoc) || is.null(pc_assoc$results)) {
        stop("`pc_assoc` must be the output of test_pc_associations() ",
             "(a list with a $results tibble)")
    }
    pcs <- as.integer(pcs)
    if (length(pcs) != 2L || anyNA(pcs)) {
        stop("`pcs` must be an integer vector of length 2")
    }
    if (any(pcs < 1L) || any(pcs > ncol(pca$x))) {
        stop("`pcs` values must be between 1 and ", ncol(pca$x))
    }
    if (!is.logical(per_variable_pcs) || length(per_variable_pcs) != 1L) {
        stop("`per_variable_pcs` must be TRUE or FALSE")
    }

    # -- Extract variable list and types from pc_assoc -------------------------
    results <- pc_assoc$results
    pc_results <- results[grepl("^PC\\d+$", results$response), ]
    var_info <- unique(pc_results[, c("variable", "type")])
    test_vars <- var_info$variable
    var_types <- stats::setNames(var_info$type, var_info$variable)

    if (length(test_vars) == 0L) {
        stop("No PC-associated variables found in `pc_assoc$results`")
    }

    missing_vars <- setdiff(test_vars, colnames(metadata))
    if (length(missing_vars) > 0L) {
        stop("Variables in `pc_assoc$results` not found in `metadata`: ",
             paste(missing_vars, collapse = ", "))
    }

    # -- Variance explained ----------------------------------------------------
    var_explained <- pca$sdev^2 / sum(pca$sdev^2) * 100

    # -- Helper: pick top 2 PCs for a variable ---------------------------------
    pick_pcs <- function(v) {
        v_rows <- pc_results[pc_results$variable == v, ]
        v_rows <- v_rows[order(v_rows$p_adj), ]
        sig_rows <- v_rows[!is.na(v_rows$p_adj) &
                           v_rows$p_adj < p_threshold, ]
        if (nrow(sig_rows) >= 2L) {
            top2 <- as.integer(sub("^PC", "", sig_rows$response[1:2]))
        } else if (nrow(sig_rows) == 1L) {
            pc1 <- as.integer(sub("^PC", "", sig_rows$response[1]))
            pc2 <- setdiff(pcs, pc1)[1]
            if (is.na(pc2)) pc2 <- pcs[2]
            top2 <- c(pc1, pc2)
        } else {
            top2 <- pcs
        }
        top2
    }

    # -- Helper: significance stars --------------------------------------------
    format_stars <- function(p) {
        if (is.na(p))       ""
        else if (p < 0.001) "***"
        else if (p < 0.01)  "**"
        else if (p < 0.05)  "*"
        else                ""
    }

    # -- Build one plot per variable -------------------------------------------
    individual_plots <- list()
    pc_data <- list()

    for (v in test_vars) {
        # Determine PC pair
        if (per_variable_pcs) {
            var_pcs <- pick_pcs(v)
        } else {
            var_pcs <- pcs
        }

        # Axis labels
        xlab <- sprintf("PC%d (%.1f%%)", var_pcs[1], var_explained[var_pcs[1]])
        ylab <- sprintf("PC%d (%.1f%%)", var_pcs[2], var_explained[var_pcs[2]])

        # Build data frame
        df <- data.frame(
            pc_x = pca$x[, var_pcs[1]],
            pc_y = pca$x[, var_pcs[2]],
            color_var = metadata[[v]],
            stringsAsFactors = FALSE
        )

        # Drop rows where color_var is NA
        complete <- !is.na(df$color_var)
        if (sum(complete) == 0L) {
            message("Skipping '", v, "': all values are NA")
            next
        }
        if (any(!complete)) {
            message("'", v, "': dropping ", sum(!complete), " NA samples")
        }
        df <- df[complete, , drop = FALSE]

        # Store pc_data
        pc_data[[v]] <- df

        # Significance stars for this variable against plotted PCs
        v_pc_rows <- pc_results[pc_results$variable == v &
                                pc_results$response %in%
                                    paste0("PC", var_pcs), ]
        min_p <- if (nrow(v_pc_rows) > 0L) {
            min(v_pc_rows$p_adj, na.rm = TRUE)
        } else {
            NA_real_
        }
        stars <- format_stars(min_p)
        title <- if (nzchar(stars)) paste(v, stars) else v

        # Build ggplot
        p <- ggplot2::ggplot(df, ggplot2::aes(x = pc_x, y = pc_y,
                                               color = color_var)) +
            ggplot2::geom_point(size = point_size) +
            ggplot2::labs(x = xlab, y = ylab, title = title, color = v) +
            ggplot2::theme_minimal()

        # Color scale depends on variable type
        if (var_types[v] == "continuous") {
            p <- p + ggplot2::scale_color_gradient(low = "#4575B4",
                                                   high = "#D73027")
        } else {
            df$color_var <- factor(df$color_var)
            n_levels <- nlevels(df$color_var)

            # Rebuild plot with factor for discrete scale
            p <- ggplot2::ggplot(df, ggplot2::aes(x = pc_x, y = pc_y,
                                                   color = color_var)) +
                ggplot2::geom_point(size = point_size) +
                ggplot2::labs(x = xlab, y = ylab, title = title, color = v) +
                ggplot2::theme_minimal()

            if (n_levels <= 8L) {
                p <- p + ggplot2::scale_color_brewer(palette = "Set2")
            }
        }

        individual_plots[[v]] <- p
    }

    if (length(individual_plots) == 0L) {
        stop("No plottable variables (all had NA-only values)")
    }

    # -- Compose multi-panel figure --------------------------------------------
    n_vars <- length(individual_plots)
    ncol_layout <- ceiling(sqrt(n_vars))
    combined <- patchwork::wrap_plots(individual_plots, ncol = ncol_layout)

    # -- Save if requested -----------------------------------------------------
    if (!is.null(plot_file)) {
        nrow_layout <- ceiling(n_vars / ncol_layout)
        w <- ncol_layout * 4
        h <- nrow_layout * 3.5
        ggplot2::ggsave(plot_file, plot = combined, width = w, height = h)
        message("Saved PCA plots to ", plot_file)
    }

    message("Done. ", n_vars, " panel(s) plotted.")
    list(
        combined_plot    = combined,
        individual_plots = individual_plots,
        pc_data          = pc_data
    )
}
