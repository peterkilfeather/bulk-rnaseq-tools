#' Build a model specification from a user-supplied formula
#'
#' Wraps a formula into the same output structure as \code{suggest_models()} so
#' it can be passed directly to \code{compare_voom_models()} and
#' \code{select_model()}.
#'
#' @param formula   A formula object (e.g., \code{~ 0 + group + age}).
#' @param metadata  data.frame used to build the design matrix via
#'                  \code{model.matrix()}.
#' @param block     Character; blocking variable column name (required because
#'                  the experiment has repeated measures).
#' @param name      Character; label for the model (default \code{"custom"}).
#' @param drop_redundant Logical; if \code{TRUE}, automatically drop
#'   non-estimable (collinear) columns from a rank-deficient design matrix with
#'   a warning.  If \code{FALSE} (default), stop with an error that names the
#'   redundant column(s).
#'
#' @return A named list matching the \code{suggest_models()} output structure:
#'   \item{covariate_ranking}{NULL (no PCA ranking performed).}
#'   \item{models}{List of length 1 containing a model specification with name,
#'         formula, design matrix, and covariates.}
#'   \item{skipped_covariates}{Empty tibble.}
#'   \item{sample_weights_flag}{FALSE (no PCA data to assess).}
#'   \item{sample_weights_vars}{Empty character vector.}
#'
#' @export
make_model <- function(formula, metadata, block, name = "custom",
                       drop_redundant = FALSE) {

    # -- Input validation ------------------------------------------------------
    if (!inherits(formula, "formula")) {
        stop("`formula` must be a formula (e.g., ~ 0 + group)")
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

    if (!is.character(name) || length(name) != 1L) {
        stop("`name` must be a single character string")
    }

    if (!is.logical(drop_redundant) || length(drop_redundant) != 1L) {
        stop("`drop_redundant` must be TRUE or FALSE")
    }

    # Check formula variables exist in metadata
    formula_vars <- all.vars(formula)
    missing_vars <- setdiff(formula_vars, colnames(metadata))
    if (length(missing_vars) > 0L) {
        stop("Formula variables not found in metadata: ",
             paste(missing_vars, collapse = ", "))
    }

    # -- Build design matrix ---------------------------------------------------
    n_samples <- nrow(metadata)

    design <- tryCatch(
        model.matrix(formula, data = metadata),
        error = function(e) {
            stop("model.matrix() failed: ", conditionMessage(e))
        }
    )

    if (nrow(design) != n_samples) {
        stop("Design matrix drops ", n_samples - nrow(design),
             " samples due to NAs in formula variables")
    }

    design_rank <- qr(design)$rank
    if (design_rank < ncol(design)) {
        redundant <- limma::nonEstimable(design)
        if (drop_redundant) {
            warning("Design matrix is rank deficient (rank ", design_rank,
                    " < ", ncol(design), " columns). ",
                    "Dropping non-estimable column(s): ",
                    paste(redundant, collapse = ", "))
            design <- design[, !colnames(design) %in% redundant, drop = FALSE]
        } else {
            stop("Design matrix is rank deficient (rank ", design_rank,
                 " < ", ncol(design), " columns). ",
                 "Non-estimable coefficient(s): ",
                 paste(redundant, collapse = ", "),
                 ". Consider removing the corresponding term(s) from your formula, ",
                 "or set drop_redundant = TRUE to drop them automatically.")
        }
    }

    # -- Identify covariates (formula terms beyond base) -----------------------
    covariates <- formula_vars

    # -- Assemble output -------------------------------------------------------
    message("Note: sample_weights_flag is set to FALSE because no PCA ",
            "association data was provided. Run suggest_models() if you need ",
            "sample-weight guidance.")

    list(
        covariate_ranking   = NULL,
        models              = list(
            list(
                name       = name,
                formula    = formula,
                design     = design,
                covariates = covariates
            )
        ),
        skipped_covariates  = tibble::tibble(
            variable = character(0),
            score    = numeric(0),
            reason   = character(0)
        ),
        sample_weights_flag = FALSE,
        sample_weights_vars = character(0)
    )
}
