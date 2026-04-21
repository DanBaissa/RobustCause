fit_mm_dml <- function(x = NULL,
                       d = NULL,
                       y = NULL,
                       data = NULL,
                       outcome = NULL,
                       treatment = NULL,
                       controls = NULL,
                       learner = c("lasso", "elastic_net", "random_forest", "hist_gradient_boosting"),
                       folds = 2L,
                       seed = 123L,
                       mm_c = 4.685,
                       max_iter = 100L,
                       tolerance = 1e-7,
                       ci_level = 0.95) {
  learner <- match.arg(learner)

  if (!is.null(data) || !is.null(outcome) || !is.null(treatment) || !is.null(controls)) {
    if (is.null(data) || !is.data.frame(data)) {
      stop("`data` must be a data.frame when using the named-column interface.", call. = FALSE)
    }
    if (!(is.character(outcome) && length(outcome) == 1L && outcome %in% names(data))) {
      stop("`outcome` must name a single column in `data`.", call. = FALSE)
    }
    if (!(is.character(treatment) && length(treatment) == 1L && treatment %in% names(data))) {
      stop("`treatment` must name a single column in `data`.", call. = FALSE)
    }
    controls <- controls %||% setdiff(names(data), c(outcome, treatment))
    missing_controls <- setdiff(controls, names(data))
    if (length(missing_controls) > 0L) {
      stop("Missing control columns in `data`: ", paste(missing_controls, collapse = ", "), call. = FALSE)
    }

    mf <- stats::na.omit(data[, unique(c(outcome, treatment, controls)), drop = FALSE])
    x <- if (length(controls) > 0L) {
      as.matrix(mf[, controls, drop = FALSE])
    } else {
      matrix(0, nrow = nrow(mf), ncol = 1)
    }
    d <- as.numeric(mf[[treatment]])
    y <- as.numeric(mf[[outcome]])
    control_names <- if (length(controls) > 0L) controls else "(none)"
  } else {
    if (is.null(x) || is.null(d) || is.null(y)) {
      stop("Provide either `x`, `d`, and `y`, or `data`, `outcome`, and `treatment`.", call. = FALSE)
    }
    x <- as.matrix(x)
    d <- as.numeric(d)
    y <- as.numeric(y)
    control_names <- colnames(x)
  }

  if (!is.numeric(x) || !is.numeric(d) || !is.numeric(y)) {
    stop("`x`, `d`, and `y` must be numeric.", call. = FALSE)
  }
  if (nrow(x) != length(d) || length(d) != length(y)) {
    stop("`x`, `d`, and `y` must have the same number of observations.", call. = FALSE)
  }

  fit <- .Call(
    rc_r_fit_mm_dml,
    unname(x),
    d,
    y,
    learner,
    as.integer(folds),
    as.integer(seed),
    as.numeric(mm_c),
    as.integer(max_iter),
    as.numeric(tolerance),
    as.numeric(ci_level)
  )

  fit$call <- match.call()
  fit$learner <- learner
  fit$folds <- as.integer(folds)
  fit$seed <- as.integer(seed)
  fit$controls <- control_names
  fit$nobs <- nrow(x)
  fit$ci_level <- ci_level
  class(fit) <- "robustcause_mm_dml"
  fit
}

summary.robustcause_mm_dml <- function(object, ...) {
  coef_table <- cbind(
    Estimate = c(object$intercept_hat, object$tau_hat),
    `Std. Error` = c(NA_real_, object$std_error),
    lower = c(NA_real_, object$ci_lower),
    upper = c(NA_real_, object$ci_upper)
  )
  rownames(coef_table) <- c("(Intercept)", "treatment")

  structure(
    list(
      call = object$call,
      coefficients = coef_table,
      learner = object$learner,
      folds = object$folds,
      nobs = object$nobs,
      converged = object$converged,
      iterations = object$iterations,
      objective = object$objective,
      kept_fraction = object$kept_fraction
    ),
    class = "summary.robustcause_mm_dml"
  )
}

print.summary.robustcause_mm_dml <- function(x, ...) {
  cat(
    sprintf(
      "RobustCause MM-DML fit | Learner: %s | Folds: %d | N: %d\n\n",
      x$learner,
      x$folds,
      x$nobs
    )
  )
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
  cat("\n")
  cat(sprintf("Converged: %s | Iterations: %s | Objective: %.6f | Kept fraction: %.3f\n",
              ifelse(isTRUE(x$converged), "TRUE", "FALSE"),
              as.character(x$iterations),
              x$objective,
              x$kept_fraction))
  invisible(x)
}

print.robustcause_mm_dml <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}
