build_adstock <- function(signal = NULL,
                          unit_ids = NULL,
                          x_preclean = NULL,
                          data = NULL,
                          signal_col = NULL,
                          unit_col = NULL,
                          controls = NULL,
                          preclean = FALSE,
                          preclean_method = c("mm", "s"),
                          nonnegative = FALSE,
                          residual_clip_multiplier = 2.5,
                          s_starts = 20L,
                          s_max_iter = 50L,
                          mm_max_iter = 50L,
                          mm_tukey_c = 4.685,
                          seed = 1234L,
                          rho = 0.85,
                          increment_method = c("plain", "huber", "tanh", "softsign", "adaptive_clip"),
                          clip_c = 1.5,
                          adaptive_upper = 0) {
  preclean_method <- match.arg(preclean_method)
  increment_method <- match.arg(increment_method)

  if (!is.null(data)) {
    if (!is.data.frame(data)) {
      stop("`data` must be a data.frame.", call. = FALSE)
    }
    if (!(is.character(signal_col) && length(signal_col) == 1L && signal_col %in% names(data))) {
      stop("`signal_col` must name a single signal column in `data`.", call. = FALSE)
    }
    if (!(is.character(unit_col) && length(unit_col) == 1L && unit_col %in% names(data))) {
      stop("`unit_col` must name a single unit-id column in `data`.", call. = FALSE)
    }

    if (is.null(controls)) {
      controls <- setdiff(names(data), c(signal_col, unit_col))
    }
    missing_controls <- setdiff(controls, names(data))
    if (length(missing_controls) > 0L) {
      stop("Missing control columns in `data`: ", paste(missing_controls, collapse = ", "), call. = FALSE)
    }

    used_cols <- unique(c(signal_col, unit_col, controls))
    mf <- stats::na.omit(data[, used_cols, drop = FALSE])
    signal <- as.numeric(mf[[signal_col]])
    unit_ids <- mf[[unit_col]]

    if (length(controls) > 0L) {
      x_preclean <- cbind("(Intercept)" = 1, as.matrix(mf[, controls, drop = FALSE]))
      control_names <- c("(Intercept)", controls)
    } else {
      x_preclean <- matrix(1, nrow = nrow(mf), ncol = 1)
      colnames(x_preclean) <- "(Intercept)"
      control_names <- "(Intercept)"
    }
  } else {
    if (is.null(signal) || is.null(unit_ids)) {
      stop("Provide either `signal` and `unit_ids`, or `data`, `signal_col`, and `unit_col`.", call. = FALSE)
    }
    signal <- as.numeric(signal)
    unit_ids <- as.integer(as.factor(unit_ids))
    if (is.null(x_preclean)) {
      x_preclean <- matrix(1, nrow = length(signal), ncol = 1)
      colnames(x_preclean) <- "(Intercept)"
    } else {
      x_preclean <- as.matrix(x_preclean)
      if (!("(Intercept)" %in% colnames(x_preclean))) {
        x_preclean <- cbind("(Intercept)" = 1, x_preclean)
      }
    }
    control_names <- colnames(x_preclean) %||% paste0("x", seq_len(ncol(x_preclean)))
  }

  if (!is.numeric(signal)) {
    stop("`signal` must be numeric.", call. = FALSE)
  }
  if (length(signal) != length(unit_ids)) {
    stop("`signal` and `unit_ids` must have the same length.", call. = FALSE)
  }
  if (!is.matrix(x_preclean) || !is.numeric(x_preclean)) {
    stop("`x_preclean` must coerce to a numeric matrix.", call. = FALSE)
  }
  if (nrow(x_preclean) != length(signal)) {
    stop("`x_preclean` must have the same number of rows as `signal`.", call. = FALSE)
  }

  fit <- .Call(
    rc_r_build_adstock,
    as.numeric(signal),
    as.integer(as.factor(unit_ids)),
    unname(x_preclean),
    list(
      preclean_enabled = isTRUE(preclean),
      nonnegative = isTRUE(nonnegative),
      use_mm = identical(preclean_method, "mm"),
      residual_clip_multiplier = as.numeric(residual_clip_multiplier),
      s_starts = as.integer(s_starts),
      s_max_iter = as.integer(s_max_iter),
      mm_max_iter = as.integer(mm_max_iter),
      mm_tukey_c = as.numeric(mm_tukey_c),
      seed = as.integer(seed),
      rho = as.numeric(rho),
      increment_method = increment_method,
      clip_c = as.numeric(clip_c),
      adaptive_upper = as.numeric(adaptive_upper)
    )
  )

  fit$call <- match.call()
  fit$control_names <- control_names
  fit$nobs <- length(signal)
  fit$signal <- signal
  fit$unit_ids <- unit_ids
  fit$x_preclean <- x_preclean
  class(fit) <- "robustcause_adstock"
  fit
}

print.robustcause_adstock <- function(x, ...) {
  cat("RobustCause adstock fit\n")
  cat("N:", x$nobs, "\n")
  cat("Preclean:", if (isTRUE(x$preclean_enabled)) x$preclean_method else "none", "\n")
  cat("Increment:", x$increment_method, "\n")
  cat("Rho:", format(signif(x$rho, 6)), "\n")
  cat("\nHead of stock:\n")
  print(utils::head(x$stock))
  invisible(x)
}

summary.robustcause_adstock <- function(object, ...) {
  delta <- object$signal - object$cleaned_signal
  structure(
    list(
      call = object$call,
      nobs = object$nobs,
      preclean_enabled = object$preclean_enabled,
      preclean_method = object$preclean_method,
      increment_method = object$increment_method,
      rho = object$rho,
      cleaning_shift_mean = mean(delta),
      cleaning_shift_max = max(abs(delta)),
      stock_summary = stats::setNames(summary(object$stock), names(summary(object$stock)))
    ),
    class = "summary.robustcause_adstock"
  )
}

print.summary.robustcause_adstock <- function(x, ...) {
  cat("RobustCause adstock summary\n\n")
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("N:", x$nobs, "\n")
  cat("Preclean:", if (isTRUE(x$preclean_enabled)) x$preclean_method else "none", "\n")
  cat("Increment:", x$increment_method, "\n")
  cat("Rho:", format(signif(x$rho, 6)), "\n")
  cat("Mean cleaning shift:", format(signif(x$cleaning_shift_mean, 6)), "\n")
  cat("Max absolute cleaning shift:", format(signif(x$cleaning_shift_max, 6)), "\n")
  cat("\nStock summary:\n")
  print(x$stock_summary)
  invisible(x)
}
