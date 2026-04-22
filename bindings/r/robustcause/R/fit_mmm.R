fit_mmm <- function(y = NULL,
                    unit_ids = NULL,
                    channels = NULL,
                    preclean_x = NULL,
                    controls = NULL,
                    data = NULL,
                    outcome = NULL,
                    unit = NULL,
                    channel_cols = NULL,
                    preclean_cols = NULL,
                    control_cols = NULL,
                    fit_method = c("huber", "ols", "s", "mm"),
                    rho = 0.85,
                    increment_method = "adaptive_clip",
                    hill_lambda = 0.10,
                    preclean = TRUE,
                    preclean_method = "mm",
                    nonnegative = TRUE,
                    residual_clip_multiplier = 2.5,
                    clip_c = 1.5,
                    adaptive_upper = 0,
                    apply_hill = TRUE,
                    huber_c = 1.345,
                    huber_max_iter = 75L,
                    s_starts = 20L,
                    s_max_iter = 50L,
                    mm_max_iter = 50L,
                    mm_tukey_c = 4.685,
                    seed = 1234L) {
  fit_method <- match.arg(fit_method)

  if (!is.null(data)) {
    if (!is.data.frame(data)) {
      stop("`data` must be a data.frame.", call. = FALSE)
    }
    if (!(is.character(outcome) && length(outcome) == 1L && outcome %in% names(data))) {
      stop("`outcome` must name a single outcome column in `data`.", call. = FALSE)
    }
    if (!(is.character(unit) && length(unit) == 1L && unit %in% names(data))) {
      stop("`unit` must name a single unit column in `data`.", call. = FALSE)
    }
    if (is.null(channel_cols) || !all(channel_cols %in% names(data))) {
      stop("`channel_cols` must name one or more media channel columns in `data`.", call. = FALSE)
    }
    preclean_cols <- preclean_cols %||% channel_cols
    control_cols <- control_cols %||% setdiff(names(data), c(outcome, unit, channel_cols))
    used_cols <- unique(c(outcome, unit, channel_cols, preclean_cols, control_cols))
    mf <- stats::na.omit(data[, used_cols, drop = FALSE])
    y <- as.numeric(mf[[outcome]])
    unit_ids <- as.integer(as.factor(mf[[unit]]))
    channels <- as.matrix(mf[, channel_cols, drop = FALSE])
    preclean_x <- cbind("(Intercept)" = 1, as.matrix(mf[, preclean_cols, drop = FALSE]))
    controls <- as.matrix(mf[, control_cols, drop = FALSE])
    channel_names <- channel_cols
    control_names <- colnames(controls)
  } else {
    if (is.null(y) || is.null(unit_ids) || is.null(channels) || is.null(preclean_x) || is.null(controls)) {
      stop("Provide either the matrix inputs or the data-frame interface.", call. = FALSE)
    }
    y <- as.numeric(y)
    unit_ids <- as.integer(as.factor(unit_ids))
    channels <- as.matrix(channels)
    preclean_x <- as.matrix(preclean_x)
    controls <- as.matrix(controls)
    channel_names <- colnames(channels) %||% paste0("channel_", seq_len(ncol(channels)))
    control_names <- colnames(controls)
  }

  if (!is.numeric(channels) || !is.numeric(preclean_x) || !is.numeric(controls)) {
    stop("`channels`, `preclean_x`, and `controls` must be numeric.", call. = FALSE)
  }
  n <- length(y)
  if (nrow(channels) != n || nrow(preclean_x) != n || nrow(controls) != n || length(unit_ids) != n) {
    stop("All MMM inputs must have the same number of observations.", call. = FALSE)
  }
  if (ncol(channels) < 1L) {
    stop("`channels` must have at least one column.", call. = FALSE)
  }

  k <- ncol(channels)
  recycle_arg <- function(value, name, choices = NULL) {
    if (length(value) == 1L) {
      value <- rep(value, k)
    }
    if (length(value) != k) {
      stop("`", name, "` must have length 1 or match the number of channels.", call. = FALSE)
    }
    if (!is.null(choices) && !all(value %in% choices)) {
      stop("`", name, "` has invalid values.", call. = FALSE)
    }
    value
  }

  rho <- as.numeric(recycle_arg(rho, "rho"))
  increment_method <- recycle_arg(increment_method, "increment_method",
                                  c("plain", "huber", "tanh", "softsign", "adaptive_clip"))
  hill_lambda <- as.numeric(recycle_arg(hill_lambda, "hill_lambda"))
  preclean <- as.logical(recycle_arg(preclean, "preclean"))
  preclean_method <- recycle_arg(preclean_method, "preclean_method", c("mm", "s"))
  nonnegative <- as.logical(recycle_arg(nonnegative, "nonnegative"))
  residual_clip_multiplier <- as.numeric(recycle_arg(residual_clip_multiplier, "residual_clip_multiplier"))
  clip_c <- as.numeric(recycle_arg(clip_c, "clip_c"))
  adaptive_upper <- as.numeric(recycle_arg(adaptive_upper, "adaptive_upper"))
  apply_hill <- as.logical(recycle_arg(apply_hill, "apply_hill"))
  preclean_seed <- as.integer(seed + seq_len(k) - 1L)

  fit <- .Call(
    rc_r_fit_mmm,
    y,
    unit_ids,
    unname(channels),
    unname(preclean_x),
    unname(controls),
    list(
      fit_method = fit_method,
      huber_c = as.numeric(huber_c),
      huber_max_iter = as.integer(huber_max_iter),
      s_starts = as.integer(s_starts),
      s_max_iter = as.integer(s_max_iter),
      mm_max_iter = as.integer(mm_max_iter),
      mm_tukey_c = as.numeric(mm_tukey_c),
      seed = as.integer(seed),
      channel_names = as.character(channel_names),
      preclean_enabled = preclean,
      preclean_method = as.character(preclean_method),
      nonnegative = nonnegative,
      residual_clip_multiplier = residual_clip_multiplier,
      preclean_seed = preclean_seed,
      rho = rho,
      increment_method = as.character(increment_method),
      clip_c = clip_c,
      adaptive_upper = adaptive_upper,
      hill_lambda = hill_lambda,
      apply_hill = apply_hill
    )
  )

  colnames(fit$cleaned_signals) <- channel_names
  colnames(fit$stocks) <- channel_names
  colnames(fit$transformed_channels) <- channel_names
  colnames(fit$design_matrix) <- c("(Intercept)", channel_names, control_names %||% paste0("control_", seq_len(ncol(controls))))
  names(fit$channel_coef) <- channel_names
  names(fit$full_coef) <- colnames(fit$design_matrix)
  fit$call <- match.call()
  fit$nobs <- n
  fit$channel_names <- channel_names
  fit$control_names <- control_names
  class(fit) <- "robustcause_mmm"
  fit
}

print.robustcause_mmm <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

summary.robustcause_mmm <- function(object, ...) {
  coef_names <- names(object$full_coef)
  roles <- rep("control", length(object$full_coef))
  roles[coef_names == "(Intercept)"] <- "intercept"
  roles[coef_names %in% object$channel_names] <- "channel"
  coefficients <- data.frame(
    term = coef_names,
    role = roles,
    estimate = unname(object$full_coef),
    row.names = NULL,
    check.names = FALSE
  )
  channel_table <- data.frame(
    channel = object$channel_names,
    estimate = unname(object$channel_coef),
    cleaned_mean = colMeans(object$cleaned_signals),
    stock_mean = colMeans(object$stocks),
    transformed_mean = colMeans(object$transformed_channels),
    row.names = NULL,
    check.names = FALSE
  )
  structure(
    list(
      call = object$call,
      ok = object$ok,
      nobs = object$nobs,
      fit_method = object$fit_method,
      coefficients = coefficients,
      channel_table = channel_table,
      design_columns = colnames(object$design_matrix),
      design_ncol = ncol(object$design_matrix),
      regression_scale = object$regression_scale
    ),
    class = "summary.robustcause_mmm"
  )
}

print.summary.robustcause_mmm <- function(x, ...) {
  cat("RobustCause MMM summary\n\n")
  cat("Call:\n")
  print(x$call)
  cat("\nFit ok:", x$ok, "\n")
  cat("N:", x$nobs, "\n")
  cat("Fit method:", x$fit_method, "\n")
  cat("Design columns:", x$design_ncol, "\n")
  cat("Regression scale:", format(signif(x$regression_scale, 6)), "\n")
  cat("\nChannel summary:\n")
  print(x$channel_table, row.names = FALSE)
  cat("\nCoefficient table:\n")
  print(x$coefficients, row.names = FALSE)
  invisible(x)
}
