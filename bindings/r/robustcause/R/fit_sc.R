fit_sc <- function(outcomes = NULL,
                   treated_unit = NULL,
                   treatment_start = NULL,
                   donors = NULL,
                   unit = NULL,
                   time = NULL,
                   outcome = NULL,
                   data = NULL,
                   predictors = NULL,
                   predictor_agg = c("mean", "last_pre"),
                   predictor_weights = NULL,
                   method = c("standard", "mm"),
                   run_placebos = TRUE,
                   na.action = c("omit", "fail"),
                   maxit = 500L,
                   tol = 1e-8,
                   ridge = 1e-10,
                   startup_maxit = 75L,
                   subproblem_maxit = 500L,
                   subproblem_tol = 1e-8,
                   l1_smoothing = 1e-6,
                   tukey_c = 4.685,
                   min_scale = 1e-8,
                   min_time_weight = 1e-8) {
  method <- match.arg(method)
  predictor_agg <- match.arg(predictor_agg)
  na.action <- match.arg(na.action)

  prepared <- .prepare_sc_inputs(
    outcomes = outcomes,
    treated_unit = treated_unit,
    treatment_start = treatment_start,
    donors = donors,
    unit = unit,
    time = time,
    outcome = outcome,
    data = data,
    predictors = predictors,
    predictor_agg = predictor_agg,
    predictor_weights = predictor_weights,
    na.action = na.action
  )

  sc_ctl <- list(
    maxit = as.integer(maxit),
    tol = as.numeric(tol),
    ridge = as.numeric(ridge)
  )
  mm_ctl <- list(
    startup_maxit = as.integer(startup_maxit),
    maxit = as.integer(maxit),
    subproblem_maxit = as.integer(subproblem_maxit),
    tol = as.numeric(tol),
    subproblem_tol = as.numeric(subproblem_tol),
    l1_smoothing = as.numeric(l1_smoothing),
    tukey_c = as.numeric(tukey_c),
    min_scale = as.numeric(min_scale),
    ridge = as.numeric(ridge),
    min_time_weight = as.numeric(min_time_weight)
  )

  fit <- if (identical(method, "mm")) {
    .Call(
      rc_r_fit_mm_sc,
      prepared$treated_pre,
      prepared$donors_pre,
      prepared$treated_post,
      prepared$donors_post,
      mm_ctl$startup_maxit,
      mm_ctl$maxit,
      mm_ctl$subproblem_maxit,
      mm_ctl$tol,
      mm_ctl$subproblem_tol,
      mm_ctl$l1_smoothing,
      mm_ctl$tukey_c,
      mm_ctl$min_scale,
      mm_ctl$ridge,
      mm_ctl$min_time_weight,
      prepared$predictor_treated,
      prepared$predictor_donors,
      prepared$predictor_weights
    )
  } else {
    .Call(
      rc_r_fit_sc,
      prepared$treated_pre,
      prepared$donors_pre,
      prepared$treated_post,
      prepared$donors_post,
      sc_ctl$maxit,
      sc_ctl$tol,
      sc_ctl$ridge,
      prepared$predictor_treated,
      prepared$predictor_donors,
      prepared$predictor_weights
    )
  }

  fit <- .normalize_native_sc_fit(fit, prepared, method)
  fit$method <- method
  fit$treated_unit <- prepared$treated_unit
  fit$donors <- prepared$donors
  fit$time_index <- prepared$time_index
  fit$pre_periods <- prepared$pre_periods
  fit$post_periods <- prepared$post_periods
  fit$predictors <- prepared$predictors
  fit$predictor_agg <- predictor_agg
  fit$call <- match.call()
  fit$control <- if (identical(method, "mm")) mm_ctl else sc_ctl
  fit$data_info <- prepared$data_info
  if (identical(method, "mm")) {
    fit$pre_rmspe <- fit$mm_rmspe
  }

  if (isTRUE(run_placebos) && ncol(prepared$outcomes_matrix) > 1L) {
    placebo <- .Call(
      rc_r_fit_sc_placebos,
      prepared$outcomes_matrix,
      as.integer(prepared$treatment_start_index),
      sc_ctl$maxit,
      sc_ctl$tol,
      sc_ctl$ridge,
      mm_ctl$startup_maxit,
      mm_ctl$subproblem_maxit,
      mm_ctl$subproblem_tol,
      mm_ctl$l1_smoothing,
      mm_ctl$tukey_c,
      mm_ctl$min_scale,
      mm_ctl$min_time_weight,
      identical(method, "mm")
    )
    fit$placebo <- .finalize_sc_placebos(placebo, prepared)
    fit$placebos <- fit$placebo
    fit$inference <- .summarize_sc_inference(fit, fit$placebo)
  } else {
    fit$placebo <- NULL
    fit$placebos <- NULL
    fit$inference <- NULL
  }

  class(fit) <- "robustcause_sc"
  fit
}

.prepare_sc_inputs <- function(outcomes,
                               treated_unit,
                               treatment_start,
                               donors,
                               unit,
                               time,
                               outcome,
                               data,
                               predictors,
                               predictor_agg,
                               predictor_weights,
                               na.action) {
  if (!is.null(data)) {
    if (!is.data.frame(data)) {
      stop("`data` must be a data.frame.", call. = FALSE)
    }
    if (!(is.character(unit) && length(unit) == 1L && unit %in% names(data))) {
      stop("`unit` must name a column in `data`.", call. = FALSE)
    }
    if (!(is.character(time) && length(time) == 1L && time %in% names(data))) {
      stop("`time` must name a column in `data`.", call. = FALSE)
    }
    if (!(is.character(outcome) && length(outcome) == 1L && outcome %in% names(data))) {
      stop("`outcome` must name a column in `data`.", call. = FALSE)
    }
    if (is.null(treated_unit) || length(treated_unit) != 1L) {
      stop("`treated_unit` must identify one treated unit for the data-frame interface.", call. = FALSE)
    }
    if (is.null(treatment_start) || length(treatment_start) != 1L) {
      stop("`treatment_start` must identify the first treated period.", call. = FALSE)
    }
    predictors <- predictors %||% character(0)
    missing_predictors <- setdiff(predictors, names(data))
    if (length(missing_predictors) > 0L) {
      stop("Missing predictor columns: ", paste(missing_predictors, collapse = ", "), call. = FALSE)
    }
    needed <- unique(c(unit, time, outcome, predictors))
    df <- data[, needed, drop = FALSE]
    if (identical(na.action, "omit")) {
      df <- stats::na.omit(df)
    } else if (any(!stats::complete.cases(df))) {
      stop("Missing values detected; use `na.action = 'omit'` to drop incomplete rows.", call. = FALSE)
    }

    unit_levels <- unique(df[[unit]])
    if (!(treated_unit %in% unit_levels)) {
      stop("`treated_unit` was not found in `data`.", call. = FALSE)
    }
    donor_units <- donors %||% setdiff(unit_levels, treated_unit)
    if (length(donor_units) < 1L) {
      stop("Need at least one donor unit.", call. = FALSE)
    }
    donor_units <- setdiff(donor_units, treated_unit)
    missing_donors <- setdiff(donor_units, unit_levels)
    if (length(missing_donors) > 0L) {
      stop("Donor units not found in `data`: ", paste(missing_donors, collapse = ", "), call. = FALSE)
    }

    time_index <- sort(unique(df[[time]]))
    if (!(treatment_start %in% time_index)) {
      stop("`treatment_start` was not found in the time index.", call. = FALSE)
    }
    wide <- reshape(
      df[, c(unit, time, outcome), drop = FALSE],
      idvar = time,
      timevar = unit,
      direction = "wide"
    )
    wide <- wide[match(time_index, wide[[time]]), , drop = FALSE]
    outcome_cols <- paste0(outcome, ".", c(treated_unit, donor_units))
    missing_outcome_cols <- setdiff(outcome_cols, names(wide))
    if (length(missing_outcome_cols) > 0L) {
      stop("Panel is unbalanced or missing outcome observations for some unit-time cells.", call. = FALSE)
    }
    outcomes_matrix <- as.matrix(wide[, outcome_cols, drop = FALSE])
    colnames(outcomes_matrix) <- c(as.character(treated_unit), as.character(donor_units))

    treatment_start_index <- match(treatment_start, time_index)
    pre_idx <- seq_len(treatment_start_index - 1L)
    post_idx <- seq.int(treatment_start_index, length(time_index))
    if (length(pre_idx) < 1L || length(post_idx) < 1L) {
      stop("Need at least one pre- and one post-treatment period.", call. = FALSE)
    }

    predictor_info <- .build_sc_predictors_from_data(
      df = df,
      unit = unit,
      time = time,
      treated_unit = treated_unit,
      donor_units = donor_units,
      time_index = time_index,
      pre_idx = pre_idx,
      predictors = predictors,
      predictor_agg = predictor_agg,
      predictor_weights = predictor_weights
    )

    return(list(
      outcomes_matrix = unname(outcomes_matrix),
      treated_pre = as.numeric(outcomes_matrix[pre_idx, 1]),
      donors_pre = unname(as.matrix(outcomes_matrix[pre_idx, -1, drop = FALSE])),
      treated_post = as.numeric(outcomes_matrix[post_idx, 1]),
      donors_post = unname(as.matrix(outcomes_matrix[post_idx, -1, drop = FALSE])),
      predictor_treated = predictor_info$predictor_treated,
      predictor_donors = predictor_info$predictor_donors,
      predictor_weights = predictor_info$predictor_weights,
      treated_unit = treated_unit,
      donors = donor_units,
      time_index = time_index,
      pre_periods = time_index[pre_idx],
      post_periods = time_index[post_idx],
      treatment_start_index = treatment_start_index - 1L,
      predictors = predictors,
      data_info = list(interface = "data.frame")
    ))
  }

  if (is.null(outcomes)) {
    stop("Provide either `outcomes` or panel `data` with `unit`, `time`, and `outcome`.", call. = FALSE)
  }
  outcomes_matrix <- as.matrix(outcomes)
  if (!is.numeric(outcomes_matrix)) {
    stop("`outcomes` must coerce to a numeric matrix.", call. = FALSE)
  }
  if (identical(na.action, "omit")) {
    keep <- stats::complete.cases(outcomes_matrix)
    outcomes_matrix <- outcomes_matrix[keep, , drop = FALSE]
  } else if (any(!is.finite(outcomes_matrix))) {
    stop("`outcomes` contains missing or non-finite values; use `na.action = 'omit'` to drop incomplete periods.", call. = FALSE)
  }
  if (ncol(outcomes_matrix) < 2L) {
    stop("`outcomes` must have at least two columns.", call. = FALSE)
  }
  col_names <- colnames(outcomes_matrix) %||% as.character(seq_len(ncol(outcomes_matrix)))
  treated_col <- .resolve_sc_unit_index(treated_unit %||% col_names[[1]], col_names)
  donor_cols <- if (is.null(donors)) setdiff(seq_len(ncol(outcomes_matrix)), treated_col) else .resolve_sc_donor_indices(donors, col_names, treated_col)
  if (length(donor_cols) < 1L) {
    stop("Need at least one donor column.", call. = FALSE)
  }
  treatment_start_index <- .resolve_sc_treatment_start(treatment_start, nrow(outcomes_matrix))
  pre_idx <- seq_len(treatment_start_index - 1L)
  post_idx <- seq.int(treatment_start_index, nrow(outcomes_matrix))
  predictor_info <- .build_sc_predictors_from_matrix(
    predictors = predictors,
    predictor_weights = predictor_weights,
    n_units = ncol(outcomes_matrix),
    treated_col = treated_col,
    donor_cols = donor_cols
  )

  list(
    outcomes_matrix = unname(as.matrix(outcomes_matrix[, c(treated_col, donor_cols), drop = FALSE])),
    treated_pre = as.numeric(outcomes_matrix[pre_idx, treated_col]),
    donors_pre = unname(as.matrix(outcomes_matrix[pre_idx, donor_cols, drop = FALSE])),
    treated_post = as.numeric(outcomes_matrix[post_idx, treated_col]),
    donors_post = unname(as.matrix(outcomes_matrix[post_idx, donor_cols, drop = FALSE])),
    predictor_treated = predictor_info$predictor_treated,
    predictor_donors = predictor_info$predictor_donors,
    predictor_weights = predictor_info$predictor_weights,
    treated_unit = col_names[[treated_col]],
    donors = col_names[donor_cols],
    time_index = seq_len(nrow(outcomes_matrix)),
    pre_periods = pre_idx,
    post_periods = post_idx,
    treatment_start_index = treatment_start_index - 1L,
    predictors = character(0),
    data_info = list(interface = "matrix")
  )
}

.resolve_sc_unit_index <- function(treated_unit, col_names) {
  if (is.numeric(treated_unit)) {
    idx <- as.integer(treated_unit[[1]])
  } else {
    idx <- match(as.character(treated_unit[[1]]), col_names)
  }
  if (is.na(idx) || idx < 1L || idx > length(col_names)) {
    stop("`treated_unit` could not be resolved.", call. = FALSE)
  }
  idx
}

.resolve_sc_donor_indices <- function(donors, col_names, treated_col) {
  idx <- if (is.numeric(donors)) as.integer(donors) else match(as.character(donors), col_names)
  if (anyNA(idx) || any(idx < 1L) || any(idx > length(col_names))) {
    stop("`donors` contains unknown units/columns.", call. = FALSE)
  }
  idx <- unique(setdiff(idx, treated_col))
  idx
}

.resolve_sc_treatment_start <- function(treatment_start, n_periods) {
  if (length(treatment_start) != 1L || is.na(treatment_start)) {
    stop("`treatment_start` must be a single period index.", call. = FALSE)
  }
  idx <- as.integer(treatment_start)
  if (idx <= 1L || idx > n_periods) {
    stop("`treatment_start` must leave at least one pre- and one post-treatment period.", call. = FALSE)
  }
  idx
}

.build_sc_predictors_from_matrix <- function(predictors, predictor_weights, n_units, treated_col, donor_cols) {
  if (is.null(predictors)) {
    return(list(
      predictor_treated = numeric(),
      predictor_donors = matrix(numeric(), nrow = 0L, ncol = length(donor_cols)),
      predictor_weights = numeric()
    ))
  }
  pred <- as.matrix(predictors)
  if (!is.numeric(pred)) {
    stop("`predictors` must coerce to a numeric matrix.", call. = FALSE)
  }
  if (ncol(pred) != n_units) {
    stop("For matrix input, `predictors` must have one column per unit in `outcomes`.", call. = FALSE)
  }
  list(
    predictor_treated = as.numeric(pred[, treated_col]),
    predictor_donors = unname(as.matrix(pred[, donor_cols, drop = FALSE])),
    predictor_weights = .normalize_sc_predictor_weights(predictor_weights, nrow(pred))
  )
}

.build_sc_predictors_from_data <- function(df,
                                           unit,
                                           time,
                                           treated_unit,
                                           donor_units,
                                           time_index,
                                           pre_idx,
                                           predictors,
                                           predictor_agg,
                                           predictor_weights) {
  if (length(predictors) == 0L) {
    return(list(
      predictor_treated = numeric(),
      predictor_donors = matrix(numeric(), nrow = 0L, ncol = length(donor_units)),
      predictor_weights = numeric()
    ))
  }
  pre_times <- time_index[pre_idx]
  pred_mat <- matrix(NA_real_, nrow = length(predictors), ncol = length(donor_units) + 1L)
  colnames(pred_mat) <- c(as.character(treated_unit), as.character(donor_units))
  rownames(pred_mat) <- predictors
  for (j in seq_along(colnames(pred_mat))) {
    unit_value <- colnames(pred_mat)[[j]]
    rows <- df[[unit]] == unit_value & df[[time]] %in% pre_times
    slice <- df[rows, predictors, drop = FALSE]
    if (nrow(slice) == 0L) {
      stop("No pre-treatment rows available to aggregate predictors for unit ", unit_value, ".", call. = FALSE)
    }
    agg <- if (identical(predictor_agg, "mean")) {
      vapply(slice, function(x) mean(x, na.rm = TRUE), numeric(1))
    } else {
      ordered <- df[rows, c(time, predictors), drop = FALSE]
      ordered <- ordered[order(ordered[[time]]), , drop = FALSE]
      as.numeric(ordered[nrow(ordered), predictors, drop = TRUE])
    }
    pred_mat[, j] <- agg
  }
  list(
    predictor_treated = as.numeric(pred_mat[, 1]),
    predictor_donors = unname(as.matrix(pred_mat[, -1, drop = FALSE])),
    predictor_weights = .normalize_sc_predictor_weights(predictor_weights, nrow(pred_mat))
  )
}

.normalize_sc_predictor_weights <- function(predictor_weights, n_predictors) {
  if (n_predictors == 0L) {
    return(numeric())
  }
  if (is.null(predictor_weights)) {
    return(rep(1, n_predictors))
  }
  weights <- as.numeric(predictor_weights)
  if (length(weights) != n_predictors || any(!is.finite(weights)) || any(weights <= 0)) {
    stop("`predictor_weights` must be a positive numeric vector matching the number of predictors.", call. = FALSE)
  }
  weights
}

.normalize_native_sc_result <- function(fit, donors, pre_periods, post_periods) {
  pre_names <- as.character(pre_periods)
  post_names <- as.character(post_periods)
  fit$weights <- stats::setNames(as.numeric(fit$weights), donors)
  fit$synthetic_pre <- stats::setNames(as.numeric(fit$synthetic_pre), pre_names)
  fit$synthetic_post <- stats::setNames(as.numeric(fit$synthetic_post), post_names)
  fit$pre_residuals <- stats::setNames(as.numeric(fit$pre_residuals), pre_names)
  fit$post_gaps <- stats::setNames(as.numeric(fit$post_gaps), post_names)
  fit
}

.build_sc_path <- function(pre_periods,
                           post_periods,
                           treated_pre,
                           treated_post,
                           synthetic_pre,
                           synthetic_post,
                           pre_residuals,
                           post_gaps) {
  data.frame(
    time = c(pre_periods, post_periods),
    period = c(rep("pre", length(pre_periods)), rep("post", length(post_periods))),
    treated = c(as.numeric(treated_pre), as.numeric(treated_post)),
    synthetic = c(as.numeric(synthetic_pre), as.numeric(synthetic_post)),
    gap = c(as.numeric(pre_residuals), as.numeric(post_gaps)),
    row.names = NULL
  )
}

.normalize_native_sc_fit <- function(fit, prepared, method) {
  fit <- .normalize_native_sc_result(
    fit = fit,
    donors = prepared$donors,
    pre_periods = prepared$pre_periods,
    post_periods = prepared$post_periods
  )
  fit$path <- .build_sc_path(
    pre_periods = prepared$pre_periods,
    post_periods = prepared$post_periods,
    treated_pre = prepared$treated_pre,
    treated_post = prepared$treated_post,
    synthetic_pre = fit$synthetic_pre,
    synthetic_post = fit$synthetic_post,
    pre_residuals = fit$pre_residuals,
    post_gaps = fit$post_gaps
  )

  if (identical(method, "mm")) {
    fit$startup_weights <- stats::setNames(as.numeric(fit$startup_weights), prepared$donors)
    fit$robust_time_weights <- stats::setNames(as.numeric(fit$robust_time_weights), as.character(prepared$pre_periods))
    fit$downweighted_periods <- as.integer(fit$downweighted_periods)
    fit$standard_fit <- .normalize_native_sc_result(
      fit = fit$standard_fit,
      donors = prepared$donors,
      pre_periods = prepared$pre_periods,
      post_periods = prepared$post_periods
    )
  }

  fit
}

.finalize_sc_placebos <- function(placebo, prepared) {
  placebo$unit_names <- c(prepared$treated_unit, prepared$donors)[placebo$unit_indices + 1L]
  rownames(placebo$weight_matrix) <- placebo$unit_names
  rownames(placebo$post_gap_matrix) <- placebo$unit_names
  colnames(placebo$post_gap_matrix) <- as.character(prepared$post_periods)
  placebo
}

.summarize_sc_inference <- function(fit, placebo) {
  treated_name <- fit$treated_unit
  idx <- match(treated_name, placebo$unit_names)
  treated_ratio <- placebo$rmspe_ratio[[idx]]
  one_sided_p <- mean(placebo$rmspe_ratio >= treated_ratio)
  avg_post_gap <- mean(abs(fit$post_gaps))
  placebo_avg_post_gap <- rowMeans(abs(placebo$post_gap_matrix))
  gap_p <- mean(placebo_avg_post_gap >= avg_post_gap)
  list(
    treated_rmspe_ratio = treated_ratio,
    placebo_p_value_rmspe_ratio = one_sided_p,
    average_absolute_post_gap = avg_post_gap,
    placebo_p_value_avg_gap = gap_p
  )
}

summary.robustcause_sc <- function(object, ...) {
  list_out <- list(
    call = object$call,
    method = object$method,
    treated_unit = object$treated_unit,
    n_donors = length(object$donors),
    pre_rmspe = object$pre_rmspe %||% object$mm_rmspe,
    post_mean_gap = mean(object$post_gaps),
    post_mean_abs_gap = mean(abs(object$post_gaps)),
    weight_herfindahl = object$weight_herfindahl,
    max_weight = object$max_weight,
    effective_donors = object$effective_donors,
    converged = object$converged,
    iterations = object$iterations,
    weights = stats::setNames(as.numeric(object$weights), object$donors),
    placebo_inference = object$inference
  )
  if (identical(object$method, "mm")) {
    list_out$standard_rmspe <- object$standard_rmspe
    list_out$startup_rmspe <- object$startup_rmspe
    list_out$effective_periods <- object$effective_periods
    list_out$downweighted_periods <- object$pre_periods[object$downweighted_periods + 1L]
  }
  class(list_out) <- "summary.robustcause_sc"
  list_out
}

print.summary.robustcause_sc <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat(sprintf("RobustCause synthetic control | Method: %s | Treated unit: %s\n\n", x$method, x$treated_unit))
  cat("Call:\n")
  print(x$call)
  cat("\nFit diagnostics:\n")
  cat(sprintf("  Pre-RMSPE: %s\n", format(signif(x$pre_rmspe, digits))))
  cat(sprintf("  Post mean gap: %s\n", format(signif(x$post_mean_gap, digits))))
  cat(sprintf("  Post mean abs gap: %s\n", format(signif(x$post_mean_abs_gap, digits))))
  cat(sprintf("  Effective donors: %s\n", format(signif(x$effective_donors, digits))))
  cat(sprintf("  Max donor weight: %s\n", format(signif(x$max_weight, digits))))
  cat(sprintf("  Converged: %s | Iterations: %s\n", x$converged, x$iterations))
  if (!is.null(x$standard_rmspe)) {
    cat(sprintf("  Standard SC RMSPE: %s | Startup RMSPE: %s\n",
                format(signif(x$standard_rmspe, digits)),
                format(signif(x$startup_rmspe, digits))))
  }
  cat("\nWeights:\n")
  print(round(x$weights, digits))
  if (!is.null(x$placebo_inference)) {
    cat("\nPlacebo inference:\n")
    cat(sprintf("  RMSPE ratio p-value: %s\n", format(signif(x$placebo_inference$placebo_p_value_rmspe_ratio, digits))))
    cat(sprintf("  Avg |gap| p-value: %s\n", format(signif(x$placebo_inference$placebo_p_value_avg_gap, digits))))
  }
  invisible(x)
}

print.robustcause_sc <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

coef.robustcause_sc <- function(object, ...) {
  object$weights
}

fitted.robustcause_sc <- function(object, ...) {
  stats::setNames(object$path$synthetic, object$path$time)
}

residuals.robustcause_sc <- function(object, ...) {
  stats::setNames(object$path$gap, object$path$time)
}

predict.robustcause_sc <- function(object,
                                   newdata = NULL,
                                   type = c("synthetic", "gap"),
                                   period = c("all", "pre", "post"),
                                   ...) {
  if (!is.null(newdata)) {
    stop("`predict.robustcause_sc` does not support `newdata`; synthetic-control predictions are tied to the fitted panel.", call. = FALSE)
  }
  type <- match.arg(type)
  period <- match.arg(period)
  out <- switch(
    type,
    synthetic = object$path$synthetic,
    gap = object$path$gap
  )
  if (identical(period, "pre")) {
    out <- out[object$path$period == "pre"]
    names(out) <- object$path$time[object$path$period == "pre"]
    return(out)
  }
  if (identical(period, "post")) {
    out <- out[object$path$period == "post"]
    names(out) <- object$path$time[object$path$period == "post"]
    return(out)
  }
  stats::setNames(out, object$path$time)
}

as.data.frame.robustcause_sc <- function(x, ...) {
  x$path
}

plot.robustcause_sc <- function(x, type = c("gaps", "weights", "placebo_gaps"), ...) {
  type <- match.arg(type)
  if (identical(type, "weights")) {
    graphics::barplot(stats::setNames(x$weights, x$donors), las = 2, ylab = "Weight", main = "Synthetic-control donor weights", ...)
    return(invisible(x))
  }
  if (identical(type, "placebo_gaps")) {
    if (is.null(x$placebo)) {
      stop("No placebo output is attached to this fit; rerun with `run_placebos = TRUE`.", call. = FALSE)
    }
    time_labels <- colnames(x$placebo$post_gap_matrix)
    time_points <- seq_along(time_labels)
    graphics::matplot(
      x = time_points,
      y = t(x$placebo$post_gap_matrix),
      type = "l",
      lty = 1,
      col = grDevices::adjustcolor("gray60", alpha.f = 0.75),
      xlab = "Post-treatment time",
      ylab = "Gap",
      main = "Placebo post-treatment gaps",
      xaxt = "n",
      ...
    )
    graphics::axis(1, at = time_points, labels = time_labels)
    treated_idx <- match(x$treated_unit, x$placebo$unit_names)
    graphics::lines(time_points, x$placebo$post_gap_matrix[treated_idx, ], lwd = 2, col = "black")
    graphics::abline(h = 0, lty = 3, col = "gray70")
    return(invisible(x))
  }
  gap <- c(x$pre_residuals, x$post_gaps)
  time_idx <- c(x$pre_periods, x$post_periods)
  numeric_like_time <- is.numeric(time_idx) || inherits(time_idx, "Date") || inherits(time_idx, "POSIXt")
  if (numeric_like_time) {
    graphics::plot(time_idx, gap, type = "l", xlab = "Time", ylab = "Gap", main = "Treatment gap path", ...)
    graphics::abline(v = min(x$post_periods), lty = 2, col = "gray60")
  } else {
    x_pos <- seq_along(time_idx)
    graphics::plot(x_pos, gap, type = "l", xlab = "Time", ylab = "Gap", main = "Treatment gap path", xaxt = "n", ...)
    graphics::axis(1, at = x_pos, labels = as.character(time_idx))
    graphics::abline(v = length(x$pre_periods) + 0.5, lty = 2, col = "gray60")
  }
  graphics::abline(h = 0, lty = 3, col = "gray70")
  invisible(x)
}
