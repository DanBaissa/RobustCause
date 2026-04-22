.prepare_mmm_inputs <- function(y = NULL,
                                unit_ids = NULL,
                                channels = NULL,
                                preclean_x = NULL,
                                controls = NULL,
                                data = NULL,
                                outcome = NULL,
                                unit = NULL,
                                channel_cols = NULL,
                                preclean_cols = NULL,
                                control_cols = NULL) {
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
    controls <- if (length(control_cols) > 0L) {
      as.matrix(mf[, control_cols, drop = FALSE])
    } else {
      matrix(numeric(), nrow = nrow(mf), ncol = 0L)
    }
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

  list(
    y = y,
    unit_ids = unit_ids,
    channels = channels,
    preclean_x = preclean_x,
    controls = controls,
    n = n,
    k = ncol(channels),
    channel_names = channel_names,
    control_names = control_names
  )
}

.recycle_mmm_arg <- function(value, k, name, choices = NULL) {
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

.mmm_robustness_defaults <- function(robustness,
                                     fit_method,
                                     increment_method,
                                     preclean,
                                     preclean_method,
                                     nonnegative,
                                     apply_hill) {
  if (is.null(robustness)) {
    return(list(
      robustness = "custom",
      fit_method = fit_method,
      increment_method = increment_method,
      preclean = preclean,
      preclean_method = preclean_method,
      nonnegative = nonnegative,
      apply_hill = apply_hill
    ))
  }

  robustness <- match.arg(robustness, c("classical", "huber", "s", "mm", "custom"))
  defaults <- switch(
    robustness,
    classical = list(fit_method = "ols", increment_method = "plain", preclean = FALSE, preclean_method = "mm", nonnegative = FALSE, apply_hill = TRUE),
    huber = list(fit_method = "huber", increment_method = "huber", preclean = FALSE, preclean_method = "mm", nonnegative = FALSE, apply_hill = TRUE),
    s = list(fit_method = "s", increment_method = "adaptive_clip", preclean = TRUE, preclean_method = "s", nonnegative = TRUE, apply_hill = TRUE),
    mm = list(fit_method = "mm", increment_method = "adaptive_clip", preclean = TRUE, preclean_method = "mm", nonnegative = TRUE, apply_hill = TRUE),
    custom = list(fit_method = fit_method, increment_method = increment_method, preclean = preclean, preclean_method = preclean_method, nonnegative = nonnegative, apply_hill = apply_hill)
  )

  list(
    robustness = robustness,
    fit_method = fit_method %||% defaults$fit_method,
    increment_method = increment_method %||% defaults$increment_method,
    preclean = preclean %||% defaults$preclean,
    preclean_method = preclean_method %||% defaults$preclean_method,
    nonnegative = nonnegative %||% defaults$nonnegative,
    apply_hill = apply_hill %||% defaults$apply_hill
  )
}

.as_channel_grid <- function(grid, k, name) {
  if (is.null(grid)) {
    return(NULL)
  }
  if (is.list(grid)) {
    if (length(grid) != k) {
      stop("`", name, "` list must have one element per channel.", call. = FALSE)
    }
    return(lapply(grid, function(x) unique(as.vector(x))))
  }
  values <- unique(as.vector(grid))
  rep(list(values), k)
}

.expand_mmm_search_grid <- function(k,
                                    rho_grid = NULL,
                                    hill_grid = NULL,
                                    max_grid_size = 2500L) {
  if (is.null(rho_grid) && is.null(hill_grid)) {
    return(NULL)
  }

  rho_grid <- rho_grid %||% rep(list(NA_real_), k)
  hill_grid <- hill_grid %||% rep(list(NA_real_), k)
  if (length(rho_grid) != k || length(hill_grid) != k) {
    stop("Search grids must align with the number of channels.", call. = FALSE)
  }

  for (j in seq_len(k)) {
    rho_grid[[j]] <- as.numeric(rho_grid[[j]])
    hill_grid[[j]] <- as.numeric(hill_grid[[j]])
    if (any(!is.finite(rho_grid[[j]])) || any(rho_grid[[j]] < 0 | rho_grid[[j]] >= 1)) {
      stop("`rho_grid` values must be finite and lie in [0, 1).", call. = FALSE)
    }
    if (any(!is.finite(hill_grid[[j]])) || any(hill_grid[[j]] < 0)) {
      stop("`hill_grid` values must be finite and nonnegative.", call. = FALSE)
    }
  }

  n_candidates <- prod(vapply(rho_grid, length, integer(1))) * prod(vapply(hill_grid, length, integer(1)))
  if (n_candidates > max_grid_size) {
    stop("MMM search grid is too large (", n_candidates, " candidates). Reduce the grids or increase `max_grid_size`.", call. = FALSE)
  }

  rho_df <- expand.grid(rho_grid, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  hill_df <- expand.grid(hill_grid, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  grid <- merge(
    cbind(.rho_key = 1L, rho_df),
    cbind(.rho_key = 1L, hill_df),
    by = ".rho_key",
    sort = FALSE
  )
  grid$.rho_key <- NULL
  names(grid) <- c(paste0("rho_", seq_len(k)), paste0("hill_", seq_len(k)))
  grid
}

.score_mmm_fit <- function(y, fit, choose_by) {
  if (!isTRUE(fit$ok)) {
    return(list(rss = Inf, mse = Inf, r2 = NA_real_, aic = Inf, bic = Inf, criterion = Inf))
  }
  fitted <- as.numeric(fit$design_matrix %*% fit$full_coef)
  resid <- y - fitted
  rss <- sum(resid^2)
  mse <- mean(resid^2)
  tss <- sum((y - mean(y))^2)
  r2 <- if (tss > 0) 1 - rss / tss else NA_real_
  k_par <- length(fit$full_coef)
  n <- length(y)
  aic <- n * log(max(rss / n, 1e-12)) + 2 * k_par
  bic <- n * log(max(rss / n, 1e-12)) + log(n) * k_par
  criterion <- switch(choose_by, mse = mse, aic = aic, bic = bic)
  list(rss = rss, mse = mse, r2 = r2, aic = aic, bic = bic, criterion = criterion)
}

.fit_mmm_fixed_impl <- function(inputs,
                                fit_method,
                                rho,
                                increment_method,
                                hill_lambda,
                                preclean,
                                preclean_method,
                                nonnegative,
                                residual_clip_multiplier,
                                clip_c,
                                adaptive_upper,
                                apply_hill,
                                huber_c,
                                huber_max_iter,
                                s_starts,
                                s_max_iter,
                                mm_max_iter,
                                mm_tukey_c,
                                seed,
                                call = NULL,
                                robustness = "custom",
                                choose_by = NULL,
                                grid_results = NULL) {
  k <- inputs$k
  fit_method <- match.arg(fit_method, c("huber", "ols", "s", "mm"))
  rho <- as.numeric(.recycle_mmm_arg(rho, k, "rho"))
  increment_method <- .recycle_mmm_arg(increment_method, k, "increment_method",
                                       c("plain", "huber", "tanh", "softsign", "adaptive_clip"))
  hill_lambda <- as.numeric(.recycle_mmm_arg(hill_lambda, k, "hill_lambda"))
  preclean <- as.logical(.recycle_mmm_arg(preclean, k, "preclean"))
  preclean_method <- .recycle_mmm_arg(preclean_method, k, "preclean_method", c("mm", "s"))
  nonnegative <- as.logical(.recycle_mmm_arg(nonnegative, k, "nonnegative"))
  residual_clip_multiplier <- as.numeric(.recycle_mmm_arg(residual_clip_multiplier, k, "residual_clip_multiplier"))
  clip_c <- as.numeric(.recycle_mmm_arg(clip_c, k, "clip_c"))
  adaptive_upper <- as.numeric(.recycle_mmm_arg(adaptive_upper, k, "adaptive_upper"))
  apply_hill <- as.logical(.recycle_mmm_arg(apply_hill, k, "apply_hill"))
  preclean_seed <- as.integer(seed + seq_len(k) - 1L)

  fit <- .Call(
    rc_r_fit_mmm,
    inputs$y,
    inputs$unit_ids,
    unname(inputs$channels),
    unname(inputs$preclean_x),
    unname(inputs$controls),
    list(
      fit_method = fit_method,
      huber_c = as.numeric(huber_c),
      huber_max_iter = as.integer(huber_max_iter),
      s_starts = as.integer(s_starts),
      s_max_iter = as.integer(s_max_iter),
      mm_max_iter = as.integer(mm_max_iter),
      mm_tukey_c = as.numeric(mm_tukey_c),
      seed = as.integer(seed),
      channel_names = as.character(inputs$channel_names),
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

  colnames(fit$cleaned_signals) <- inputs$channel_names
  colnames(fit$stocks) <- inputs$channel_names
  colnames(fit$transformed_channels) <- inputs$channel_names
  colnames(fit$design_matrix) <- c("(Intercept)", inputs$channel_names, inputs$control_names %||% if (ncol(inputs$controls) > 0L) paste0("control_", seq_len(ncol(inputs$controls))) else character())
  names(fit$channel_coef) <- inputs$channel_names
  names(fit$full_coef) <- colnames(fit$design_matrix)
  if (!is.finite(fit$regression_scale)) {
    resid <- inputs$y - as.numeric(fit$design_matrix %*% fit$full_coef)
    fit$regression_scale <- sqrt(sum(resid^2) / max(1, inputs$n - ncol(fit$design_matrix)))
  }
  fit$call <- call %||% match.call()
  fit$nobs <- inputs$n
  fit$channel_names <- inputs$channel_names
  fit$control_names <- inputs$control_names
  fit$robustness <- robustness
  fit$choose_by <- choose_by
  fit$grid_results <- grid_results
  fit$searched <- !is.null(grid_results)
  fit$channel_config <- data.frame(
    channel = inputs$channel_names,
    rho = rho,
    hill_lambda = hill_lambda,
    increment_method = increment_method,
    preclean = preclean,
    preclean_method = ifelse(preclean, preclean_method, "none"),
    apply_hill = apply_hill,
    row.names = NULL,
    check.names = FALSE
  )
  class(fit) <- "robustcause_mmm"
  fit
}

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
                    robustness = NULL,
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
                    rho_grid = NULL,
                    hill_grid = NULL,
                    choose_by = c("bic", "aic", "mse"),
                    max_grid_size = 2500L,
                    huber_c = 1.345,
                    huber_max_iter = 75L,
                    s_starts = 20L,
                    s_max_iter = 50L,
                    mm_max_iter = 50L,
                    mm_tukey_c = 4.685,
                    seed = 1234L) {
  fit_method_missing <- missing(fit_method)
  increment_method_missing <- missing(increment_method)
  preclean_missing <- missing(preclean)
  preclean_method_missing <- missing(preclean_method)
  nonnegative_missing <- missing(nonnegative)
  apply_hill_missing <- missing(apply_hill)
  choose_by <- match.arg(choose_by)
  fit_method <- match.arg(fit_method)

  inputs <- .prepare_mmm_inputs(
    y = y,
    unit_ids = unit_ids,
    channels = channels,
    preclean_x = preclean_x,
    controls = controls,
    data = data,
    outcome = outcome,
    unit = unit,
    channel_cols = channel_cols,
    preclean_cols = preclean_cols,
    control_cols = control_cols
  )

  if (is.null(robustness)) {
    defaults <- list(
      robustness = "custom",
      fit_method = fit_method,
      increment_method = increment_method,
      preclean = preclean,
      preclean_method = preclean_method,
      nonnegative = nonnegative,
      apply_hill = apply_hill
    )
  } else {
    defaults <- .mmm_robustness_defaults(
      robustness = robustness,
      fit_method = if (fit_method_missing) NULL else fit_method,
      increment_method = if (increment_method_missing) NULL else increment_method,
      preclean = if (preclean_missing) NULL else preclean,
      preclean_method = if (preclean_method_missing) NULL else preclean_method,
      nonnegative = if (nonnegative_missing) NULL else nonnegative,
      apply_hill = if (apply_hill_missing) NULL else apply_hill
    )
  }
  fit_method <- defaults$fit_method
  increment_method <- defaults$increment_method
  preclean <- defaults$preclean
  preclean_method <- defaults$preclean_method
  nonnegative <- defaults$nonnegative
  apply_hill <- defaults$apply_hill

  fixed_rho <- as.numeric(.recycle_mmm_arg(rho, inputs$k, "rho"))
  fixed_hill <- as.numeric(.recycle_mmm_arg(hill_lambda, inputs$k, "hill_lambda"))
  rho_grid <- .as_channel_grid(rho_grid, inputs$k, "rho_grid")
  hill_grid <- .as_channel_grid(hill_grid, inputs$k, "hill_grid")
  if (!is.null(rho_grid) && is.null(hill_grid)) {
    hill_grid <- lapply(seq_len(inputs$k), function(j) fixed_hill[[j]])
  }
  if (is.null(rho_grid) && !is.null(hill_grid)) {
    rho_grid <- lapply(seq_len(inputs$k), function(j) fixed_rho[[j]])
  }
  search_grid <- .expand_mmm_search_grid(inputs$k, rho_grid = rho_grid, hill_grid = hill_grid, max_grid_size = max_grid_size)

  if (is.null(search_grid)) {
    return(.fit_mmm_fixed_impl(
      inputs = inputs,
      fit_method = fit_method,
      rho = fixed_rho,
      increment_method = increment_method,
      hill_lambda = fixed_hill,
      preclean = preclean,
      preclean_method = preclean_method,
      nonnegative = nonnegative,
      residual_clip_multiplier = residual_clip_multiplier,
      clip_c = clip_c,
      adaptive_upper = adaptive_upper,
      apply_hill = apply_hill,
      huber_c = huber_c,
      huber_max_iter = huber_max_iter,
      s_starts = s_starts,
      s_max_iter = s_max_iter,
      mm_max_iter = mm_max_iter,
      mm_tukey_c = mm_tukey_c,
      seed = seed,
      call = match.call(),
      robustness = defaults$robustness,
      choose_by = choose_by
    ))
  }

  results <- vector("list", nrow(search_grid))
  best_idx <- NA_integer_
  best_criterion <- Inf
  best_fit <- NULL

  for (i in seq_len(nrow(search_grid))) {
    rho_i <- as.numeric(search_grid[i, paste0("rho_", seq_len(inputs$k)), drop = TRUE])
    hill_i <- as.numeric(search_grid[i, paste0("hill_", seq_len(inputs$k)), drop = TRUE])
    fit_i <- .fit_mmm_fixed_impl(
      inputs = inputs,
      fit_method = fit_method,
      rho = rho_i,
      increment_method = increment_method,
      hill_lambda = hill_i,
      preclean = preclean,
      preclean_method = preclean_method,
      nonnegative = nonnegative,
      residual_clip_multiplier = residual_clip_multiplier,
      clip_c = clip_c,
      adaptive_upper = adaptive_upper,
      apply_hill = apply_hill,
      huber_c = huber_c,
      huber_max_iter = huber_max_iter,
      s_starts = s_starts,
      s_max_iter = s_max_iter,
      mm_max_iter = mm_max_iter,
      mm_tukey_c = mm_tukey_c,
      seed = seed,
      robustness = defaults$robustness,
      choose_by = choose_by
    )
    score_i <- .score_mmm_fit(inputs$y, fit_i, choose_by = choose_by)
    results[[i]] <- data.frame(
      candidate = i,
      criterion = score_i$criterion,
      bic = score_i$bic,
      aic = score_i$aic,
      mse = score_i$mse,
      rss = score_i$rss,
      r2 = score_i$r2,
      row.names = NULL
    )
    for (j in seq_len(inputs$k)) {
      results[[i]][[paste0("rho_", inputs$channel_names[[j]])]] <- rho_i[[j]]
      results[[i]][[paste0("hill_", inputs$channel_names[[j]])]] <- hill_i[[j]]
    }
    if (is.finite(score_i$criterion) && score_i$criterion < best_criterion) {
      best_criterion <- score_i$criterion
      best_idx <- i
      best_fit <- fit_i
    }
  }

  if (is.na(best_idx) || is.null(best_fit)) {
    stop("No valid MMM candidates were fit.", call. = FALSE)
  }

  grid_results <- do.call(rbind, results)
  grid_results <- grid_results[order(grid_results$criterion, grid_results$bic, grid_results$aic), , drop = FALSE]
  rownames(grid_results) <- NULL
  best_fit$grid_results <- grid_results
  best_fit$searched <- TRUE
  best_fit$call <- match.call()
  best_fit
}

print.robustcause_mmm <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

summary.robustcause_mmm <- function(object, top_n = 6L, ...) {
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
  if (!is.null(object$channel_config)) {
    channel_table <- merge(channel_table, object$channel_config, by = "channel", sort = FALSE)
    channel_table <- channel_table[match(object$channel_names, channel_table$channel), , drop = FALSE]
    rownames(channel_table) <- NULL
  }
  top_grid <- NULL
  if (!is.null(object$grid_results)) {
    top_n <- max(1L, min(as.integer(top_n), nrow(object$grid_results)))
    top_grid <- utils::head(object$grid_results, top_n)
  }
  structure(
    list(
      call = object$call,
      ok = object$ok,
      nobs = object$nobs,
      robustness = object$robustness %||% "custom",
      fit_method = object$fit_method,
      searched = isTRUE(object$searched),
      choose_by = object$choose_by,
      coefficients = coefficients,
      channel_table = channel_table,
      design_columns = colnames(object$design_matrix),
      design_ncol = ncol(object$design_matrix),
      regression_scale = object$regression_scale,
      top_grid = top_grid
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
  cat("Robustness:", x$robustness, "\n")
  cat("Fit method:", x$fit_method, "\n")
  cat("Design columns:", x$design_ncol, "\n")
  cat("Regression scale:", format(signif(x$regression_scale, 6)), "\n")
  if (isTRUE(x$searched)) {
    cat("Model search:", x$choose_by, "\n")
  }
  cat("\nChannel summary:\n")
  print(x$channel_table, row.names = FALSE)
  if (!is.null(x$top_grid)) {
    cat("\nTop grid candidates:\n")
    print(x$top_grid, row.names = FALSE)
  }
  cat("\nCoefficient table:\n")
  print(x$coefficients, row.names = FALSE)
  invisible(x)
}
