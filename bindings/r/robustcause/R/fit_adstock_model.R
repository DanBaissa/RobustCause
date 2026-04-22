.hill_transform_adstock <- function(x, lambda) {
  if (lambda <= 0) {
    return(as.numeric(x))
  }
  as.numeric(x) / (1 + lambda * as.numeric(x))
}

.as_optional_matrix <- function(x, n, name) {
  if (is.null(x)) {
    return(matrix(numeric(), nrow = n, ncol = 0L))
  }
  x <- as.matrix(x)
  if (!is.numeric(x)) {
    stop("`", name, "` must be numeric.", call. = FALSE)
  }
  if (nrow(x) != n) {
    stop("`", name, "` must have the same number of rows as the outcome.", call. = FALSE)
  }
  x
}

.resolve_adstock_level <- function(robustness,
                                   preclean_method,
                                   increment_method,
                                   fit_method) {
  robustness <- match.arg(robustness, c("classical", "huber", "s", "mm"))
  defaults <- switch(
    robustness,
    classical = list(preclean_method = "none", increment_method = "plain", fit_method = "ols"),
    huber = list(preclean_method = "none", increment_method = "huber", fit_method = "huber"),
    s = list(preclean_method = "s", increment_method = "adaptive_clip", fit_method = "s"),
    mm = list(preclean_method = "mm", increment_method = "adaptive_clip", fit_method = "mm")
  )
  list(
    robustness = robustness,
    preclean_method = preclean_method %||% defaults$preclean_method,
    increment_method = increment_method %||% defaults$increment_method,
    fit_method = fit_method %||% defaults$fit_method
  )
}

.fit_final_adstock_model <- function(y,
                                     transformed_stock,
                                     controls,
                                     fit_method,
                                     huber_c,
                                     huber_max_iter,
                                     s_starts,
                                     s_max_iter,
                                     mm_max_iter,
                                     mm_tukey_c,
                                     seed) {
  controls <- .as_optional_matrix(controls, length(y), "controls")
  X <- cbind(stock = as.numeric(transformed_stock), controls)
  colnames(X) <- c("stock", colnames(controls) %||% if (ncol(controls) > 0L) paste0("control_", seq_len(ncol(controls))) else character())

  if (identical(fit_method, "ols")) {
    X_design <- cbind("(Intercept)" = 1, X)
    fit <- stats::lm.fit(x = X_design, y = y)
    resid <- as.numeric(fit$residuals)
    fitted <- as.numeric(y - resid)
    sigma2 <- sum(resid^2) / max(1, length(y) - ncol(X_design))
    XtX_inv <- tryCatch(solve(crossprod(X_design)), error = function(e) NULL)
    se <- if (is.null(XtX_inv)) rep(NA_real_, ncol(X_design)) else sqrt(pmax(diag(XtX_inv) * sigma2, 0))
    coef <- as.numeric(fit$coefficients)
    names(coef) <- colnames(X_design)
    return(list(
      fit = fit,
      class = "ols",
      coef = coef,
      fitted = fitted,
      resid = resid,
      scale = sqrt(max(sigma2, 0)),
      se = se
    ))
  }

  if (identical(fit_method, "huber")) {
    fit <- fit_rlm(
      X,
      y,
      method = "m",
      psi = "huber",
      tuning = huber_c,
      maxit = huber_max_iter,
      add_intercept = TRUE
    )
    vc <- vcov(fit, hc_type = "HC3")
    return(list(
      fit = fit,
      class = "robustcause_rlm",
      coef = fit$coef,
      fitted = fit$fitted,
      resid = fit$resid,
      scale = fit$scale,
      se = sqrt(pmax(diag(vc), 0))
    ))
  }

  if (identical(fit_method, "mm")) {
    fit <- fit_rlm(
      X,
      y,
      method = "mm",
      tuning = mm_tukey_c,
      maxit = mm_max_iter,
      add_intercept = TRUE
    )
    vc <- vcov(fit, hc_type = "HC3")
    return(list(
      fit = fit,
      class = "robustcause_rlm",
      coef = fit$coef,
      fitted = fit$fitted,
      resid = fit$resid,
      scale = fit$scale,
      se = sqrt(pmax(diag(vc), 0))
    ))
  }

  fit <- fit_s_estimator(
    X,
    y,
    n_starts = s_starts,
    max_refine = s_max_iter,
    add_intercept = TRUE,
    seed = seed
  )
  return(list(
    fit = fit,
    class = "robustcause_s",
    coef = fit$coef,
    fitted = fit$fitted,
    resid = fit$resid,
    scale = fit$scale,
    se = rep(NA_real_, length(fit$coef))
  ))
}

.score_adstock_candidate <- function(y, model_fit, choose_by) {
  resid <- as.numeric(model_fit$resid)
  rss <- sum(resid^2)
  mse <- mean(resid^2)
  tss <- sum((y - mean(y))^2)
  r2 <- if (tss > 0) 1 - rss / tss else NA_real_
  k <- length(model_fit$coef)
  n <- length(y)
  aic <- n * log(max(rss / n, 1e-12)) + 2 * k
  bic <- n * log(max(rss / n, 1e-12)) + log(n) * k
  criterion <- switch(
    choose_by,
    mse = mse,
    aic = aic,
    bic = bic
  )
  list(rss = rss, mse = mse, r2 = r2, aic = aic, bic = bic, criterion = criterion)
}

fit_adstock_model <- function(y = NULL,
                              signal = NULL,
                              unit_ids = NULL,
                              controls = NULL,
                              data = NULL,
                              outcome = NULL,
                              signal_col = NULL,
                              unit_col = NULL,
                              control_cols = NULL,
                              robustness = c("classical", "huber", "s", "mm"),
                              preclean_method = NULL,
                              increment_method = NULL,
                              fit_method = NULL,
                              rho_grid = seq(0.1, 0.9, by = 0.05),
                              hill_grid = c(0, 0.05, 0.10, 0.20),
                              choose_by = c("bic", "aic", "mse"),
                              nonnegative = FALSE,
                              residual_clip_multiplier = 2.5,
                              clip_c = 1.5,
                              adaptive_upper = 0,
                              s_starts = 20L,
                              s_max_iter = 50L,
                              mm_max_iter = 50L,
                              mm_tukey_c = 4.685,
                              huber_c = 1.345,
                              huber_max_iter = 75L,
                              seed = 1234L) {
  choose_by <- match.arg(choose_by)
  level <- .resolve_adstock_level(
    robustness = robustness,
    preclean_method = preclean_method,
    increment_method = increment_method,
    fit_method = fit_method
  )

  if (!is.null(data)) {
    if (!is.data.frame(data)) {
      stop("`data` must be a data.frame.", call. = FALSE)
    }
    if (!(is.character(outcome) && length(outcome) == 1L && outcome %in% names(data))) {
      stop("`outcome` must name a single outcome column in `data`.", call. = FALSE)
    }
    if (!(is.character(signal_col) && length(signal_col) == 1L && signal_col %in% names(data))) {
      stop("`signal_col` must name a single signal column in `data`.", call. = FALSE)
    }
    if (!(is.character(unit_col) && length(unit_col) == 1L && unit_col %in% names(data))) {
      stop("`unit_col` must name a single unit-id column in `data`.", call. = FALSE)
    }
    control_cols <- control_cols %||% setdiff(names(data), c(outcome, signal_col, unit_col))
    missing_controls <- setdiff(control_cols, names(data))
    if (length(missing_controls) > 0L) {
      stop("Missing `control_cols` in `data`: ", paste(missing_controls, collapse = ", "), call. = FALSE)
    }
    used_cols <- unique(c(outcome, signal_col, unit_col, control_cols))
    mf <- stats::na.omit(data[, used_cols, drop = FALSE])
    y <- as.numeric(mf[[outcome]])
    signal <- as.numeric(mf[[signal_col]])
    unit_ids <- mf[[unit_col]]
    controls <- if (length(control_cols) > 0L) as.matrix(mf[, control_cols, drop = FALSE]) else matrix(numeric(), nrow = nrow(mf), ncol = 0L)
    control_names <- colnames(controls)
  } else {
    if (is.null(y) || is.null(signal) || is.null(unit_ids)) {
      stop("Provide either the vector interface or the data-frame interface.", call. = FALSE)
    }
    y <- as.numeric(y)
    signal <- as.numeric(signal)
    controls <- .as_optional_matrix(controls, length(y), "controls")
    control_names <- colnames(controls)
  }

  unit_ids <- as.integer(as.factor(unit_ids))
  if (length(y) != length(signal) || length(y) != length(unit_ids)) {
    stop("`y`, `signal`, and `unit_ids` must have the same length.", call. = FALSE)
  }

  rho_grid <- sort(unique(as.numeric(rho_grid)))
  hill_grid <- sort(unique(as.numeric(hill_grid)))
  if (length(rho_grid) < 1L || any(!is.finite(rho_grid)) || any(rho_grid < 0 | rho_grid >= 1)) {
    stop("`rho_grid` must contain finite values in [0, 1).", call. = FALSE)
  }
  if (length(hill_grid) < 1L || any(!is.finite(hill_grid)) || any(hill_grid < 0)) {
    stop("`hill_grid` must contain finite nonnegative values.", call. = FALSE)
  }

  preclean <- !identical(level$preclean_method, "none")
  x_preclean <- if (ncol(controls) > 0L) cbind("(Intercept)" = 1, controls) else matrix(1, nrow = length(y), ncol = 1L, dimnames = list(NULL, "(Intercept)"))
  candidate_grid <- expand.grid(rho = rho_grid, hill_lambda = hill_grid, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

  candidate_rows <- vector("list", nrow(candidate_grid))
  best_idx <- NA_integer_
  best_criterion <- Inf
  best_build <- NULL
  best_model <- NULL

  for (i in seq_len(nrow(candidate_grid))) {
    rho <- candidate_grid$rho[[i]]
    hill_lambda <- candidate_grid$hill_lambda[[i]]

    build <- build_adstock(
      signal = signal,
      unit_ids = unit_ids,
      x_preclean = x_preclean,
      preclean = preclean,
      preclean_method = if (preclean) level$preclean_method else "mm",
      nonnegative = nonnegative,
      residual_clip_multiplier = residual_clip_multiplier,
      s_starts = s_starts,
      s_max_iter = s_max_iter,
      mm_max_iter = mm_max_iter,
      mm_tukey_c = mm_tukey_c,
      seed = seed,
      rho = rho,
      increment_method = level$increment_method,
      clip_c = clip_c,
      adaptive_upper = adaptive_upper
    )

    transformed_stock <- .hill_transform_adstock(build$stock, hill_lambda)
    model_fit <- .fit_final_adstock_model(
      y = y,
      transformed_stock = transformed_stock,
      controls = controls,
      fit_method = level$fit_method,
      huber_c = huber_c,
      huber_max_iter = huber_max_iter,
      s_starts = s_starts,
      s_max_iter = s_max_iter,
      mm_max_iter = mm_max_iter,
      mm_tukey_c = mm_tukey_c,
      seed = seed
    )
    score <- .score_adstock_candidate(y, model_fit, choose_by = choose_by)

    candidate_rows[[i]] <- data.frame(
      rho = rho,
      hill_lambda = hill_lambda,
      criterion = score$criterion,
      bic = score$bic,
      aic = score$aic,
      mse = score$mse,
      rss = score$rss,
      r2 = score$r2,
      row.names = NULL
    )

    if (is.finite(score$criterion) && score$criterion < best_criterion) {
      best_criterion <- score$criterion
      best_idx <- i
      best_build <- build
      best_model <- model_fit
    }
  }

  if (is.na(best_idx)) {
    stop("No valid adstock candidates were fit.", call. = FALSE)
  }

  grid_results <- do.call(rbind, candidate_rows)
  grid_results <- grid_results[order(grid_results$criterion, grid_results$bic, grid_results$aic), , drop = FALSE]
  rownames(grid_results) <- NULL

  top_candidate <- candidate_grid[best_idx, , drop = FALSE]
  coef <- as.numeric(best_model$coef)
  coef_names <- names(best_model$coef) %||% c("(Intercept)", "stock", control_names %||% if (ncol(controls) > 0L) paste0("control_", seq_len(ncol(controls))) else character())
  names(coef) <- coef_names

  out <- list(
    call = match.call(),
    nobs = length(y),
    robustness = level$robustness,
    preclean_method = level$preclean_method,
    increment_method = level$increment_method,
    fit_method = level$fit_method,
    choose_by = choose_by,
    best_rho = top_candidate$rho[[1]],
    best_hill_lambda = top_candidate$hill_lambda[[1]],
    signal = signal,
    cleaned_signal = best_build$cleaned_signal,
    stock = best_build$stock,
    transformed_stock = .hill_transform_adstock(best_build$stock, top_candidate$hill_lambda[[1]]),
    y = y,
    controls = controls,
    control_names = control_names,
    coefficients = coef,
    coefficient_se = best_model$se,
    final_fit = best_model$fit,
    final_fit_class = best_model$class,
    grid_results = grid_results
  )
  class(out) <- "robustcause_adstock_model"
  out
}

print.robustcause_adstock_model <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

summary.robustcause_adstock_model <- function(object, top_n = 6L, ...) {
  coef_table <- data.frame(
    term = names(object$coefficients),
    estimate = unname(object$coefficients),
    std_error = object$coefficient_se,
    row.names = NULL,
    check.names = FALSE
  )
  top_n <- max(1L, min(as.integer(top_n), nrow(object$grid_results)))
  structure(
    list(
      call = object$call,
      nobs = object$nobs,
      robustness = object$robustness,
      preclean_method = object$preclean_method,
      increment_method = object$increment_method,
      fit_method = object$fit_method,
      choose_by = object$choose_by,
      best_rho = object$best_rho,
      best_hill_lambda = object$best_hill_lambda,
      signal_summary = stats::setNames(summary(object$signal), names(summary(object$signal))),
      cleaned_signal_summary = stats::setNames(summary(object$cleaned_signal), names(summary(object$cleaned_signal))),
      stock_summary = stats::setNames(summary(object$stock), names(summary(object$stock))),
      transformed_stock_summary = stats::setNames(summary(object$transformed_stock), names(summary(object$transformed_stock))),
      coefficients = coef_table,
      top_grid = utils::head(object$grid_results, top_n)
    ),
    class = "summary.robustcause_adstock_model"
  )
}

print.summary.robustcause_adstock_model <- function(x, ...) {
  cat("RobustCause adstock model summary\n\n")
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("N:", x$nobs, "\n")
  cat("Robustness:", x$robustness, "\n")
  cat("Preclean:", x$preclean_method, "\n")
  cat("Increment:", x$increment_method, "\n")
  cat("Final fit:", x$fit_method, "\n")
  cat("Selection criterion:", x$choose_by, "\n")
  cat("Best rho:", format(signif(x$best_rho, 6)), "\n")
  cat("Best hill lambda:", format(signif(x$best_hill_lambda, 6)), "\n")
  cat("\nSignal summary:\n")
  print(x$signal_summary)
  cat("\nCleaned signal summary:\n")
  print(x$cleaned_signal_summary)
  cat("\nStock summary:\n")
  print(x$stock_summary)
  cat("\nTransformed stock summary:\n")
  print(x$transformed_stock_summary)
  cat("\nTop grid candidates:\n")
  print(x$top_grid, row.names = FALSE)
  cat("\nCoefficient table:\n")
  print(x$coefficients, row.names = FALSE)
  invisible(x)
}
