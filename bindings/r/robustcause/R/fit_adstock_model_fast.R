#' Fit a cached robust adstock response model.
#'
#' This is a faster variant of `fit_adstock_model()` for simulation studies.
#' With `cache_preclean = TRUE`, S/MM precleaning is performed once and the
#' cleaned signal is reused across the rho / Hill candidate grid.
fit_adstock_model_fast <- function(y = NULL,
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
                                   seed = 1234L,
                                   cache_preclean = TRUE,
                                   n_cores = 1L,
                                   return_preclean_diagnostics = TRUE) {
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
    controls <- if (length(control_cols) > 0L) {
      as.matrix(mf[, control_cols, drop = FALSE])
    } else {
      matrix(numeric(), nrow = nrow(mf), ncol = 0L)
    }
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

  n_cores <- as.integer(n_cores)
  if (length(n_cores) != 1L || is.na(n_cores) || n_cores < 1L) {
    stop("`n_cores` must be a positive integer.", call. = FALSE)
  }

  preclean <- !identical(level$preclean_method, "none")
  cached_preclean <- isTRUE(cache_preclean) && preclean
  x_preclean <- if (ncol(controls) > 0L) {
    cbind("(Intercept)" = 1, controls)
  } else {
    matrix(1, nrow = length(y), ncol = 1L, dimnames = list(NULL, "(Intercept)"))
  }

  if (!isTRUE(cache_preclean) && preclean && n_cores > 1L) {
    warning(
      "`cache_preclean = FALSE` with `n_cores > 1` repeats robust ",
      "precleaning inside each candidate fit.",
      call. = FALSE
    )
  }

  preclean_fit <- NULL
  candidate_signal <- signal
  candidate_preclean <- preclean
  candidate_preclean_method <- if (preclean) level$preclean_method else "mm"

  if (cached_preclean) {
    preclean_fit <- build_adstock(
      signal = signal,
      unit_ids = unit_ids,
      x_preclean = x_preclean,
      preclean = TRUE,
      preclean_method = level$preclean_method,
      nonnegative = nonnegative,
      residual_clip_multiplier = residual_clip_multiplier,
      s_starts = s_starts,
      s_max_iter = s_max_iter,
      mm_max_iter = mm_max_iter,
      mm_tukey_c = mm_tukey_c,
      seed = seed,
      rho = rho_grid[[1]],
      increment_method = "plain",
      clip_c = clip_c,
      adaptive_upper = adaptive_upper
    )
    candidate_signal <- preclean_fit$cleaned_signal
    candidate_preclean <- FALSE
    candidate_preclean_method <- "mm"
  }

  candidate_grid <- expand.grid(
    rho = rho_grid,
    hill_lambda = hill_grid,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  fit_one_candidate <- function(i) {
    rho <- candidate_grid$rho[[i]]
    hill_lambda <- candidate_grid$hill_lambda[[i]]

    tryCatch({
      build <- build_adstock(
        signal = candidate_signal,
        unit_ids = unit_ids,
        x_preclean = x_preclean,
        preclean = candidate_preclean,
        preclean_method = candidate_preclean_method,
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

      data.frame(
        candidate = i,
        rho = rho,
        hill_lambda = hill_lambda,
        criterion = score$criterion,
        bic = score$bic,
        aic = score$aic,
        mse = score$mse,
        rss = score$rss,
        r2 = score$r2,
        candidate_error = NA_character_,
        row.names = NULL
      )
    }, error = function(e) {
      data.frame(
        candidate = i,
        rho = rho,
        hill_lambda = hill_lambda,
        criterion = Inf,
        bic = Inf,
        aic = Inf,
        mse = Inf,
        rss = Inf,
        r2 = NA_real_,
        candidate_error = conditionMessage(e),
        row.names = NULL
      )
    })
  }

  candidate_rows <- .adstock_fast_lapply(
    X = seq_len(nrow(candidate_grid)),
    FUN = fit_one_candidate,
    n_cores = n_cores
  )

  grid_results_with_index <- do.call(rbind, candidate_rows)
  grid_results_with_index <- grid_results_with_index[
    order(grid_results_with_index$criterion,
          grid_results_with_index$bic,
          grid_results_with_index$aic),
    ,
    drop = FALSE
  ]
  rownames(grid_results_with_index) <- NULL

  valid <- is.finite(grid_results_with_index$criterion)
  if (!any(valid)) {
    stop("No valid adstock candidates were fit.", call. = FALSE)
  }

  best_idx <- grid_results_with_index$candidate[which(valid)[[1]]]
  top_candidate <- candidate_grid[best_idx, , drop = FALSE]

  best_build <- build_adstock(
    signal = candidate_signal,
    unit_ids = unit_ids,
    x_preclean = x_preclean,
    preclean = candidate_preclean,
    preclean_method = candidate_preclean_method,
    nonnegative = nonnegative,
    residual_clip_multiplier = residual_clip_multiplier,
    s_starts = s_starts,
    s_max_iter = s_max_iter,
    mm_max_iter = mm_max_iter,
    mm_tukey_c = mm_tukey_c,
    seed = seed,
    rho = top_candidate$rho[[1]],
    increment_method = level$increment_method,
    clip_c = clip_c,
    adaptive_upper = adaptive_upper
  )

  transformed_stock <- .hill_transform_adstock(
    best_build$stock,
    top_candidate$hill_lambda[[1]]
  )
  best_model <- .fit_final_adstock_model(
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

  cleaned_signal <- if (cached_preclean) {
    preclean_fit$cleaned_signal
  } else {
    best_build$cleaned_signal
  }

  grid_results <- grid_results_with_index[, setdiff(names(grid_results_with_index), "candidate"), drop = FALSE]
  if ("candidate_error" %in% names(grid_results) && all(is.na(grid_results$candidate_error))) {
    grid_results$candidate_error <- NULL
  }

  coef <- as.numeric(best_model$coef)
  coef_names <- names(best_model$coef) %||% c(
    "(Intercept)",
    "stock",
    control_names %||% if (ncol(controls) > 0L) paste0("control_", seq_len(ncol(controls))) else character()
  )
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
    cleaned_signal = cleaned_signal,
    stock = best_build$stock,
    transformed_stock = transformed_stock,
    y = y,
    controls = controls,
    control_names = control_names,
    coefficients = coef,
    coefficient_se = best_model$se,
    final_fit = best_model$fit,
    final_fit_class = best_model$class,
    grid_results = grid_results,
    cache_preclean = isTRUE(cache_preclean),
    cached_preclean = cached_preclean,
    n_cores = n_cores
  )

  if (isTRUE(return_preclean_diagnostics)) {
    cleaning_shift <- signal - cleaned_signal
    out$preclean_diagnostics <- data.frame(
      preclean_enabled = preclean,
      preclean_cached = cached_preclean,
      preclean_method = if (preclean) level$preclean_method else "none",
      signal_mean = mean(signal),
      cleaned_signal_mean = mean(cleaned_signal),
      cleaning_shift_mean = mean(cleaning_shift),
      cleaning_shift_abs_mean = mean(abs(cleaning_shift)),
      cleaning_shift_sd = stats::sd(cleaning_shift),
      cleaning_shift_max = max(abs(cleaning_shift)),
      row.names = NULL
    )
  }

  class(out) <- "robustcause_adstock_model"
  out
}

.adstock_fast_lapply <- function(X, FUN, n_cores = 1L) {
  n_cores <- as.integer(n_cores)
  if (n_cores <= 1L || length(X) <= 1L) {
    return(lapply(X, FUN))
  }

  n_cores <- min(n_cores, length(X))

  if (.Platform$OS.type == "unix") {
    return(parallel::mclapply(
      X = X,
      FUN = FUN,
      mc.cores = n_cores,
      mc.set.seed = TRUE
    ))
  }

  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  ok <- tryCatch({
    parallel::clusterEvalQ(cl, {
      library(robustcause)
      NULL
    })
    TRUE
  }, error = function(e) FALSE)

  if (!ok) {
    warning(
      "Could not initialize worker processes for package-level parallelism; ",
      "falling back to sequential candidate evaluation.",
      call. = FALSE
    )
    return(lapply(X, FUN))
  }

  parallel::parLapply(cl = cl, X = X, fun = FUN)
}
