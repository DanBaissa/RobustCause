simulate_sc_panel <- function(n_donors = 20L,
                              n_pre = 30L,
                              n_post = 10L,
                              n_factors = 2L,
                              true_sparsity = min(5L, n_donors),
                              factor_rho = 0.7,
                              noise_sd = 1,
                              treated_noise_sd = noise_sd,
                              treatment_effect = 1,
                              effect_shape = c("constant", "ramp", "delayed", "temporary", "none"),
                              contamination = c(
                                "none",
                                "treated_pre_spike",
                                "donor_pre_spike",
                                "treated_and_donor_pre_spike",
                                "treated_post_spike",
                                "heavy_tails",
                                "donor_pool_mismatch"
                              ),
                              contamination_rate = 0.05,
                              contamination_magnitude = 8,
                              contamination_sd = 1,
                              donor_mismatch_strength = 0,
                              heavy_tail_df = 3,
                              seed = NULL) {
  effect_shape <- match.arg(effect_shape)
  contamination <- match.arg(contamination)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  n_donors <- as.integer(n_donors)
  n_pre <- as.integer(n_pre)
  n_post <- as.integer(n_post)
  n_factors <- as.integer(n_factors)
  true_sparsity <- as.integer(min(true_sparsity, n_donors))
  n_time <- n_pre + n_post
  time_index <- seq_len(n_time)
  treatment_start <- n_pre + 1L

  if (n_donors < 2L) {
    stop("`n_donors` must be at least 2.", call. = FALSE)
  }
  if (n_pre < 2L || n_post < 1L) {
    stop("Need at least two pre-periods and one post-period.", call. = FALSE)
  }
  if (true_sparsity < 1L) {
    stop("`true_sparsity` must be at least 1.", call. = FALSE)
  }

  factors <- matrix(0, nrow = n_time, ncol = n_factors)
  factors[1, ] <- stats::rnorm(n_factors)
  if (n_time > 1L) {
    for (tt in 2:n_time) {
      factors[tt, ] <- factor_rho * factors[tt - 1L, ] + stats::rnorm(n_factors)
    }
  }

  donor_loadings <- matrix(stats::rnorm(n_donors * n_factors), nrow = n_donors, ncol = n_factors)
  donor_intercepts <- stats::rnorm(n_donors, sd = 0.5)
  donor_noise <- matrix(stats::rnorm(n_time * n_donors, sd = noise_sd), nrow = n_time, ncol = n_donors)
  donors_clean <- factors %*% t(donor_loadings)
  donors_clean <- sweep(donors_clean, 2L, donor_intercepts, "+") + donor_noise
  colnames(donors_clean) <- paste0("donor_", seq_len(n_donors))

  active_donors <- sort(sample(seq_len(n_donors), size = true_sparsity, replace = FALSE))
  raw_weights <- stats::rexp(true_sparsity)
  true_weights <- rep(0, n_donors)
  true_weights[active_donors] <- raw_weights / sum(raw_weights)
  names(true_weights) <- colnames(donors_clean)

  treated_counterfactual <- as.numeric(donors_clean %*% true_weights) +
    stats::rnorm(n_time, sd = treated_noise_sd)

  if (identical(contamination, "donor_pool_mismatch") || donor_mismatch_strength > 0) {
    extra_factor <- as.numeric(stats::filter(stats::rnorm(n_time), filter = factor_rho, method = "recursive"))
    extra_factor[!is.finite(extra_factor)] <- 0
    treated_counterfactual <- treated_counterfactual + donor_mismatch_strength * extra_factor
  }

  post_effect <- switch(
    effect_shape,
    constant = rep(treatment_effect, n_post),
    ramp = seq(0, treatment_effect, length.out = n_post),
    delayed = c(rep(0, floor(n_post / 2)), rep(treatment_effect, n_post - floor(n_post / 2))),
    temporary = c(rep(treatment_effect, ceiling(n_post / 2)), rep(0, n_post - ceiling(n_post / 2))),
    none = rep(0, n_post)
  )
  treatment_path <- c(rep(0, n_pre), post_effect)

  treated_observed <- treated_counterfactual + treatment_path
  outcomes_clean <- cbind(treated = treated_observed, donors_clean)
  outcomes_observed <- outcomes_clean

  contamination_info <- data.frame(
    row = integer(),
    col = integer(),
    unit = character(),
    period = character(),
    shock = numeric(),
    stringsAsFactors = FALSE
  )

  add_shock <- function(rows, cols) {
    if (length(rows) == 0L || length(cols) == 0L) {
      return(invisible(NULL))
    }
    grid <- expand.grid(row = rows, col = cols, KEEP.OUT.ATTRS = FALSE)
    n_target <- max(1L, floor(contamination_rate * nrow(grid)))
    n_target <- min(n_target, nrow(grid))
    picked <- grid[sample(seq_len(nrow(grid)), size = n_target, replace = FALSE), , drop = FALSE]
    shock <- sample(c(-1, 1), nrow(picked), replace = TRUE) *
      (contamination_magnitude + stats::rnorm(nrow(picked), sd = contamination_sd))
    for (ii in seq_len(nrow(picked))) {
      outcomes_observed[picked$row[ii], picked$col[ii]] <<- outcomes_observed[picked$row[ii], picked$col[ii]] + shock[ii]
    }
    new_info <- data.frame(
      row = picked$row,
      col = picked$col,
      unit = colnames(outcomes_observed)[picked$col],
      period = ifelse(picked$row < treatment_start, "pre", "post"),
      shock = shock,
      stringsAsFactors = FALSE
    )
    contamination_info <<- rbind(contamination_info, new_info)
    invisible(NULL)
  }

  pre_rows <- seq_len(n_pre)
  post_rows <- seq.int(treatment_start, n_time)
  treated_col <- 1L
  donor_cols <- seq.int(2L, n_donors + 1L)

  if (identical(contamination, "treated_pre_spike")) {
    add_shock(pre_rows, treated_col)
  } else if (identical(contamination, "donor_pre_spike")) {
    add_shock(pre_rows, donor_cols)
  } else if (identical(contamination, "treated_and_donor_pre_spike")) {
    add_shock(pre_rows, treated_col)
    add_shock(pre_rows, donor_cols)
  } else if (identical(contamination, "treated_post_spike")) {
    add_shock(post_rows, treated_col)
  } else if (identical(contamination, "heavy_tails")) {
    heavy_noise <- matrix(stats::rt(n_time * (n_donors + 1L), df = heavy_tail_df), nrow = n_time)
    heavy_noise <- heavy_noise * noise_sd
    outcomes_observed <- outcomes_observed + heavy_noise
  }

  panel <- data.frame(
    unit = rep(colnames(outcomes_observed), each = n_time),
    time = rep(time_index, times = ncol(outcomes_observed)),
    outcome = as.numeric(outcomes_observed),
    stringsAsFactors = FALSE
  )

  structure(
    list(
      panel = panel,
      outcomes = outcomes_observed,
      outcomes_clean = outcomes_clean,
      donor_outcomes_clean = donors_clean,
      contamination_info = contamination_info,
      truth = list(
        donor_names = colnames(donors_clean),
        true_weights = true_weights,
        active_donors = names(true_weights)[true_weights > 0],
        treated_counterfactual = treated_counterfactual,
        treatment_path = treatment_path,
        treatment_effect_post = post_effect,
        treatment_start = treatment_start,
        n_pre = n_pre,
        n_post = n_post
      ),
      design = list(
        n_donors = n_donors,
        n_pre = n_pre,
        n_post = n_post,
        n_factors = n_factors,
        true_sparsity = true_sparsity,
        factor_rho = factor_rho,
        noise_sd = noise_sd,
        treated_noise_sd = treated_noise_sd,
        treatment_effect = treatment_effect,
        effect_shape = effect_shape,
        contamination = contamination,
        contamination_rate = contamination_rate,
        contamination_magnitude = contamination_magnitude,
        donor_mismatch_strength = donor_mismatch_strength,
        heavy_tail_df = heavy_tail_df
      )
    ),
    class = "robustcause_sc_simulation"
  )
}

run_sc_benchmark_once <- function(seed,
                                  methods = c("standard", "mm"),
                                  run_placebos = FALSE,
                                  return_fits = FALSE,
                                  maxit = 500L,
                                  tol = 1e-8,
                                  ridge = 1e-10,
                                  startup_maxit = 75L,
                                  subproblem_maxit = 500L,
                                  subproblem_tol = 1e-8,
                                  l1_smoothing = 1e-6,
                                  tukey_c = 4.685,
                                  min_scale = 1e-8,
                                  min_time_weight = 1e-8,
                                  ...) {
  sim <- simulate_sc_panel(seed = seed, ...)
  truth <- sim$truth
  donor_names <- truth$donor_names
  true_weights <- truth$true_weights
  true_post_effect <- truth$treatment_effect_post

  rows <- lapply(methods, function(method) {
    fit <- tryCatch(
      fit_sc(
        outcomes = sim$outcomes,
        treated_unit = "treated",
        treatment_start = truth$treatment_start,
        method = method,
        run_placebos = run_placebos,
        maxit = maxit,
        tol = tol,
        ridge = ridge,
        startup_maxit = startup_maxit,
        subproblem_maxit = subproblem_maxit,
        subproblem_tol = subproblem_tol,
        l1_smoothing = l1_smoothing,
        tukey_c = tukey_c,
        min_scale = min_scale,
        min_time_weight = min_time_weight
      ),
      error = function(e) e
    )

    base <- data.frame(
      seed = seed,
      method = method,
      status = if (inherits(fit, "error")) "error" else "ok",
      error_message = if (inherits(fit, "error")) conditionMessage(fit) else NA_character_,
      stringsAsFactors = FALSE
    )

    for (nm in names(sim$design)) {
      base[[nm]] <- sim$design[[nm]]
    }

    if (inherits(fit, "error")) {
      base$gap_rmse <- NA_real_
      base$gap_bias <- NA_real_
      base$gap_mae <- NA_real_
      base$avg_estimated_effect <- NA_real_
      base$avg_true_effect <- mean(true_post_effect)
      base$pre_rmspe <- NA_real_
      base$weight_rmse <- NA_real_
      base$weight_l1 <- NA_real_
      base$max_abs_weight_error <- NA_real_
      base$converged <- NA
      base$iterations <- NA_integer_
      base$effective_donors <- NA_real_
      base$max_weight <- NA_real_
      base$placebo_p_value_rmspe_ratio <- NA_real_
      base$placebo_p_value_avg_gap <- NA_real_
      base$downweighted_periods_n <- NA_integer_
      base$robust_time_weight_mean <- NA_real_
      return(base)
    }

    estimated_gap <- as.numeric(fit$post_gaps)
    gap_error <- estimated_gap - true_post_effect
    est_weights <- rep(0, length(donor_names))
    names(est_weights) <- donor_names
    est_weights[names(fit$weights)] <- as.numeric(fit$weights)
    weight_error <- est_weights - true_weights

    base$gap_rmse <- sqrt(mean(gap_error^2))
    base$gap_bias <- mean(gap_error)
    base$gap_mae <- mean(abs(gap_error))
    base$avg_estimated_effect <- mean(estimated_gap)
    base$avg_true_effect <- mean(true_post_effect)
    base$pre_rmspe <- fit$pre_rmspe
    base$weight_rmse <- sqrt(mean(weight_error^2))
    base$weight_l1 <- sum(abs(weight_error))
    base$max_abs_weight_error <- max(abs(weight_error))
    base$converged <- isTRUE(fit$converged)
    base$iterations <- fit$iterations
    base$effective_donors <- fit$effective_donors
    base$max_weight <- fit$max_weight
    base$placebo_p_value_rmspe_ratio <- if (!is.null(fit$inference)) fit$inference$placebo_p_value_rmspe_ratio else NA_real_
    base$placebo_p_value_avg_gap <- if (!is.null(fit$inference)) fit$inference$placebo_p_value_avg_gap else NA_real_
    base$downweighted_periods_n <- if (!is.null(fit$downweighted_periods)) length(fit$downweighted_periods) else NA_integer_
    base$robust_time_weight_mean <- if (!is.null(fit$robust_time_weights)) mean(fit$robust_time_weights) else NA_real_

    if (isTRUE(return_fits)) {
      base$fit <- I(list(fit))
      base$simulation <- I(list(sim))
    }

    base
  })

  do.call(rbind, rows)
}

make_sc_benchmark_design <- function(contamination = c(
                                 "none",
                                 "treated_pre_spike",
                                 "donor_pre_spike",
                                 "treated_and_donor_pre_spike",
                                 "heavy_tails",
                                 "donor_pool_mismatch"
                               ),
                               n_donors = c(10L, 25L),
                               n_pre = c(15L, 30L, 60L),
                               n_post = 10L,
                               n_factors = 2L,
                               true_sparsity = 5L,
                               treatment_effect = c(0, 1),
                               effect_shape = c("constant", "ramp"),
                               contamination_rate = c(0, 0.05, 0.10),
                               contamination_magnitude = c(4, 8),
                               donor_mismatch_strength = c(0, 1),
                               noise_sd = 1) {
  design <- expand.grid(
    contamination = contamination,
    n_donors = n_donors,
    n_pre = n_pre,
    n_post = n_post,
    n_factors = n_factors,
    true_sparsity = true_sparsity,
    treatment_effect = treatment_effect,
    effect_shape = effect_shape,
    contamination_rate = contamination_rate,
    contamination_magnitude = contamination_magnitude,
    donor_mismatch_strength = donor_mismatch_strength,
    noise_sd = noise_sd,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  keep_none <- design$contamination == "none" &
    design$contamination_rate == 0 &
    design$contamination_magnitude == contamination_magnitude[[1]] &
    design$donor_mismatch_strength == 0
  keep_mismatch <- design$contamination == "donor_pool_mismatch" &
    design$contamination_rate == 0 &
    design$contamination_magnitude == contamination_magnitude[[1]] &
    design$donor_mismatch_strength > 0
  keep_other <- !(design$contamination %in% c("none", "donor_pool_mismatch")) &
    design$contamination_rate > 0 &
    design$donor_mismatch_strength == 0

  design <- design[keep_none | keep_mismatch | keep_other, , drop = FALSE]
  rownames(design) <- NULL
  design$design_id <- seq_len(nrow(design))
  design
}

run_sc_benchmark_grid <- function(design,
                                  reps = 100L,
                                  methods = c("standard", "mm"),
                                  run_placebos = FALSE,
                                  n_cores = 1L,
                                  seed = 100000L,
                                  out_dir = NULL,
                                  chunk_id = NULL,
                                  return_fits = FALSE) {
  if (!is.data.frame(design)) {
    stop("`design` must be a data.frame, usually from `make_sc_benchmark_design()`.", call. = FALSE)
  }
  reps <- as.integer(reps)
  n_cores <- as.integer(n_cores)
  if (reps < 1L || n_cores < 1L) {
    stop("`reps` and `n_cores` must be positive integers.", call. = FALSE)
  }

  if (is.null(design$design_id)) {
    design$design_id <- seq_len(nrow(design))
  }

  if (!is.null(out_dir)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    label <- if (is.null(chunk_id)) "all" else as.character(chunk_id)
    out_file <- file.path(out_dir, paste0("sc_benchmark_", label, ".rds"))
    if (file.exists(out_file)) {
      return(readRDS(out_file))
    }
  } else {
    out_file <- NULL
  }

  task_grid <- expand.grid(
    design_row = seq_len(nrow(design)),
    replication = seq_len(reps),
    KEEP.OUT.ATTRS = FALSE
  )
  task_grid$seed <- seed + seq_len(nrow(task_grid))

  simulate_arg_names <- setdiff(names(formals(simulate_sc_panel)), "seed")

  run_task <- function(ii) {
    task <- task_grid[ii, , drop = FALSE]
    row <- design[task$design_row, , drop = FALSE]
    args <- as.list(row[intersect(names(row), simulate_arg_names)])
    args <- lapply(args, function(x) x[[1]])

    ans <- do.call(
      run_sc_benchmark_once,
      c(
        list(
          seed = task$seed,
          methods = methods,
          run_placebos = run_placebos,
          return_fits = return_fits
        ),
        args
      )
    )
    ans$design_id <- row$design_id[[1]]
    ans$replication <- task$replication
    ans
  }

  pieces <- .sc_benchmark_lapply(seq_len(nrow(task_grid)), run_task, n_cores = n_cores)
  results <- do.call(rbind, pieces)
  rownames(results) <- NULL

  if (!is.null(out_file)) {
    saveRDS(results, out_file)
  }

  results
}

summarize_sc_benchmark <- function(results,
                                   group_vars = c("contamination", "n_donors", "n_pre", "treatment_effect", "method")) {
  if (!is.data.frame(results)) {
    stop("`results` must be a data.frame returned by `run_sc_benchmark_grid()`.", call. = FALSE)
  }
  group_vars <- intersect(group_vars, names(results))
  if (length(group_vars) == 0L) {
    group_vars <- "method"
  }

  split_key <- interaction(results[group_vars], drop = TRUE, sep = " | ")
  groups <- split(results, split_key)

  summaries <- lapply(groups, function(x) {
    key <- x[1L, group_vars, drop = FALSE]
    data.frame(
      key,
      n = nrow(x),
      ok_rate = mean(x$status == "ok"),
      gap_rmse = mean(x$gap_rmse, na.rm = TRUE),
      gap_bias = mean(x$gap_bias, na.rm = TRUE),
      gap_mae = mean(x$gap_mae, na.rm = TRUE),
      weight_rmse = mean(x$weight_rmse, na.rm = TRUE),
      pre_rmspe = mean(x$pre_rmspe, na.rm = TRUE),
      avg_effect_hat = mean(x$avg_estimated_effect, na.rm = TRUE),
      avg_effect_true = mean(x$avg_true_effect, na.rm = TRUE),
      convergence_rate = mean(x$converged, na.rm = TRUE),
      mean_iterations = mean(x$iterations, na.rm = TRUE),
      mean_effective_donors = mean(x$effective_donors, na.rm = TRUE),
      mean_downweighted_periods = mean(x$downweighted_periods_n, na.rm = TRUE),
      row.names = NULL,
      check.names = FALSE
    )
  })

  out <- do.call(rbind, summaries)
  rownames(out) <- NULL
  out
}

.sc_benchmark_lapply <- function(X, FUN, n_cores = 1L) {
  n_cores <- as.integer(n_cores)
  if (n_cores <= 1L || length(X) <= 1L) {
    return(lapply(X, FUN))
  }
  n_cores <- min(n_cores, length(X))

  if (.Platform$OS.type == "unix") {
    return(parallel::mclapply(X, FUN, mc.cores = n_cores, mc.set.seed = TRUE))
  }

  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterEvalQ(cl, library(robustcause))
  parallel::parLapply(cl, X, FUN)
}
