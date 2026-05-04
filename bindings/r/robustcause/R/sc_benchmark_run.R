run_sc_benchmark_once <- function(seed,
                                  methods = c("standard", "mm"),
                                  run_placebos = FALSE,
                                  return_fits = FALSE,
                                  predictor_lambda = 1,
                                  predictor_scale = c("mad", "sd", "none"),
                                  robust_donors = FALSE,
                                  donor_penalty_lambda = 1,
                                  donor_tukey_c = 4.685,
                                  min_donor_weight = 1e-4,
                                  max_donor_penalty = 1e4,
                                  ...) {
  predictor_scale <- match.arg(predictor_scale)
  sim <- simulate_sc_panel(seed = seed, ...)
  truth <- sim$truth
  true_weights <- truth$true_weights
  true_post_effect <- truth$treatment_effect_post

  rows <- lapply(methods, function(method) {
    fit <- tryCatch(
      fit_sc(
        outcomes = sim$outcomes,
        treated_unit = "treated",
        treatment_start = truth$treatment_start,
        predictors = sim$predictors,
        predictor_lambda = predictor_lambda,
        predictor_scale = predictor_scale,
        method = method,
        robust_donors = robust_donors,
        donor_penalty_lambda = donor_penalty_lambda,
        donor_tukey_c = donor_tukey_c,
        min_donor_weight = min_donor_weight,
        max_donor_penalty = max_donor_penalty,
        run_placebos = run_placebos
      ),
      error = function(e) e
    )

    base <- data.frame(
      seed = seed,
      method = method,
      robust_donors = isTRUE(robust_donors),
      donor_penalty_lambda = donor_penalty_lambda,
      status = if (inherits(fit, "error")) "error" else "ok",
      error_message = if (inherits(fit, "error")) conditionMessage(fit) else NA_character_,
      stringsAsFactors = FALSE
    )

    for (nm in names(sim$design)) base[[nm]] <- sim$design[[nm]]

    if (inherits(fit, "error")) {
      base$gap_rmse <- NA_real_
      base$gap_bias <- NA_real_
      base$gap_mae <- NA_real_
      base$avg_estimated_effect <- NA_real_
      base$avg_true_effect <- mean(true_post_effect)
      base$pre_rmspe <- NA_real_
      base$predictor_rmse <- NA_real_
      base$mean_robust_predictor_weight <- NA_real_
      base$min_robust_predictor_weight <- NA_real_
      base$mean_robust_donor_weight <- NA_real_
      base$min_robust_donor_weight <- NA_real_
      base$max_donor_penalty_fit <- NA_real_
      base$weight_rmse <- NA_real_
      base$weight_l1 <- NA_real_
      base$placebo_p_value_rmspe_ratio <- NA_real_
      base$placebo_p_value_avg_gap <- NA_real_
      return(base)
    }

    est_weights <- rep(0, length(true_weights))
    names(est_weights) <- names(true_weights)
    est_weights[names(fit$weights)] <- as.numeric(fit$weights)
    weight_error <- est_weights - true_weights
    gap_error <- as.numeric(fit$post_gaps) - true_post_effect

    base$gap_rmse <- sqrt(mean(gap_error^2))
    base$gap_bias <- mean(gap_error)
    base$gap_mae <- mean(abs(gap_error))
    base$avg_estimated_effect <- mean(fit$post_gaps)
    base$avg_true_effect <- mean(true_post_effect)
    base$pre_rmspe <- fit$pre_rmspe
    base$predictor_rmse <- fit$predictor_rmse
    base$mean_robust_predictor_weight <- if (!is.null(fit$robust_predictor_weights) && length(fit$robust_predictor_weights) > 0L) mean(fit$robust_predictor_weights) else NA_real_
    base$min_robust_predictor_weight <- if (!is.null(fit$robust_predictor_weights) && length(fit$robust_predictor_weights) > 0L) min(fit$robust_predictor_weights) else NA_real_
    base$mean_robust_donor_weight <- if (!is.null(fit$robust_donor_weights) && length(fit$robust_donor_weights) > 0L) mean(fit$robust_donor_weights) else NA_real_
    base$min_robust_donor_weight <- if (!is.null(fit$robust_donor_weights) && length(fit$robust_donor_weights) > 0L) min(fit$robust_donor_weights) else NA_real_
    base$max_donor_penalty_fit <- if (!is.null(fit$donor_penalties) && length(fit$donor_penalties) > 0L) max(fit$donor_penalties) else NA_real_
    base$weight_rmse <- sqrt(mean(weight_error^2))
    base$weight_l1 <- sum(abs(weight_error))
    base$placebo_p_value_rmspe_ratio <- if (!is.null(fit$inference)) fit$inference$placebo_p_value_rmspe_ratio else NA_real_
    base$placebo_p_value_avg_gap <- if (!is.null(fit$inference)) fit$inference$placebo_p_value_avg_gap else NA_real_

    if (isTRUE(return_fits)) {
      base$fit <- I(list(fit))
      base$simulation <- I(list(sim))
    }
    base
  })

  do.call(rbind, rows)
}

run_sc_benchmark_grid <- function(design,
                                  reps = 100L,
                                  methods = c("standard", "mm"),
                                  run_placebos = FALSE,
                                  n_cores = 1L,
                                  seed = 100000L,
                                  out_dir = NULL,
                                  chunk_id = NULL,
                                  return_fits = FALSE,
                                  predictor_lambda = 1,
                                  predictor_scale = c("mad", "sd", "none"),
                                  robust_donors = FALSE,
                                  donor_penalty_lambda = 1,
                                  donor_tukey_c = 4.685,
                                  min_donor_weight = 1e-4,
                                  max_donor_penalty = 1e4) {
  predictor_scale <- match.arg(predictor_scale)
  if (!is.data.frame(design)) stop("`design` must be a data.frame.", call. = FALSE)
  if (is.null(design$design_id)) design$design_id <- seq_len(nrow(design))

  if (!is.null(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    label <- if (is.null(chunk_id)) "all" else as.character(chunk_id)
    out_file <- file.path(out_dir, paste0("sc_benchmark_", label, ".rds"))
    if (file.exists(out_file)) return(readRDS(out_file))
  } else {
    out_file <- NULL
  }

  task_grid <- expand.grid(design_row = seq_len(nrow(design)), replication = seq_len(as.integer(reps)), KEEP.OUT.ATTRS = FALSE)
  task_grid$seed <- seed + seq_len(nrow(task_grid))
  sim_args <- setdiff(names(formals(simulate_sc_panel)), "seed")

  run_task <- function(ii) {
    task <- task_grid[ii, , drop = FALSE]
    row <- design[task$design_row, , drop = FALSE]
    args <- as.list(row[intersect(names(row), sim_args)])
    args <- lapply(args, function(x) x[[1]])
    ans <- do.call(
      run_sc_benchmark_once,
      c(
        list(
          seed = task$seed,
          methods = methods,
          run_placebos = run_placebos,
          return_fits = return_fits,
          predictor_lambda = predictor_lambda,
          predictor_scale = predictor_scale,
          robust_donors = robust_donors,
          donor_penalty_lambda = donor_penalty_lambda,
          donor_tukey_c = donor_tukey_c,
          min_donor_weight = min_donor_weight,
          max_donor_penalty = max_donor_penalty
        ),
        args
      )
    )
    ans$design_id <- row$design_id[[1]]
    ans$replication <- task$replication
    ans
  }

  pieces <- sc_parallel_lapply(seq_len(nrow(task_grid)), run_task, n_cores = n_cores)
  results <- do.call(rbind, pieces)
  rownames(results) <- NULL
  if (!is.null(out_file)) saveRDS(results, out_file)
  results
}

sc_parallel_lapply <- function(X, FUN, n_cores = 1L) {
  n_cores <- as.integer(n_cores)
  if (n_cores <= 1L || length(X) <= 1L) return(lapply(X, FUN))
  n_cores <- min(n_cores, length(X))
  if (.Platform$OS.type == "unix") return(parallel::mclapply(X, FUN, mc.cores = n_cores, mc.set.seed = TRUE))
  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterEvalQ(cl, library(robustcause))
  parallel::parLapply(cl, X, FUN)
}
