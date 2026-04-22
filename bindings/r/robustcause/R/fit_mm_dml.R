fit_mm_dml <- function(x = NULL,
                       d = NULL,
                       y = NULL,
                       outcome = NULL,
                       treatment = NULL,
                       controls = NULL,
                       data = NULL,
                       learner = c("lasso", "elastic_net", "random_forest", "hist_gradient_boosting"),
                       se_type = c("analytic", "bootstrap"),
                       bootstrap_replications = 200L,
                       bootstrap_seed = NULL,
                       progress = interactive(),
                       n_cores = 1L,
                       folds = 2L,
                       fold_type = c("random"),
                       seed = 123L,
                       mm_c = 4.685,
                       max_iter = 100L,
                       tolerance = 1e-7,
                       ci_level = 0.95) {
  learner <- match.arg(learner)
  se_type <- match.arg(se_type)
  fold_type <- match.arg(fold_type)
  cl <- match.call()

  if (!is.null(data) || !is.null(outcome) || !is.null(treatment) || !is.null(controls)) {
    return(.fit_mm_dml_data(
      outcome = outcome,
      treatment = treatment,
      controls = controls,
      data = data,
      learner = learner,
      se_type = se_type,
      bootstrap_replications = bootstrap_replications,
      bootstrap_seed = bootstrap_seed,
      progress = progress,
      n_cores = n_cores,
      folds = folds,
      fold_type = fold_type,
      seed = seed,
      mm_c = mm_c,
      max_iter = max_iter,
      tolerance = tolerance,
      ci_level = ci_level,
      original_call = cl
    ))
  }

  .fit_mm_dml_matrix(
    x = x,
    d = d,
    y = y,
    learner = learner,
    se_type = se_type,
    bootstrap_replications = bootstrap_replications,
    bootstrap_seed = bootstrap_seed,
    progress = progress,
    n_cores = n_cores,
    folds = folds,
    fold_type = fold_type,
    seed = seed,
    mm_c = mm_c,
    max_iter = max_iter,
    tolerance = tolerance,
    ci_level = ci_level,
    original_call = cl
  )
}

.fit_mm_dml_data <- function(outcome,
                             treatment,
                             controls,
                             data,
                             learner,
                             se_type,
                             bootstrap_replications,
                             bootstrap_seed,
                             progress,
                             n_cores,
                             folds,
                             fold_type,
                             seed,
                             mm_c,
                             max_iter,
                             tolerance,
                             ci_level,
                             original_call) {
  if (is.null(data) || !is.data.frame(data)) {
    stop("`data` must be a data.frame.", call. = FALSE)
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

  needed <- unique(c(outcome, treatment, controls))
  df <- stats::na.omit(data[, needed, drop = FALSE])
  x <- if (length(controls) > 0L) {
    as.matrix(df[, controls, drop = FALSE])
  } else {
    matrix(0, nrow = nrow(df), ncol = 1)
  }

  fit <- .fit_mm_dml_matrix(
    x = x,
    d = as.numeric(df[[treatment]]),
    y = as.numeric(df[[outcome]]),
    learner = learner,
    se_type = se_type,
    bootstrap_replications = bootstrap_replications,
    bootstrap_seed = bootstrap_seed,
    progress = progress,
    n_cores = n_cores,
    folds = folds,
    fold_type = fold_type,
    seed = seed,
    mm_c = mm_c,
    max_iter = max_iter,
    tolerance = tolerance,
    ci_level = ci_level,
    original_call = original_call
  )

  fit$controls <- if (length(controls) > 0L) controls else "(none)"
  fit$treatment_name <- treatment
  fit$treatment_names <- treatment
  fit
}

.fit_mm_dml_matrix <- function(x,
                               d,
                               y,
                               learner,
                               se_type,
                               bootstrap_replications,
                               bootstrap_seed,
                               progress,
                               n_cores,
                               folds,
                               fold_type,
                               seed,
                               mm_c,
                               max_iter,
                               tolerance,
                               ci_level,
                               original_call) {
  if (!identical(fold_type, "random")) {
    stop("Only `fold_type = 'random'` is currently supported in robustcause.", call. = FALSE)
  }
  if (is.null(x) || is.null(d) || is.null(y)) {
    stop("Provide either `x`, `d`, and `y`, or `data`, `outcome`, and `treatment`.", call. = FALSE)
  }

  x <- as.matrix(x)
  d <- as.numeric(d)
  y <- as.numeric(y)
  if (!is.numeric(x) || !is.numeric(d) || !is.numeric(y)) {
    stop("`x`, `d`, and `y` must be numeric.", call. = FALSE)
  }
  if (nrow(x) != length(d) || length(d) != length(y)) {
    stop("`x`, `d`, and `y` must have the same number of observations.", call. = FALSE)
  }
  if (folds < 2L) {
    stop("`folds` must be at least 2.", call. = FALSE)
  }

  fit <- .run_mm_dml_core(
    x = x,
    d = d,
    y = y,
    learner = learner,
    folds = folds,
    seed = seed,
    mm_c = mm_c,
    max_iter = max_iter,
    tolerance = tolerance,
    ci_level = ci_level
  )

  fit$backend <- "cpp"
  fit$learner_name <- learner
  fit$call <- original_call %||% match.call()
  fit$treatment_name <- deparse(substitute(d))
  fit$treatment_names <- fit$treatment_name
  fit$nobs <- nrow(x)
  fit$folds <- as.integer(folds)
  fit$fold_type <- fold_type
  fit$progress <- isTRUE(progress)
  fit$n_cores <- as.integer(n_cores)
  fit$se_type <- se_type
  fit$bootstrap_replications <- if (identical(se_type, "bootstrap")) as.integer(bootstrap_replications) else 0L
  fit$bootstrap_estimates <- NULL
  fit$controls <- colnames(x)

  if (identical(se_type, "bootstrap")) {
    boot <- .bootstrap_mm_dml_matrix(
      x = x,
      d = d,
      y = y,
      learner = learner,
      folds = folds,
      progress = progress,
      seed = seed,
      mm_c = mm_c,
      max_iter = max_iter,
      tolerance = tolerance,
      ci_level = ci_level,
      bootstrap_replications = bootstrap_replications,
      bootstrap_seed = bootstrap_seed %||% (seed + 10000L),
      n_cores = n_cores
    )
    fit$std_error <- unname(boot$std_error[["d"]])
    fit$ci_lower <- unname(boot$ci_lower[["d"]])
    fit$ci_upper <- unname(boot$ci_upper[["d"]])
    fit$bootstrap_estimates <- boot$estimates
    fit$bootstrap_intercept_se <- unname(boot$std_error[["(Intercept)"]])
  } else {
    fit$bootstrap_intercept_se <- NA_real_
  }

  class(fit) <- "robustcause_mm_dml"
  fit
}

.run_mm_dml_core <- function(x,
                             d,
                             y,
                             learner,
                             folds,
                             seed,
                             mm_c,
                             max_iter,
                             tolerance,
                             ci_level) {
  .Call(
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
}

.bootstrap_mm_dml_matrix <- function(x,
                                     d,
                                     y,
                                     learner,
                                     folds,
                                     progress,
                                     seed,
                                     mm_c,
                                     max_iter,
                                     tolerance,
                                     ci_level,
                                     bootstrap_replications,
                                     bootstrap_seed,
                                     n_cores) {
  set.seed(bootstrap_seed)
  n <- length(y)
  draws <- lapply(seq_len(bootstrap_replications), function(...) sample.int(n, size = n, replace = TRUE))
  bootstrap_fun <- function(b) {
    fit_b <- try(
      .run_mm_dml_core(
        x = x[draws[[b]], , drop = FALSE],
        d = d[draws[[b]]],
        y = y[draws[[b]]],
        learner = learner,
        folds = folds,
        seed = seed + b,
        mm_c = mm_c,
        max_iter = max_iter,
        tolerance = tolerance,
        ci_level = ci_level
      ),
      silent = TRUE
    )
    if (inherits(fit_b, "try-error")) {
      return(c(NA_real_, NA_real_))
    }
    c(fit_b$intercept_hat, fit_b$tau_hat)
  }

  worker_count <- max(1L, as.integer(n_cores %||% 1L))
  if (worker_count <= 1L || bootstrap_replications <= 1L) {
    estimates_list <- .run_with_progress(seq_len(bootstrap_replications), bootstrap_fun, progress = progress, label = "Bootstrap")
  } else {
    cl <- parallel::makePSOCKcluster(worker_count)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterCall(cl, function(paths) .libPaths(paths), .libPaths())
    parallel::clusterEvalQ(cl, {
      suppressPackageStartupMessages(library(robustcause))
      NULL
    })
    parallel::clusterExport(
      cl,
      varlist = c("x", "d", "y", "draws", "learner", "folds", "seed", "mm_c", "max_iter", "tolerance", "ci_level"),
      envir = environment()
    )

    chunk_size <- max(1L, ceiling(bootstrap_replications / max(1L, worker_count * 8L)))
    chunks <- split(seq_len(bootstrap_replications), ceiling(seq_len(bootstrap_replications) / chunk_size))
    start_time <- Sys.time()
    completed <- 0L
    pb <- NULL
    if (isTRUE(progress)) {
      pb <- utils::txtProgressBar(min = 0, max = bootstrap_replications, style = 3)
      on.exit(close(pb), add = TRUE)
    }
    estimates_list <- vector("list", bootstrap_replications)
    for (chunk in chunks) {
      chunk_results <- parallel::parLapply(cl, chunk, function(b) {
        fit_b <- try(
          robustcause:::`.run_mm_dml_core`(
            x = x[draws[[b]], , drop = FALSE],
            d = d[draws[[b]]],
            y = y[draws[[b]]],
            learner = learner,
            folds = folds,
            seed = seed + b,
            mm_c = mm_c,
            max_iter = max_iter,
            tolerance = tolerance,
            ci_level = ci_level
          ),
          silent = TRUE
        )
        if (inherits(fit_b, "try-error")) {
          return(c(NA_real_, NA_real_))
        }
        c(fit_b$intercept_hat, fit_b$tau_hat)
      })
      estimates_list[completed + seq_along(chunk_results)] <- chunk_results
      completed <- completed + length(chunk_results)
      .update_progress_bar(pb, completed, bootstrap_replications, start_time, "Bootstrap")
    }
  }

  estimates <- do.call(rbind, estimates_list)
  colnames(estimates) <- c("(Intercept)", "d")
  estimates <- estimates[stats::complete.cases(estimates), , drop = FALSE]
  if (nrow(estimates) == 0L && worker_count > 1L) {
    estimates_list <- .run_with_progress(seq_len(bootstrap_replications), bootstrap_fun, progress = FALSE, label = "Bootstrap")
    estimates <- do.call(rbind, estimates_list)
    colnames(estimates) <- c("(Intercept)", "d")
    estimates <- estimates[stats::complete.cases(estimates), , drop = FALSE]
  }
  alpha <- (1 - ci_level) / 2

  list(
    estimates = estimates,
    std_error = as.list(setNames(apply(estimates, 2, stats::sd), colnames(estimates))),
    ci_lower = as.list(setNames(apply(estimates, 2, stats::quantile, probs = alpha), colnames(estimates))),
    ci_upper = as.list(setNames(apply(estimates, 2, stats::quantile, probs = 1 - alpha), colnames(estimates)))
  )
}

.run_with_progress <- function(X, fun, progress = interactive(), label = "Bootstrap") {
  total <- length(X)
  start_time <- Sys.time()
  pb <- NULL
  if (isTRUE(progress)) {
    pb <- utils::txtProgressBar(min = 0, max = total, style = 3)
    on.exit(close(pb), add = TRUE)
  }
  results <- vector("list", total)
  for (i in seq_along(X)) {
    results[[i]] <- fun(X[[i]])
    .update_progress_bar(pb, i, total, start_time, label)
  }
  results
}

.update_progress_bar <- function(pb, completed, total, start_time, label = "Bootstrap") {
  if (is.null(pb)) {
    return(invisible(NULL))
  }
  utils::setTxtProgressBar(pb, completed)
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  eta <- if (completed > 0L && completed < total) elapsed / completed * (total - completed) else 0
  message(
    sprintf(
      "\r%s: %d/%d | elapsed %s | eta %s",
      label,
      completed,
      total,
      .format_eta_seconds(elapsed),
      .format_eta_seconds(eta)
    ),
    appendLF = completed >= total
  )
  invisible(NULL)
}

.format_eta_seconds <- function(seconds) {
  seconds <- max(0, round(seconds))
  mins <- seconds %/% 60
  secs <- seconds %% 60
  sprintf("%02d:%02d", mins, secs)
}

summary.robustcause_mm_dml <- function(object, ...) {
  residuals <- as.numeric(object$stage_residuals)
  d_res <- as.numeric(object$d_residuals)
  n <- length(residuals)
  p <- 2L
  df_residual <- max(0L, n - p)
  sigma <- if (df_residual > 0L) sqrt(sum(residuals^2) / df_residual) else NA_real_
  intercept_se <- if (identical(object$se_type, "bootstrap") && is.finite(object$bootstrap_intercept_se)) {
    object$bootstrap_intercept_se
  } else if (length(d_res) > 0L && df_residual > 0L) {
    sigma * sqrt(1 / n + mean(d_res)^2 / sum((d_res - mean(d_res))^2))
  } else {
    NA_real_
  }

  coef_est <- c(object$intercept_hat, object$tau_hat)
  coef_se <- c(intercept_se, object$std_error)
  coef_t <- coef_est / coef_se
  coef_p <- 2 * stats::pt(abs(coef_t), df = df_residual, lower.tail = FALSE)
  coef_table <- cbind(
    Estimate = coef_est,
    `Std. Error` = coef_se,
    `t value` = coef_t,
    `Pr(>|t|)` = coef_p
  )
  rownames(coef_table) <- c("(Intercept)", object$treatment_name %||% "d")

  structure(
    list(
      call = object$call,
      coefficients = coef_table,
      learner = object$learner_name,
      backend = object$backend %||% "cpp",
      iterations = object$iterations,
      converged = object$converged,
      nobs = object$nobs,
      folds = object$folds,
      fold_type = object$fold_type %||% "random",
      se_type = object$se_type %||% "analytic",
      bootstrap_replications = object$bootstrap_replications %||% 0L,
      n_cores = object$n_cores %||% 1L
    ),
    class = "summary.robustcause_mm_dml"
  )
}

print.summary.robustcause_mm_dml <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat(
    sprintf(
      "MM-DML fit | Learner: %s | Folds: %d | Backend: %s | N: %d\n\n",
      x$learner,
      x$folds,
      x$backend,
      x$nobs
    )
  )
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  stats::printCoefmat(x$coefficients, digits = digits, signif.stars = TRUE, na.print = "NA")
  cat("\n")
  cat(sprintf("Fold assignment: %s\n", x$fold_type))
  if (identical(x$se_type, "bootstrap")) {
    cat(sprintf("SE method: bootstrap | Bootstrap replications: %d | Cores: %d\n", x$bootstrap_replications, x$n_cores))
  } else {
    cat("SE method: analytic\n")
  }
  cat(sprintf("MM iterations: %s | Converged: %s\n", as.character(x$iterations), ifelse(isTRUE(x$converged), "TRUE", "FALSE")))
  invisible(x)
}

print.robustcause_mm_dml <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}
