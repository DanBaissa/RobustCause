make_custom_learner <- function(model_y, model_t = model_y, label = "custom") {
  structure(
    list(
      model_y = model_y,
      model_t = model_t,
      label = label
    ),
    class = "robustcause_custom_learner"
  )
}

fit_mm_dml <- function(x = NULL,
                       d = NULL,
                       y = NULL,
                       outcome = NULL,
                       treatment = NULL,
                       controls = NULL,
                       second_stage_controls = NULL,
                       data = NULL,
                       learner = c("lasso", "elastic_net", "random_forest", "hist_gradient_boosting"),
                       learner_params = NULL,
                       se_type = c("analytic", "bootstrap"),
                       bootstrap_replications = 200L,
                       bootstrap_seed = NULL,
                       progress = interactive(),
                       n_cores = 1L,
                       folds = 2L,
                       fold_type = c("random", "spatial_block"),
                       block_id = NULL,
                       spatial_k = NULL,
                       spatial_by = NULL,
                       seed = 123L,
                       mm_c = 4.685,
                       max_iter = 100L,
                       tolerance = 1e-7,
                       ci_level = 0.95) {
  se_type <- match.arg(se_type)
  fold_type <- match.arg(fold_type)
  learner <- .normalize_mm_dml_learner(learner)
  learner_params <- .normalize_mm_dml_learner_params(learner, learner_params)
  cl <- match.call()

  if (!is.null(data) || !is.null(outcome) || !is.null(treatment) || !is.null(controls) || !is.null(second_stage_controls)) {
    return(.fit_mm_dml_data(
      outcome = outcome,
      treatment = treatment,
      controls = controls,
      second_stage_controls = second_stage_controls,
      data = data,
      learner = learner,
      learner_params = learner_params,
      se_type = se_type,
      bootstrap_replications = bootstrap_replications,
      bootstrap_seed = bootstrap_seed,
      progress = progress,
      n_cores = n_cores,
      folds = folds,
      fold_type = fold_type,
      block_id = block_id,
      spatial_k = spatial_k,
      spatial_by = spatial_by,
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
    learner_params = learner_params,
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

.fit_mm_dml_matrix <- function(x,
                               d,
                               y,
                               learner,
                               learner_params,
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
  if (is.null(x) || is.null(d) || is.null(y)) {
    stop("Provide either `x`, `d`, and `y`, or `data`, `outcome`, and `treatment`.", call. = FALSE)
  }

  x <- as.matrix(x)
  d <- as.matrix(d)
  y <- as.numeric(y)
  if (!is.numeric(x) || !is.numeric(d) || !is.numeric(y)) {
    stop("`x`, `d`, and `y` must be numeric.", call. = FALSE)
  }
  if (nrow(x) != nrow(d) || nrow(d) != length(y)) {
    stop("`x`, `d`, and `y` must have the same number of observations.", call. = FALSE)
  }
  if (folds < 2L) {
    stop("`folds` must be at least 2.", call. = FALSE)
  }
  if (identical(fold_type, "spatial_block")) {
    stop("`fold_type = 'spatial_block'` is currently supported only for the data-frame interface.", call. = FALSE)
  }

  if (.supports_native_mm_dml(learner, d, fold_type)) {
    fit <- .run_native_mm_dml_matrix(
      x = x,
      d = as.numeric(d[, 1]),
      y = y,
      learner = learner,
      learner_params = learner_params,
      se_type = se_type,
      bootstrap_replications = bootstrap_replications,
      bootstrap_seed = bootstrap_seed,
      progress = progress,
      n_cores = n_cores,
      folds = folds,
      seed = seed,
      mm_c = mm_c,
      max_iter = max_iter,
      tolerance = tolerance,
      ci_level = ci_level,
      original_call = original_call
    )
    fit$controls <- colnames(x)
    fit$treatment_name <- colnames(d)[1] %||% "d"
    fit$treatment_names <- fit$treatment_name
    return(fit)
  }

  result <- .run_mm_dml_joint_stage(
    x = x,
    d = d,
    y = y,
    learner = learner,
    learner_params = learner_params,
    folds = folds,
    fold_type = fold_type,
    block_id = NULL,
    second_stage_controls = NULL,
    seed = seed,
    mm_c = mm_c,
    max_iter = max_iter,
    tolerance = tolerance,
    ci_level = ci_level,
    call = original_call
  )

  if (identical(se_type, "bootstrap")) {
    boot <- .bootstrap_mm_dml_joint(
      x = x,
      d = d,
      y = y,
      learner = learner,
      learner_params = learner_params,
      folds = folds,
      fold_type = fold_type,
      block_id = NULL,
      second_stage_controls = NULL,
      seed = seed,
      mm_c = mm_c,
      max_iter = max_iter,
      tolerance = tolerance,
      ci_level = ci_level,
      bootstrap_replications = bootstrap_replications,
      bootstrap_seed = bootstrap_seed %||% (seed + 10000L),
      progress = progress,
      n_cores = n_cores
    )
    result <- .apply_bootstrap_to_joint_result(result, boot, primary_term = result$treatment_names[[1]])
  } else {
    result$bootstrap_estimates <- NULL
    result$bootstrap_replications <- 0L
    result$bootstrap_ci <- NULL
  }

  result$se_type <- se_type
  result$fold_type <- fold_type
  result$n_cores <- as.integer(n_cores)
  result$progress <- isTRUE(progress)
  result$controls <- colnames(x)
  class(result) <- "robustcause_mm_dml"
  result
}

.fit_mm_dml_data <- function(outcome,
                             treatment,
                             controls,
                             second_stage_controls,
                             data,
                             learner,
                             learner_params,
                             se_type,
                             bootstrap_replications,
                             bootstrap_seed,
                             progress,
                             n_cores,
                             folds,
                             fold_type,
                             block_id,
                             spatial_k,
                             spatial_by,
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
  if (is.null(treatment) || !all(treatment %in% names(data))) {
    stop("`treatment` must contain one or more column names in `data`.", call. = FALSE)
  }
  controls <- controls %||% setdiff(names(data), c(outcome, treatment))
  second_stage_controls <- second_stage_controls %||% character(0)

  block_column <- .normalize_block_id_argument(block_id, data)
  if (identical(fold_type, "spatial_block") && is.null(block_column)) {
    stop("`block_id` must name a column in `data` when `fold_type = 'spatial_block'`.", call. = FALSE)
  }

  needed <- unique(c(outcome, treatment, controls, second_stage_controls, block_column))
  missing_cols <- setdiff(needed, names(data))
  if (length(missing_cols) > 0L) {
    stop("Missing columns in `data`: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  df <- stats::na.omit(data[, needed, drop = FALSE])
  n <- nrow(df)
  if (n < max(10L, folds + 2L)) {
    stop("Not enough complete observations after dropping missing values.", call. = FALSE)
  }

  y_vec <- as.numeric(df[[outcome]])
  d_mat <- as.matrix(df[, treatment, drop = FALSE])
  x_mat <- if (length(controls) > 0L) {
    stats::model.matrix(stats::reformulate(controls), data = df)[, -1, drop = FALSE]
  } else {
    matrix(0, nrow = n, ncol = 0)
  }
  stage_controls_mat <- if (length(second_stage_controls) > 0L) {
    stats::model.matrix(stats::reformulate(second_stage_controls), data = df)[, -1, drop = FALSE]
  } else {
    NULL
  }

  if (.supports_native_mm_dml(learner, d_mat, fold_type) && length(treatment) == 1L && is.null(stage_controls_mat)) {
    fit <- .run_native_mm_dml_matrix(
      x = if (ncol(x_mat) == 0L) matrix(0, nrow = n, ncol = 1L) else x_mat,
      d = as.numeric(d_mat[, 1]),
      y = y_vec,
      learner = learner,
      learner_params = learner_params,
      se_type = se_type,
      bootstrap_replications = bootstrap_replications,
      bootstrap_seed = bootstrap_seed,
      progress = progress,
      n_cores = n_cores,
      folds = folds,
      seed = seed,
      mm_c = mm_c,
      max_iter = max_iter,
      tolerance = tolerance,
      ci_level = ci_level,
      original_call = original_call
    )
    fit$call <- original_call
    fit$controls <- if (length(controls) > 0L) controls else "(none)"
    fit$treatment_name <- treatment[[1]]
    fit$treatment_names <- treatment
    class(fit) <- "robustcause_mm_dml"
    return(fit)
  }

  result <- .run_mm_dml_joint_stage(
    x = if (ncol(x_mat) == 0L) matrix(0, nrow = n, ncol = 1L) else x_mat,
    d = d_mat,
    y = y_vec,
    learner = learner,
    learner_params = learner_params,
    folds = folds,
    fold_type = fold_type,
    block_id = if (!is.null(block_column)) df[[block_column]] else NULL,
    second_stage_controls = stage_controls_mat,
    seed = seed,
    mm_c = mm_c,
    max_iter = max_iter,
    tolerance = tolerance,
    ci_level = ci_level,
    call = original_call
  )

  result$controls <- if (length(controls) > 0L) controls else "(none)"
  result$second_stage_controls <- colnames(stage_controls_mat) %||% character(0)
  result$fold_type <- fold_type
  result$block_id <- if (!is.null(block_column)) df[[block_column]] else NULL
  result$spatial_k <- spatial_k
  result$spatial_by <- spatial_by

  if (identical(se_type, "bootstrap")) {
    boot <- .bootstrap_mm_dml_joint(
      x = if (ncol(x_mat) == 0L) matrix(0, nrow = n, ncol = 1L) else x_mat,
      d = d_mat,
      y = y_vec,
      learner = learner,
      learner_params = learner_params,
      folds = folds,
      fold_type = fold_type,
      block_id = if (!is.null(block_column)) df[[block_column]] else NULL,
      second_stage_controls = stage_controls_mat,
      seed = seed,
      mm_c = mm_c,
      max_iter = max_iter,
      tolerance = tolerance,
      ci_level = ci_level,
      bootstrap_replications = bootstrap_replications,
      bootstrap_seed = bootstrap_seed %||% (seed + 20000L),
      progress = progress,
      n_cores = n_cores
    )
    result <- .apply_bootstrap_to_joint_result(result, boot, primary_term = treatment[[1]])
  } else {
    result$bootstrap_estimates <- NULL
    result$bootstrap_replications <- 0L
    result$bootstrap_ci <- NULL
  }

  result$se_type <- se_type
  result$n_cores <- as.integer(n_cores)
  result$progress <- isTRUE(progress)
  class(result) <- "robustcause_mm_dml"
  result
}

.supports_native_mm_dml <- function(learner, d, fold_type) {
  is.character(learner) &&
    length(learner) == 1L &&
    learner %in% c("lasso", "elastic_net", "random_forest", "hist_gradient_boosting") &&
    identical(fold_type, "random") &&
    ncol(as.matrix(d)) == 1L
}

.run_native_mm_dml_matrix <- function(x,
                                      d,
                                      y,
                                      learner,
                                      learner_params,
                                      se_type,
                                      bootstrap_replications,
                                      bootstrap_seed,
                                      progress,
                                      n_cores,
                                      folds,
                                      seed,
                                      mm_c,
                                      max_iter,
                                      tolerance,
                                      ci_level,
                                      original_call) {
  fit <- .Call(
    rc_r_fit_mm_dml,
    unname(as.matrix(x)),
    as.numeric(d),
    as.numeric(y),
    learner,
    learner_params,
    as.integer(folds),
    as.integer(seed),
    as.numeric(mm_c),
    as.integer(max_iter),
    as.numeric(tolerance),
    as.numeric(ci_level)
  )

  fit$backend <- "cpp"
  fit$learner_name <- learner
  fit$learner_params <- learner_params
  fit$call <- original_call %||% match.call()
  fit$treatment_name <- "d"
  fit$treatment_names <- "d"
  fit$nobs <- nrow(as.matrix(x))
  fit$folds <- as.integer(folds)
  fit$fold_type <- "random"
  fit$progress <- isTRUE(progress)
  fit$n_cores <- as.integer(n_cores)
  fit$se_type <- se_type
  fit$bootstrap_replications <- if (identical(se_type, "bootstrap")) as.integer(bootstrap_replications) else 0L
  fit$bootstrap_estimates <- NULL
  fit$bootstrap_intercept_se <- NA_real_

  if (identical(se_type, "bootstrap")) {
    boot <- .bootstrap_mm_dml_native(
      x = as.matrix(x),
      d = as.numeric(d),
      y = as.numeric(y),
      learner = learner,
      learner_params = learner_params,
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
  }

  class(fit) <- "robustcause_mm_dml"
  fit
}

.run_mm_dml_joint_stage <- function(x,
                                    d,
                                    y,
                                    learner,
                                    learner_params,
                                    folds,
                                    fold_type,
                                    block_id,
                                    second_stage_controls,
                                    seed,
                                    mm_c,
                                    max_iter,
                                    tolerance,
                                    ci_level,
                                    call) {
  residuals <- .crossfit_joint_residuals(
    x = x,
    d = d,
    y = y,
    learner = learner,
    learner_params = learner_params,
    folds = folds,
    fold_type = fold_type,
    block_id = block_id,
    seed = seed
  )

  design <- cbind(residuals$d_res, second_stage_controls)
  if (is.null(colnames(design))) {
    colnames(design) <- c(
      paste0("d", seq_len(ncol(residuals$d_res))),
      if (!is.null(second_stage_controls)) paste0("z", seq_len(ncol(second_stage_controls))) else character(0)
    )
  }
  stage_fit <- fit_rlm(
    x = design,
    y = residuals$y_res,
    method = "mm",
    psi = "tukey_bisquare",
    tuning = mm_c,
    maxit = max_iter,
    tol = tolerance,
    add_intercept = TRUE
  )
  stage_fit$call <- call
  stage_summary <- summary(stage_fit, hc_type = "HC3", level = ci_level)
  stage_vcov <- vcov(stage_fit, hc_type = "HC3")
  dimnames(stage_vcov) <- list(stage_fit$coef_names, stage_fit$coef_names)
  coef_table <- stage_summary$coefficients[, c("Estimate", "Std. Error", "z value", "Pr(>|z|)"), drop = FALSE]
  colnames(coef_table)[3] <- "t value"
  beta <- stage_fit$coef
  names(beta) <- stage_fit$coef_names
  stage_design <- cbind("(Intercept)" = 1, unname(as.matrix(design)))

  treatment_terms <- colnames(residuals$d_res)
  if (is.null(treatment_terms)) {
    treatment_terms <- paste0("d", seq_len(ncol(residuals$d_res)))
  }
  primary_term <- treatment_terms[[1]]
  ci <- confint(stage_fit, hc_type = "HC3", level = ci_level)

  list(
    call = call,
    learner_name = .learner_label(learner),
    learner = learner,
    learner_params = learner_params,
    backend = "r-joint",
    folds = as.integer(folds),
    nobs = nrow(as.matrix(x)),
    treatment_name = if (length(treatment_terms) == 1L) treatment_terms[[1]] else paste(treatment_terms, collapse = ", "),
    treatment_names = treatment_terms,
    y_residuals = residuals$y_res,
    d_residuals = residuals$d_res,
    fold_id = residuals$fold_id,
    stage_residuals = residuals$y_res - as.numeric(stage_design %*% beta),
    stage_fit = stage_fit,
    stage_vcov = stage_vcov,
    coefficients = coef_table,
    ci_level = ci_level,
    converged = isTRUE(stage_fit$converged),
    iterations = stage_fit$iterations,
    tau_hat = unname(beta[primary_term]),
    intercept_hat = unname(beta[["(Intercept)"]]),
    std_error = unname(coef_table[primary_term, "Std. Error"]),
    ci_lower = unname(ci[primary_term, 1]),
    ci_upper = unname(ci[primary_term, 2]),
    objective = mean((residuals$y_res - as.numeric(stage_design %*% beta))^2),
    kept_fraction = 1.0,
    residual_outcome_scale = stats::sd(residuals$y_res),
    residual_treatment_scale = if (is.matrix(residuals$d_res)) stats::sd(residuals$d_res[, 1]) else stats::sd(residuals$d_res)
  )
}

.bootstrap_mm_dml_native <- function(x,
                                     d,
                                     y,
                                     learner,
                                     learner_params,
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
  estimates_list <- .parallel_bootstrap_apply(
    X = seq_len(bootstrap_replications),
    n_cores = n_cores,
    progress = progress,
    label = "Bootstrap",
    fun = function(b) {
      fit_b <- try(
        .Call(
          rc_r_fit_mm_dml,
          unname(as.matrix(x[draws[[b]], , drop = FALSE])),
          as.numeric(d[draws[[b]]]),
          as.numeric(y[draws[[b]]]),
          learner,
          learner_params,
          as.integer(folds),
          as.integer(seed + b),
          as.numeric(mm_c),
          as.integer(max_iter),
          as.numeric(tolerance),
          as.numeric(ci_level)
        ),
        silent = TRUE
      )
      if (inherits(fit_b, "try-error")) {
        return(c(NA_real_, NA_real_))
      }
      c(fit_b$intercept_hat, fit_b$tau_hat)
    }
  )
  estimates <- do.call(rbind, estimates_list)
  colnames(estimates) <- c("(Intercept)", "d")
  estimates <- estimates[stats::complete.cases(estimates), , drop = FALSE]
  if (nrow(estimates) == 0L) {
    stop("All native MM-DML bootstrap replications failed.", call. = FALSE)
  }
  alpha <- (1 - ci_level) / 2

  list(
    estimates = estimates,
    std_error = as.list(setNames(apply(estimates, 2, stats::sd), colnames(estimates))),
    ci_lower = as.list(setNames(apply(estimates, 2, stats::quantile, probs = alpha), colnames(estimates))),
    ci_upper = as.list(setNames(apply(estimates, 2, stats::quantile, probs = 1 - alpha), colnames(estimates)))
  )
}

.bootstrap_mm_dml_joint <- function(x,
                                    d,
                                    y,
                                    learner,
                                    learner_params,
                                    folds,
                                    fold_type,
                                    block_id,
                                    second_stage_controls,
                                    seed,
                                    mm_c,
                                    max_iter,
                                    tolerance,
                                    ci_level,
                                    bootstrap_replications,
                                    bootstrap_seed,
                                    progress,
                                    n_cores) {
  set.seed(bootstrap_seed)
  n <- length(y)
  draws <- lapply(seq_len(bootstrap_replications), function(...) sample.int(n, size = n, replace = TRUE))
  coef_names <- NULL

  estimates_list <- .parallel_bootstrap_apply(
    X = seq_len(bootstrap_replications),
    n_cores = n_cores,
    progress = progress,
    label = "Bootstrap",
    fun = function(b) {
      fit_b <- try(
        .run_mm_dml_joint_stage(
          x = x[draws[[b]], , drop = FALSE],
          d = d[draws[[b]], , drop = FALSE],
          y = y[draws[[b]]],
          learner = learner,
          learner_params = learner_params,
          folds = folds,
          fold_type = fold_type,
          block_id = if (is.null(block_id)) NULL else block_id[draws[[b]]],
          second_stage_controls = if (is.null(second_stage_controls)) NULL else second_stage_controls[draws[[b]], , drop = FALSE],
          seed = seed + b,
          mm_c = mm_c,
          max_iter = max_iter,
          tolerance = tolerance,
          ci_level = ci_level,
          call = NULL
        ),
        silent = TRUE
      )
      if (inherits(fit_b, "try-error")) {
        return(NA)
      }
      beta <- fit_b$stage_fit$coef
      names(beta) <- fit_b$stage_fit$coef_names
      beta
    }
  )

  estimates_list <- estimates_list[vapply(estimates_list, function(x) !(length(x) == 1L && is.na(x)), logical(1))]
  if (length(estimates_list) == 0L) {
    estimates_list <- .run_with_progress(
      X = seq_len(bootstrap_replications),
      fun = function(b) {
        fit_b <- try(
          .run_mm_dml_joint_stage(
            x = x[draws[[b]], , drop = FALSE],
            d = d[draws[[b]], , drop = FALSE],
            y = y[draws[[b]]],
            learner = learner,
            learner_params = learner_params,
            folds = folds,
            fold_type = fold_type,
            block_id = if (is.null(block_id)) NULL else block_id[draws[[b]]],
            second_stage_controls = if (is.null(second_stage_controls)) NULL else second_stage_controls[draws[[b]], , drop = FALSE],
            seed = seed + b,
            mm_c = mm_c,
            max_iter = max_iter,
            tolerance = tolerance,
            ci_level = ci_level,
            call = NULL
          ),
          silent = TRUE
        )
        if (inherits(fit_b, "try-error")) {
          return(NA)
        }
        beta <- fit_b$stage_fit$coef
        names(beta) <- fit_b$stage_fit$coef_names
        beta
      },
      progress = FALSE,
      label = "Bootstrap"
    )
    estimates_list <- estimates_list[vapply(estimates_list, function(x) !(length(x) == 1L && is.na(x)), logical(1))]
    if (length(estimates_list) == 0L) {
      stop("All MM-DML bootstrap replications failed.", call. = FALSE)
    }
  }
  coef_names <- unique(unlist(lapply(estimates_list, names), use.names = FALSE))
  estimates <- do.call(rbind, lapply(estimates_list, function(beta) {
    out <- setNames(rep(NA_real_, length(coef_names)), coef_names)
    out[names(beta)] <- beta
    out
  }))

  alpha <- (1 - ci_level) / 2
  std_error <- setNames(as.list(apply(estimates, 2, stats::sd, na.rm = TRUE)), colnames(estimates))
  lower <- setNames(as.list(apply(estimates, 2, stats::quantile, probs = alpha, na.rm = TRUE)), colnames(estimates))
  upper <- setNames(as.list(apply(estimates, 2, stats::quantile, probs = 1 - alpha, na.rm = TRUE)), colnames(estimates))
  ci <- lapply(colnames(estimates), function(nm) c(lower[[nm]], upper[[nm]]))
  names(ci) <- colnames(estimates)
  list(estimates = estimates, std_error = std_error, ci = ci)
}

.apply_bootstrap_to_joint_result <- function(result, boot, primary_term) {
  coef_names <- rownames(result$coefficients)
  p_col <- grep("^Pr\\(", colnames(result$coefficients), value = TRUE)[1]
  stat_col <- if ("t value" %in% colnames(result$coefficients)) "t value" else grep("value$", colnames(result$coefficients), value = TRUE)[1]
  for (nm in intersect(names(boot$std_error), coef_names)) {
    se <- boot$std_error[[nm]]
    result$coefficients[nm, "Std. Error"] <- se
    result$coefficients[nm, stat_col] <- result$coefficients[nm, "Estimate"] / se
    result$coefficients[nm, p_col] <- 2 * stats::pnorm(abs(result$coefficients[nm, stat_col]), lower.tail = FALSE)
  }
  result$bootstrap_estimates <- boot$estimates
  result$bootstrap_replications <- nrow(boot$estimates)
  result$bootstrap_ci <- boot$ci
  if (primary_term %in% rownames(result$coefficients)) {
    result$std_error <- unname(result$coefficients[primary_term, "Std. Error"])
    result$ci_lower <- unname(boot$ci[[primary_term]][1])
    result$ci_upper <- unname(boot$ci[[primary_term]][2])
  }
  result
}

.parallel_bootstrap_apply <- function(X, n_cores, fun, progress = interactive(), label = "Bootstrap") {
  worker_count <- max(1L, as.integer(n_cores %||% 1L))
  total <- length(X)
  if (worker_count <= 1L || total <= 1L) {
    return(.run_with_progress(X, fun, progress = progress, label = label))
  }

  cl <- parallel::makePSOCKcluster(worker_count)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterCall(cl, function(paths) .libPaths(paths), .libPaths())
  parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages(library(robustcause))
    NULL
  })

  chunk_size <- max(1L, ceiling(total / max(1L, worker_count * 8L)))
  chunks <- split(X, ceiling(seq_along(X) / chunk_size))
  start_time <- Sys.time()
  completed <- 0L
  pb <- NULL
  if (isTRUE(progress)) {
    pb <- utils::txtProgressBar(min = 0, max = total, style = 3)
    on.exit(close(pb), add = TRUE)
  }
  results <- vector("list", total)
  for (chunk in chunks) {
    chunk_results <- parallel::parLapply(cl, chunk, fun)
    results[completed + seq_along(chunk_results)] <- chunk_results
    completed <- completed + length(chunk_results)
    .update_progress_bar(pb, completed, total, start_time, label)
  }
  results
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

.crossfit_joint_residuals <- function(x,
                                      d,
                                      y,
                                      learner,
                                      learner_params,
                                      folds = 2L,
                                      fold_type = c("random", "spatial_block"),
                                      block_id = NULL,
                                      seed = 123L) {
  n <- length(y)
  d <- as.matrix(d)
  if (is.null(colnames(d))) {
    colnames(d) <- if (ncol(d) == 1L) "d" else paste0("d", seq_len(ncol(d)))
  }
  fold_id <- .make_fold_id(
    n = n,
    folds = folds,
    fold_type = fold_type,
    block_id = block_id,
    seed = seed
  )
  y_hat <- rep(NA_real_, n)
  d_hat <- matrix(NA_real_, nrow = n, ncol = ncol(d), dimnames = list(NULL, colnames(d)))
  for (fold in seq_len(folds)) {
    train_idx <- which(fold_id != fold)
    test_idx <- which(fold_id == fold)
    y_hat[test_idx] <- .predict_nuisance(
      x_train = x[train_idx, , drop = FALSE],
      y_train = y[train_idx],
      x_test = x[test_idx, , drop = FALSE],
      learner = learner,
      learner_params = learner_params,
      seed = seed + fold,
      role = "y"
    )
    for (j in seq_len(ncol(d))) {
      d_hat[test_idx, j] <- .predict_nuisance(
        x_train = x[train_idx, , drop = FALSE],
        y_train = d[train_idx, j],
        x_test = x[test_idx, , drop = FALSE],
        learner = learner,
        learner_params = learner_params,
        seed = seed + 100L * j + fold,
        role = "t"
      )
    }
  }
  list(
    y_res = as.numeric(y - y_hat),
    d_res = d - d_hat,
    fold_id = fold_id
  )
}

.make_fold_id <- function(n, folds, fold_type = c("random", "spatial_block"), block_id = NULL, seed = 123L) {
  fold_type <- match.arg(fold_type)
  if (folds < 2L) {
    stop("`folds` must be at least 2.", call. = FALSE)
  }
  if (identical(fold_type, "random")) {
    set.seed(seed)
    return(sample(rep(seq_len(folds), length.out = n)))
  }
  if (is.null(block_id)) {
    stop("`block_id` must be provided when `fold_type = 'spatial_block'`.", call. = FALSE)
  }
  if (length(block_id) != n) {
    stop("`block_id` must have the same length as the number of observations.", call. = FALSE)
  }
  block_factor <- as.factor(block_id)
  if (anyNA(block_factor)) {
    stop("`block_id` cannot contain missing values.", call. = FALSE)
  }
  unique_blocks <- levels(block_factor)
  if (length(unique_blocks) < folds) {
    stop("Need at least as many unique spatial blocks as folds.", call. = FALSE)
  }
  set.seed(seed)
  block_folds <- sample(rep(seq_len(folds), length.out = length(unique_blocks)))
  names(block_folds) <- unique_blocks
  unname(block_folds[as.character(block_factor)])
}

.normalize_block_id_argument <- function(block_id, data) {
  if (is.null(block_id)) {
    return(NULL)
  }
  if (is.character(block_id) && length(block_id) == 1L && block_id %in% names(data)) {
    return(block_id)
  }
  NULL
}

.predict_nuisance <- function(x_train, y_train, x_test, learner, learner_params, seed, role = c("y", "t")) {
  role <- match.arg(role)
  x_train_df <- .as_feature_frame(x_train)
  x_test_df <- .as_feature_frame(x_test, template = x_train_df)
  if (ncol(x_train_df) == 0L) {
    return(rep(mean(y_train), nrow(x_test_df)))
  }

  if (inherits(learner, "robustcause_custom_learner") || (is.list(learner) && !is.character(learner))) {
    return(.predict_custom_nuisance(x_train_df, y_train, x_test_df, learner, seed = seed, role = role))
  }

  learner_name <- .normalize_learner_name(learner)
  predictor <- .nuisance_fitters(learner_params)[[learner_name]]
  if (is.null(predictor)) {
    stop("Unsupported learner: ", learner_name, call. = FALSE)
  }
  predictor(x_train_df, y_train, x_test_df, seed)
}

.predict_custom_nuisance <- function(x_train, y_train, x_test, learner, seed, role) {
  spec <- .resolve_custom_learner_spec(learner, role)
  if (is.function(spec)) {
    out <- spec(x_train = x_train, y_train = y_train, x_test = x_test, seed = seed, role = role)
  } else if (is.list(spec) && is.function(spec$fit) && is.function(spec$predict)) {
    model <- spec$fit(x_train = x_train, y_train = y_train, seed = seed, role = role)
    out <- spec$predict(model = model, x_test = x_test, seed = seed, role = role)
  } else {
    stop("Custom learner specs must be prediction functions or lists with `fit` and `predict` functions.", call. = FALSE)
  }
  out <- as.numeric(out)
  if (length(out) != nrow(x_test)) {
    stop("Custom learner returned ", length(out), " predictions for ", nrow(x_test), " rows.", call. = FALSE)
  }
  out
}

.resolve_custom_learner_spec <- function(learner, role) {
  if (inherits(learner, "robustcause_custom_learner")) {
    return(if (identical(role, "y")) learner$model_y else learner$model_t)
  }
  if (is.list(learner) && all(c("model_y", "model_t") %in% names(learner))) {
    return(if (identical(role, "y")) learner$model_y else learner$model_t)
  }
  if (is.function(learner)) {
    return(learner)
  }
  learner
}

.as_feature_frame <- function(x, template = NULL) {
  if (is.data.frame(x)) {
    out <- x
  } else {
    out <- as.data.frame(x)
  }
  if (!is.null(template)) {
    names(out) <- names(template)
  } else if (ncol(out) > 0L && all(grepl("^V[0-9]+$", names(out)))) {
    names(out) <- paste0("x", seq_len(ncol(out)))
  }
  out
}

.normalize_mm_dml_learner <- function(learner) {
  if (inherits(learner, "robustcause_custom_learner")) {
    return(learner)
  }
  if (is.character(learner)) {
    if (length(learner) != 1L) {
      stop("`learner` must be a single string or a custom learner object.", call. = FALSE)
    }
    return(.normalize_learner_name(learner))
  }
  if (is.function(learner)) {
    return(make_custom_learner(learner))
  }
  if (is.list(learner)) {
    if (all(c("model_y", "model_t") %in% names(learner))) {
      class(learner) <- unique(c("robustcause_custom_learner", class(learner)))
      learner$label <- learner$label %||% "custom"
      return(learner)
    }
    if (all(c("fit", "predict") %in% names(learner))) {
      return(make_custom_learner(learner))
    }
  }
  stop("Unsupported `learner`. Use a named learner string or `make_custom_learner(...)`.", call. = FALSE)
}

.normalize_mm_dml_learner_params <- function(learner, learner_params) {
  if (is.null(learner_params)) {
    return(list())
  }
  if (!is.list(learner_params)) {
    stop("`learner_params` must be NULL or a named list.", call. = FALSE)
  }
  if (length(learner_params) == 0L) {
    return(list())
  }
  if (is.null(names(learner_params)) || anyNA(names(learner_params)) || any(names(learner_params) == "")) {
    stop("`learner_params` must be a named list.", call. = FALSE)
  }

  if (!is.character(learner) || length(learner) != 1L) {
    return(learner_params)
  }

  defaults <- switch(
    learner,
    lasso = list(
      cv_folds = 3L,
      n_lambda = 100L,
      lambda_min_ratio = 1e-4,
      standardize = TRUE,
      lambda_grid = NULL,
      use_lambda_1se = FALSE,
      max_iter = 250L,
      tolerance = 1e-6
    ),
    elastic_net = list(
      cv_folds = 3L,
      n_lambda = 100L,
      lambda_min_ratio = 1e-4,
      standardize = TRUE,
      lambda_grid = NULL,
      use_lambda_1se = FALSE,
      max_iter = 250L,
      tolerance = 1e-6,
      l1_ratios = c(0.1, 0.5, 0.9, 0.95, 1.0)
    ),
    random_forest = list(
      n_estimators = 300L,
      max_depth = 8L,
      min_samples_split = 10L,
      min_samples_leaf = 5L,
      max_features = NULL,
      max_features_fraction = NULL,
      sample_fraction = 0.632,
      bootstrap = TRUE,
      replacement = TRUE,
      split_candidates = 16L,
      compute_oob = TRUE,
      compute_importance = TRUE
    ),
    hist_gradient_boosting = list(
      max_depth = 6L,
      min_samples_leaf = 5L,
      max_features = NULL,
      max_iter = 200L,
      learning_rate = 0.05
    ),
    list()
  )
  modifyList(defaults, learner_params)
}

.normalize_learner_name <- function(learner) {
  aliases <- c(
    rf = "random_forest",
    random_forest = "random_forest",
    ols = "ols",
    lasso = "lasso",
    elastic_net = "elastic_net",
    xgboost = "xgboost",
    gam = "gam",
    spatial_gam = "spatial_gam",
    hist_gradient_boosting = "hist_gradient_boosting"
  )
  out <- unname(aliases[[learner]])
  if (is.null(out)) {
    stop("Unsupported learner: ", learner, call. = FALSE)
  }
  out
}

.learner_label <- function(learner) {
  if (inherits(learner, "robustcause_custom_learner")) {
    return(learner$label %||% "custom")
  }
  if (is.character(learner)) {
    return(learner)
  }
  "custom"
}

.nuisance_fitters <- function(learner_params = NULL) {
  list(
    ols = function(x_train, y_train, x_test, seed) {
      fit <- stats::lm.fit(x = cbind(1, data.matrix(x_train)), y = y_train)
      as.numeric(cbind(1, data.matrix(x_test)) %*% fit$coefficients)
    },
    lasso = function(x_train, y_train, x_test, seed) {
      if (!requireNamespace("glmnet", quietly = TRUE)) {
        stop("The `glmnet` package is required for learner = 'lasso'.", call. = FALSE)
      }
      params <- modifyList(list(
        cv_folds = 3L,
        n_lambda = 100L,
        lambda_min_ratio = 1e-4,
        standardize = TRUE,
        lambda_grid = NULL,
        use_lambda_1se = FALSE
      ), learner_params %||% list())
      fit <- glmnet::cv.glmnet(
        x = data.matrix(x_train),
        y = y_train,
        alpha = 1,
        family = "gaussian",
        nfolds = params$cv_folds,
        nlambda = params$n_lambda,
        lambda.min.ratio = params$lambda_min_ratio,
        standardize = isTRUE(params$standardize),
        lambda = params$lambda_grid
      )
      as.numeric(stats::predict(
        fit,
        newx = data.matrix(x_test),
        s = if (isTRUE(params$use_lambda_1se)) "lambda.1se" else "lambda.min"
      ))
    },
    elastic_net = function(x_train, y_train, x_test, seed) {
      if (!requireNamespace("glmnet", quietly = TRUE)) {
        stop("The `glmnet` package is required for learner = 'elastic_net'.", call. = FALSE)
      }
      params <- modifyList(list(
        cv_folds = 3L,
        n_lambda = 100L,
        lambda_min_ratio = 1e-4,
        standardize = TRUE,
        lambda_grid = NULL,
        alpha = 0.5,
        use_lambda_1se = FALSE
      ), learner_params %||% list())
      fit <- glmnet::cv.glmnet(
        x = data.matrix(x_train),
        y = y_train,
        alpha = params$alpha,
        family = "gaussian",
        nfolds = params$cv_folds,
        nlambda = params$n_lambda,
        lambda.min.ratio = params$lambda_min_ratio,
        standardize = isTRUE(params$standardize),
        lambda = params$lambda_grid
      )
      as.numeric(stats::predict(
        fit,
        newx = data.matrix(x_test),
        s = if (isTRUE(params$use_lambda_1se)) "lambda.1se" else "lambda.min"
      ))
    },
    random_forest = function(x_train, y_train, x_test, seed) {
      if (!requireNamespace("ranger", quietly = TRUE)) {
        stop("The `ranger` package is required for learner = 'random_forest'.", call. = FALSE)
      }
      params <- modifyList(list(
        n_estimators = 300L,
        min_samples_leaf = 5L,
        max_depth = NULL,
        max_features = NULL,
        sample_fraction = NULL,
        replacement = TRUE
      ), learner_params %||% list())
      fit <- ranger::ranger(
        x = x_train,
        y = y_train,
        num.trees = params$n_estimators,
        min.node.size = params$min_samples_leaf,
        max.depth = params$max_depth,
        mtry = params$max_features %||% max(1L, floor(sqrt(ncol(x_train)))),
        sample.fraction = params$sample_fraction,
        replace = isTRUE(params$replacement),
        seed = seed,
        num.threads = 1
      )
      as.numeric(stats::predict(fit, data = x_test)$predictions)
    },
    hist_gradient_boosting = function(x_train, y_train, x_test, seed) {
      if (!requireNamespace("xgboost", quietly = TRUE)) {
        return(.nuisance_fitters(learner_params)$random_forest(x_train, y_train, x_test, seed))
      }
      dtrain <- xgboost::xgb.DMatrix(data = data.matrix(x_train), label = y_train)
      dtest <- xgboost::xgb.DMatrix(data = data.matrix(x_test))
      fit <- xgboost::xgb.train(
        params = list(
          objective = "reg:squarederror",
          eta = 0.05,
          max_depth = 6,
          subsample = 0.8,
          colsample_bytree = 0.8,
          nthread = 1
        ),
        data = dtrain,
        nrounds = 200,
        verbose = 0
      )
      as.numeric(stats::predict(fit, dtest))
    },
    xgboost = function(x_train, y_train, x_test, seed) {
      .nuisance_fitters(learner_params)$hist_gradient_boosting(x_train, y_train, x_test, seed)
    },
    gam = function(x_train, y_train, x_test, seed) {
      if (!requireNamespace("mgcv", quietly = TRUE)) {
        stop("The `mgcv` package is required for learner = 'gam'.", call. = FALSE)
      }
      fit <- mgcv::gam(.build_gam_formula(x_train), data = cbind(y = y_train, x_train), method = "REML")
      as.numeric(stats::predict(fit, newdata = x_test))
    },
    spatial_gam = function(x_train, y_train, x_test, seed) {
      if (!all(c("X", "Y") %in% names(x_train))) {
        return(.nuisance_fitters(learner_params)$gam(x_train, y_train, x_test, seed))
      }
      if (!requireNamespace("mgcv", quietly = TRUE)) {
        stop("The `mgcv` package is required for learner = 'spatial_gam'.", call. = FALSE)
      }
      other_vars <- setdiff(names(x_train), c("X", "Y"))
      rhs <- c(
        "s(X, Y, bs = 'tp')",
        .build_gam_rhs(x_train[, other_vars, drop = FALSE])
      )
      fit <- mgcv::gam(
        stats::as.formula(paste("y ~", paste(rhs[nzchar(rhs)], collapse = " + "))),
        data = cbind(y = y_train, x_train),
        method = "REML"
      )
      as.numeric(stats::predict(fit, newdata = x_test))
    }
  )
}

.build_gam_formula <- function(x_train) {
  rhs <- .build_gam_rhs(x_train)
  stats::as.formula(paste("y ~", paste(rhs[nzchar(rhs)], collapse = " + ")))
}

.build_gam_rhs <- function(x_train) {
  if (ncol(x_train) == 0L) {
    return("1")
  }
  vapply(names(x_train), function(nm) {
    column <- x_train[[nm]]
    if (is.factor(column) || is.character(column)) {
      return(nm)
    }
    unique_count <- length(unique(stats::na.omit(column)))
    if (unique_count >= 10L) {
      sprintf("s(%s, bs = 'cs')", nm)
    } else {
      nm
    }
  }, character(1))
}

summary.robustcause_mm_dml <- function(object, ...) {
  if (!is.null(object$coefficients)) {
    coef_table <- object$coefficients
  } else {
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
    treatment_label <- object$treatment_name %||% object$treatment_names %||% "d"
    if (length(treatment_label) != 1L || is.na(treatment_label) || !nzchar(treatment_label)) {
      treatment_label <- "d"
    }
    coef_est <- c(object$intercept_hat, object$tau_hat)
    coef_se <- c(intercept_se, object$std_error)
    coef_t <- coef_est / coef_se
    coef_p <- 2 * stats::pt(abs(coef_t), df = df_residual, lower.tail = FALSE)
    coef_table <- matrix(
      c(coef_est, coef_se, coef_t, coef_p),
      nrow = 2L,
      dimnames = list(
        c("(Intercept)", treatment_label),
        c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
      )
    )
  }

  structure(
    list(
      call = object$call,
      coefficients = coef_table,
      learner = object$learner_name %||% .learner_label(object$learner),
      backend = object$backend %||% "cpp",
      iterations = object$iterations,
      converged = object$converged,
      nobs = object$nobs,
      folds = object$folds,
      fold_type = object$fold_type %||% "random",
      se_type = object$se_type %||% "analytic",
      bootstrap_replications = object$bootstrap_replications %||% 0L,
      n_cores = object$n_cores %||% 1L,
      learner_details = object$learner_details %||% NULL
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
  if (!is.null(x$learner_details)) {
    detail_lines <- .format_mm_dml_learner_details(x$learner_details)
    if (length(detail_lines) > 0L) {
      cat("Learner details:\n")
      cat(paste0("  ", detail_lines, collapse = "\n"), "\n")
    }
  }
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

plot_marginal_effect <- function(fit,
                                 focal,
                                 interaction,
                                 moderator_values = c(0, 1),
                                 moderator_labels = NULL,
                                 ci_level = 0.95,
                                 xlab = NULL,
                                 ylab = NULL,
                                 main = NULL,
                                 pch = 19,
                                 col = "black",
                                 ci_col = "black",
                                 ci_lty = 1,
                                 ci_lwd = 1.2,
                                 ...) {
  if (!inherits(fit, "robustcause_mm_dml")) {
    stop("`fit` must be an object returned by `fit_mm_dml()`.", call. = FALSE)
  }
  coefficients <- .extract_mm_dml_coefficients(fit)
  missing_terms <- setdiff(c(focal, interaction), rownames(coefficients))
  if (length(missing_terms) > 0L) {
    stop("Coefficient(s) not found in fit: ", paste(missing_terms, collapse = ", "), call. = FALSE)
  }
  moderator_values <- as.numeric(moderator_values)
  if (length(moderator_values) < 1L) {
    stop("`moderator_values` must contain at least one numeric value.", call. = FALSE)
  }
  if (!is.null(moderator_labels) && length(moderator_labels) != length(moderator_values)) {
    stop("`moderator_labels` must have the same length as `moderator_values`.", call. = FALSE)
  }

  effect_df <- data.frame(
    moderator_value = moderator_values,
    marginal_effect = coefficients[focal, "Estimate"] + moderator_values * coefficients[interaction, "Estimate"]
  )

  if (!is.null(fit$bootstrap_estimates) &&
      is.matrix(fit$bootstrap_estimates) &&
      all(c(focal, interaction) %in% colnames(fit$bootstrap_estimates))) {
    alpha <- (1 - ci_level) / 2
    draws <- fit$bootstrap_estimates[, focal] + moderator_values %o% fit$bootstrap_estimates[, interaction]
    effect_df$ci_lower <- apply(draws, 1, stats::quantile, probs = alpha, na.rm = TRUE)
    effect_df$ci_upper <- apply(draws, 1, stats::quantile, probs = 1 - alpha, na.rm = TRUE)
    effect_df$se <- apply(draws, 1, stats::sd, na.rm = TRUE)
    interval_method <- "bootstrap"
  } else {
    vcov_mat <- .extract_mm_dml_vcov(fit, coefficients)
    var_focal <- vcov_mat[focal, focal]
    var_interaction <- vcov_mat[interaction, interaction]
    cov_term <- vcov_mat[focal, interaction]
    effect_df$se <- sqrt(pmax(0, var_focal + moderator_values^2 * var_interaction + 2 * moderator_values * cov_term))
    z <- stats::qnorm(1 - (1 - ci_level) / 2)
    effect_df$ci_lower <- effect_df$marginal_effect - z * effect_df$se
    effect_df$ci_upper <- effect_df$marginal_effect + z * effect_df$se
    interval_method <- "analytic"
  }

  x_positions <- seq_along(moderator_values)
  xlab <- xlab %||% "Moderator value"
  ylab <- ylab %||% paste("Marginal effect of", focal)
  main <- main %||% sprintf("Marginal effect of %s", focal)
  ylim <- range(c(effect_df$ci_lower, effect_df$ci_upper), na.rm = TRUE)
  if (!all(is.finite(ylim))) {
    ylim <- range(effect_df$marginal_effect, na.rm = TRUE)
  }

  graphics::plot(
    x_positions,
    effect_df$marginal_effect,
    ylim = ylim,
    xlab = xlab,
    ylab = ylab,
    main = main,
    xaxt = "n",
    pch = pch,
    col = col,
    ...
  )
  axis_labels <- moderator_labels %||% format(moderator_values, trim = TRUE)
  graphics::axis(1, at = x_positions, labels = axis_labels)
  graphics::arrows(
    x0 = x_positions,
    y0 = effect_df$ci_lower,
    x1 = x_positions,
    y1 = effect_df$ci_upper,
    angle = 90,
    code = 3,
    length = 0.05,
    col = ci_col,
    lty = ci_lty,
    lwd = ci_lwd
  )
  graphics::abline(h = 0, lty = 2, col = "gray60")

  effect_df$focal <- focal
  effect_df$interaction <- interaction
  effect_df$interval_method <- interval_method
  invisible(effect_df)
}

.extract_mm_dml_coefficients <- function(fit) {
  if (!is.null(fit$coefficients)) {
    return(fit$coefficients)
  }
  coef_table <- cbind(
    Estimate = c(fit$intercept_hat, fit$tau_hat),
    `Std. Error` = c(fit$bootstrap_intercept_se %||% NA_real_, fit$std_error)
  )
  rownames(coef_table) <- c("(Intercept)", fit$treatment_name %||% "d")
  coef_table
}

.extract_mm_dml_vcov <- function(fit, coefficients) {
  coef_names <- rownames(coefficients)
  if (!is.null(fit$stage_vcov)) {
    return(fit$stage_vcov[coef_names, coef_names, drop = FALSE])
  }
  if (all(c("(Intercept)", fit$treatment_name %||% "d") %in% coef_names) && length(coef_names) == 2L) {
    se_vec <- coefficients[, "Std. Error"]
    names(se_vec) <- coef_names
    out <- matrix(0, nrow = 2, ncol = 2, dimnames = list(coef_names, coef_names))
    diag(out) <- se_vec^2
    return(out)
  }
  stop("Analytic covariance matrix is not available for the requested marginal effect.", call. = FALSE)
}

.format_mm_dml_learner_details <- function(details) {
  if (!is.list(details)) {
    return(character(0))
  }
  out <- character(0)
  for (role in intersect(c("y", "d"), names(details))) {
    info <- details[[role]]
    if (!is.list(info)) {
      next
    }
    parts <- character(0)
    if (!is.null(info$selected_lambda) && is.finite(info$selected_lambda)) {
      parts <- c(parts, sprintf("lambda=%.4g", info$selected_lambda))
    }
    if (!is.null(info$lambda_min) && is.finite(info$lambda_min) && !identical(info$lambda_min, info$selected_lambda)) {
      parts <- c(parts, sprintf("lambda.min=%.4g", info$lambda_min))
    }
    if (!is.null(info$lambda_1se) && is.finite(info$lambda_1se)) {
      parts <- c(parts, sprintf("lambda.1se=%.4g", info$lambda_1se))
    }
    if (!is.null(info$cv_error) && is.finite(info$cv_error)) {
      parts <- c(parts, sprintf("cv_mse=%.4g", info$cv_error))
    }
    if (!is.null(info$nonzero) && is.finite(info$nonzero)) {
      parts <- c(parts, sprintf("nonzero=%d", as.integer(info$nonzero)))
    }
    if (!is.null(info$oob_error) && is.finite(info$oob_error)) {
      parts <- c(parts, sprintf("oob_mse=%.4g", info$oob_error))
    }
    if (!is.null(info$num_trees) && is.finite(info$num_trees)) {
      parts <- c(parts, sprintf("trees=%d", as.integer(info$num_trees)))
    }
    if (length(parts) > 0L) {
      out <- c(out, sprintf("%s-model: %s", role, paste(parts, collapse = " | ")))
    }
  }
  out
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}
