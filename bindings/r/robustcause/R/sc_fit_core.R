sc_name_native_result <- function(fit, prepared) {
  fit$weights <- stats::setNames(as.numeric(fit$weights), prepared$donors)

  actual <- sc_compute_actual_fit(fit$weights, prepared)
  fit$pre_residuals <- stats::setNames(actual$pre_residuals, prepared$pre_periods)
  fit$post_gaps <- stats::setNames(actual$post_gaps, prepared$post_periods)
  fit$synthetic_pre <- stats::setNames(actual$synthetic_pre, prepared$pre_periods)
  fit$synthetic_post <- stats::setNames(actual$synthetic_post, prepared$post_periods)
  fit$pre_rmspe <- sc_rmspe(actual$pre_residuals)

  fit$predictor_residuals <- stats::setNames(actual$predictor_residuals, prepared$predictor_names)
  fit$predictor_rmse <- if (length(actual$predictor_residuals) == 0L) NA_real_ else sc_rmspe(actual$predictor_residuals)

  if (!is.null(fit$startup_weights)) {
    fit$startup_weights <- stats::setNames(as.numeric(fit$startup_weights), prepared$donors)
  }
  if (!is.null(fit$robust_time_weights)) {
    all_weights <- as.numeric(fit$robust_time_weights)
    n_time <- length(prepared$pre_periods)
    n_pred <- length(prepared$predictor_names)

    fit$robust_time_weights <- stats::setNames(all_weights[seq_len(n_time)], prepared$pre_periods)
    if (n_pred > 0L) {
      fit$robust_predictor_weights <- stats::setNames(all_weights[n_time + seq_len(n_pred)], prepared$predictor_names)
    } else {
      fit$robust_predictor_weights <- numeric()
    }
  }

  fit$has_predictors <- length(prepared$predictor_names) > 0L
  fit$predictor_lambda <- prepared$predictor_lambda
  fit
}

sc_compute_actual_fit <- function(weights, prepared) {
  weights <- as.numeric(weights)
  synthetic_pre <- as.numeric(prepared$donors_pre %*% weights)
  synthetic_post <- as.numeric(prepared$donors_post %*% weights)
  predictor_synthetic <- if (length(prepared$predictor_names) == 0L) {
    numeric()
  } else {
    as.numeric(prepared$predictor_donors %*% weights)
  }

  list(
    synthetic_pre = synthetic_pre,
    synthetic_post = synthetic_post,
    pre_residuals = prepared$treated_pre - synthetic_pre,
    post_gaps = prepared$treated_post - synthetic_post,
    predictor_residuals = prepared$predictor_treated - predictor_synthetic
  )
}

sc_augmented_pre_data <- function(prepared) {
  y <- prepared$treated_pre
  x <- prepared$donors_pre

  if (length(prepared$predictor_names) == 0L || prepared$predictor_lambda <= 0) {
    return(list(y = y, x = x))
  }

  row_weights <- sqrt(prepared$predictor_lambda * prepared$predictor_weights)
  y_pred <- prepared$predictor_treated * row_weights
  x_pred <- prepared$predictor_donors * row_weights

  list(
    y = c(y, y_pred),
    x = rbind(x, x_pred)
  )
}

sc_fit_standard_core <- function(prepared,
                                 maxit = 1000L,
                                 tol = 1e-8,
                                 ridge = 1e-10) {
  aug <- sc_augmented_pre_data(prepared)
  empty_post <- matrix(numeric(), nrow = 0L, ncol = ncol(prepared$donors_pre))

  fit <- .Call(
    rc_r_fit_sc,
    aug$y,
    aug$x,
    numeric(),
    empty_post,
    as.integer(maxit),
    as.numeric(tol),
    as.numeric(ridge),
    prepared$predictor_treated,
    prepared$predictor_donors,
    prepared$predictor_weights
  )

  sc_name_native_result(fit, prepared)
}

sc_fit_mm_core <- function(prepared,
                           maxit = 1000L,
                           mm_max_iter = 25L,
                           subproblem_maxit = maxit,
                           tol = 1e-8,
                           subproblem_tol = tol,
                           l1_smoothing = 1e-6,
                           tukey_c = 4.685,
                           min_scale = 1e-8,
                           ridge = 1e-10,
                           min_time_weight = 1e-8) {
  aug <- sc_augmented_pre_data(prepared)
  empty_post <- matrix(numeric(), nrow = 0L, ncol = ncol(prepared$donors_pre))

  fit <- .Call(
    rc_r_fit_mm_sc,
    aug$y,
    aug$x,
    numeric(),
    empty_post,
    as.integer(maxit),
    as.integer(mm_max_iter),
    as.integer(subproblem_maxit),
    as.numeric(tol),
    as.numeric(subproblem_tol),
    as.numeric(l1_smoothing),
    as.numeric(tukey_c),
    as.numeric(min_scale),
    as.numeric(ridge),
    as.numeric(min_time_weight),
    prepared$predictor_treated,
    prepared$predictor_donors,
    prepared$predictor_weights
  )

  sc_name_native_result(fit, prepared)
}

sc_fit_one_from_prepared <- function(prepared,
                                    method = c("standard", "mm"),
                                    maxit = 1000L,
                                    mm_max_iter = 25L,
                                    tukey_c = 4.685,
                                    min_time_weight = 1e-8,
                                    tol = 1e-8,
                                    ridge = 1e-10) {
  method <- match.arg(method)
  if (identical(method, "mm")) {
    sc_fit_mm_core(
      prepared = prepared,
      maxit = maxit,
      mm_max_iter = mm_max_iter,
      tukey_c = tukey_c,
      min_time_weight = min_time_weight,
      tol = tol,
      ridge = ridge
    )
  } else {
    sc_fit_standard_core(
      prepared = prepared,
      maxit = maxit,
      tol = tol,
      ridge = ridge
    )
  }
}
