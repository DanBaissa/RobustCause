sc_name_native_result <- function(fit, prepared) {
  fit$weights <- stats::setNames(as.numeric(fit$weights), prepared$donors)
  fit$pre_residuals <- stats::setNames(as.numeric(fit$pre_residuals), prepared$pre_periods)
  fit$post_gaps <- stats::setNames(as.numeric(fit$post_gaps), prepared$post_periods)
  fit$synthetic_pre <- stats::setNames(as.numeric(fit$synthetic_pre), prepared$pre_periods)
  fit$synthetic_post <- stats::setNames(as.numeric(fit$synthetic_post), prepared$post_periods)

  if (!is.null(fit$startup_weights)) {
    fit$startup_weights <- stats::setNames(as.numeric(fit$startup_weights), prepared$donors)
  }
  if (!is.null(fit$robust_time_weights)) {
    fit$robust_time_weights <- stats::setNames(as.numeric(fit$robust_time_weights), prepared$pre_periods)
  }

  fit
}

sc_fit_standard_core <- function(prepared,
                                 maxit = 1000L,
                                 tol = 1e-8,
                                 ridge = 1e-10) {
  fit <- .Call(
    rc_r_fit_sc,
    prepared$treated_pre,
    prepared$donors_pre,
    prepared$treated_post,
    prepared$donors_post,
    as.integer(maxit),
    as.numeric(tol),
    as.numeric(ridge),
    numeric(),
    matrix(numeric(), nrow = 0L, ncol = ncol(prepared$donors_pre)),
    numeric()
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
  fit <- .Call(
    rc_r_fit_mm_sc,
    prepared$treated_pre,
    prepared$donors_pre,
    prepared$treated_post,
    prepared$donors_post,
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
    numeric(),
    matrix(numeric(), nrow = 0L, ncol = ncol(prepared$donors_pre)),
    numeric()
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
