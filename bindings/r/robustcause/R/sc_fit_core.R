sc_loss <- function(z, y, x, period_weights = NULL) {
  w <- sc_softmax(z)
  r <- as.numeric(y - x %*% w)
  if (is.null(period_weights)) {
    mean(r^2)
  } else {
    sum(period_weights * r^2) / sum(period_weights)
  }
}

sc_fit_weighted <- function(y, x, period_weights = NULL, maxit = 1000L) {
  x <- as.matrix(x)
  y <- as.numeric(y)
  n_donors <- ncol(x)
  z0 <- rep(0, n_donors)
  opt <- stats::optim(
    par = z0,
    fn = sc_loss,
    y = y,
    x = x,
    period_weights = period_weights,
    method = "BFGS",
    control = list(maxit = as.integer(maxit))
  )
  sc_softmax(opt$par)
}

sc_fit_standard_core <- function(prepared, maxit = 1000L) {
  w <- sc_fit_weighted(prepared$treated_pre, prepared$donors_pre, maxit = maxit)
  pre_resid <- as.numeric(prepared$treated_pre - prepared$donors_pre %*% w)
  post_gap <- as.numeric(prepared$treated_post - prepared$donors_post %*% w)
  list(
    weights = stats::setNames(w, prepared$donors),
    pre_residuals = stats::setNames(pre_resid, prepared$pre_periods),
    post_gaps = stats::setNames(post_gap, prepared$post_periods),
    synthetic_pre = stats::setNames(as.numeric(prepared$donors_pre %*% w), prepared$pre_periods),
    synthetic_post = stats::setNames(as.numeric(prepared$donors_post %*% w), prepared$post_periods),
    pre_rmspe = sc_rmspe(pre_resid),
    converged = TRUE,
    iterations = NA_integer_
  )
}

sc_fit_mm_core <- function(prepared,
                           maxit = 1000L,
                           mm_max_iter = 25L,
                           tukey_c = 4.685,
                           min_time_weight = 1e-8) {
  start <- sc_fit_standard_core(prepared, maxit = maxit)
  w <- as.numeric(start$weights)
  time_weights <- rep(1, length(prepared$treated_pre))

  for (iter in seq_len(as.integer(mm_max_iter))) {
    resid <- as.numeric(prepared$treated_pre - prepared$donors_pre %*% w)
    scale <- stats::median(abs(resid - stats::median(resid))) / 0.6745
    time_weights_new <- sc_tukey_weights(resid, scale, c = tukey_c, min_weight = min_time_weight)
    w_new <- sc_fit_weighted(prepared$treated_pre, prepared$donors_pre, period_weights = time_weights_new, maxit = maxit)
    if (max(abs(w_new - w)) < 1e-6) {
      w <- w_new
      time_weights <- time_weights_new
      break
    }
    w <- w_new
    time_weights <- time_weights_new
  }

  pre_resid <- as.numeric(prepared$treated_pre - prepared$donors_pre %*% w)
  post_gap <- as.numeric(prepared$treated_post - prepared$donors_post %*% w)
  list(
    weights = stats::setNames(w, prepared$donors),
    startup_weights = start$weights,
    robust_time_weights = stats::setNames(time_weights, prepared$pre_periods),
    downweighted_periods = which(time_weights < 0.5),
    pre_residuals = stats::setNames(pre_resid, prepared$pre_periods),
    post_gaps = stats::setNames(post_gap, prepared$post_periods),
    synthetic_pre = stats::setNames(as.numeric(prepared$donors_pre %*% w), prepared$pre_periods),
    synthetic_post = stats::setNames(as.numeric(prepared$donors_post %*% w), prepared$post_periods),
    pre_rmspe = sc_rmspe(pre_resid),
    standard_rmspe = start$pre_rmspe,
    converged = TRUE,
    iterations = iter
  )
}
