fit_sc <- function(outcomes = NULL,
                   treated_unit = NULL,
                   treatment_start = NULL,
                   donors = NULL,
                   method = c("standard", "mm"),
                   run_placebos = FALSE,
                   maxit = 1000L,
                   mm_max_iter = 25L,
                   tukey_c = 4.685,
                   min_time_weight = 1e-8,
                   na.action = c("omit", "fail")) {
  method <- match.arg(method)
  na.action <- match.arg(na.action)

  prepared <- sc_prepare_matrix(
    outcomes = outcomes,
    treated_unit = treated_unit,
    treatment_start = treatment_start,
    donors = donors,
    na.action = na.action
  )

  fit <- sc_fit_one_from_prepared(
    prepared = prepared,
    method = method,
    maxit = maxit,
    mm_max_iter = mm_max_iter,
    tukey_c = tukey_c,
    min_time_weight = min_time_weight
  )

  fit$method <- method
  fit$treated_unit <- prepared$treated_unit
  fit$donors <- prepared$donors
  fit$pre_periods <- prepared$pre_periods
  fit$post_periods <- prepared$post_periods
  fit$call <- match.call()
  fit$effective_donors <- 1 / sum(as.numeric(fit$weights)^2)
  fit$max_weight <- max(as.numeric(fit$weights))
  fit$weight_herfindahl <- sum(as.numeric(fit$weights)^2)
  fit$path <- data.frame(
    time = c(prepared$pre_periods, prepared$post_periods),
    period = c(rep("pre", length(prepared$pre_periods)), rep("post", length(prepared$post_periods))),
    treated = c(prepared$treated_pre, prepared$treated_post),
    synthetic = c(as.numeric(fit$synthetic_pre), as.numeric(fit$synthetic_post)),
    gap = c(as.numeric(fit$pre_residuals), as.numeric(fit$post_gaps)),
    row.names = NULL
  )

  if (isTRUE(run_placebos)) {
    fit$placebo <- sc_run_placebos(
      prepared = prepared,
      method = method,
      maxit = maxit,
      mm_max_iter = mm_max_iter,
      tukey_c = tukey_c,
      min_time_weight = min_time_weight
    )
    fit$inference <- list(
      placebo_p_value_rmspe_ratio = fit$placebo$placebo_p_value_rmspe_ratio,
      placebo_p_value_avg_gap = fit$placebo$placebo_p_value_avg_gap,
      treated_rmspe_ratio = fit$placebo$treated_rmspe_ratio,
      treated_avg_abs_post_gap = fit$placebo$treated_avg_abs_post_gap
    )
  } else {
    fit$placebo <- NULL
    fit$inference <- NULL
  }

  class(fit) <- "robustcause_sc"
  fit
}

summary.robustcause_sc <- function(object, ...) {
  out <- list(
    call = object$call,
    method = object$method,
    treated_unit = object$treated_unit,
    n_donors = length(object$donors),
    pre_rmspe = object$pre_rmspe,
    post_mean_gap = mean(object$post_gaps),
    post_mean_abs_gap = mean(abs(object$post_gaps)),
    effective_donors = object$effective_donors,
    max_weight = object$max_weight,
    converged = object$converged,
    iterations = object$iterations,
    weights = object$weights,
    inference = object$inference
  )
  if (!is.null(object$robust_time_weights)) {
    out$mean_robust_time_weight <- mean(object$robust_time_weights)
    out$downweighted_periods <- object$downweighted_periods
  }
  class(out) <- "summary.robustcause_sc"
  out
}

print.summary.robustcause_sc <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat(sprintf("RobustCause synthetic control | Method: %s | Treated unit: %s\n\n", x$method, x$treated_unit))
  cat(sprintf("Pre-RMSPE: %s\n", format(signif(x$pre_rmspe, digits))))
  cat(sprintf("Post mean gap: %s\n", format(signif(x$post_mean_gap, digits))))
  cat(sprintf("Post mean abs gap: %s\n", format(signif(x$post_mean_abs_gap, digits))))
  cat(sprintf("Effective donors: %s\n", format(signif(x$effective_donors, digits))))
  cat(sprintf("Max donor weight: %s\n", format(signif(x$max_weight, digits))))
  cat("\nWeights:\n")
  print(round(x$weights, digits))
  if (!is.null(x$inference)) {
    cat("\nPlacebo diagnostics:\n")
    cat(sprintf("  RMSPE ratio p-value: %s\n", format(signif(x$inference$placebo_p_value_rmspe_ratio, digits))))
    cat(sprintf("  Avg |gap| p-value: %s\n", format(signif(x$inference$placebo_p_value_avg_gap, digits))))
  }
  invisible(x)
}

print.robustcause_sc <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

coef.robustcause_sc <- function(object, ...) object$weights
fitted.robustcause_sc <- function(object, ...) stats::setNames(object$path$synthetic, object$path$time)
residuals.robustcause_sc <- function(object, ...) stats::setNames(object$path$gap, object$path$time)
as.data.frame.robustcause_sc <- function(x, ...) x$path
