fit_sc <- function(outcomes = NULL,
                   treated_unit = NULL,
                   treatment_start = NULL,
                   donors = NULL,
                   predictors = NULL,
                   predictor_weights = NULL,
                   predictor_lambda = 1,
                   predictor_scale = c("mad", "sd", "none"),
                   method = c("standard", "mm"),
                   robust_donors = FALSE,
                   donor_penalty_lambda = 1,
                   donor_tukey_c = 4.685,
                   min_donor_weight = 1e-4,
                   max_donor_penalty = 1e4,
                   run_placebos = FALSE,
                   maxit = 1000L,
                   mm_max_iter = 25L,
                   tukey_c = 4.685,
                   min_time_weight = 1e-8,
                   na.action = c("omit", "fail")) {
  method <- match.arg(method)
  predictor_scale <- match.arg(predictor_scale)
  na.action <- match.arg(na.action)

  prepared <- sc_prepare_matrix(
    outcomes = outcomes,
    treated_unit = treated_unit,
    treatment_start = treatment_start,
    donors = donors,
    predictors = predictors,
    predictor_weights = predictor_weights,
    predictor_lambda = predictor_lambda,
    predictor_scale = predictor_scale,
    na.action = na.action
  )

  fit <- sc_fit_one_from_prepared(
    prepared = prepared,
    method = method,
    maxit = maxit,
    mm_max_iter = mm_max_iter,
    tukey_c = tukey_c,
    min_time_weight = min_time_weight,
    robust_donors = robust_donors,
    donor_penalty_lambda = donor_penalty_lambda,
    donor_tukey_c = donor_tukey_c,
    min_donor_weight = min_donor_weight,
    max_donor_penalty = max_donor_penalty
  )

  fit$method <- method
  fit$robust_donors <- isTRUE(robust_donors)
  fit$donor_penalty_lambda <- donor_penalty_lambda
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
      min_time_weight = min_time_weight,
      robust_donors = robust_donors,
      donor_penalty_lambda = donor_penalty_lambda,
      donor_tukey_c = donor_tukey_c,
      min_donor_weight = min_donor_weight,
      max_donor_penalty = max_donor_penalty
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
    robust_donors = object$robust_donors,
    donor_penalty_lambda = object$donor_penalty_lambda,
    treated_unit = object$treated_unit,
    n_donors = length(object$donors),
    pre_rmspe = object$pre_rmspe,
    predictor_rmse = object$predictor_rmse,
    post_mean_gap = mean(object$post_gaps),
    post_mean_abs_gap = mean(abs(object$post_gaps)),
    effective_donors = object$effective_donors,
    max_weight = object$max_weight,
    converged = object$converged,
    iterations = object$iterations,
    weights = object$weights,
    predictor_residuals = object$predictor_residuals,
    donor_diagnostics = object$donor_diagnostics,
    donor_scores = object$donor_scores,
    robust_donor_weights = object$robust_donor_weights,
    donor_penalties = object$donor_penalties,
    inference = object$inference
  )
  if (!is.null(object$robust_time_weights)) {
    out$mean_robust_time_weight <- mean(object$robust_time_weights)
    out$downweighted_periods <- object$downweighted_periods
  }
  if (!is.null(object$robust_predictor_weights) && length(object$robust_predictor_weights) > 0L) {
    out$mean_robust_predictor_weight <- mean(object$robust_predictor_weights)
    out$downweighted_predictors <- names(object$robust_predictor_weights)[object$robust_predictor_weights < 0.5]
    out$robust_predictor_weights <- object$robust_predictor_weights
  }
  if (!is.null(object$robust_donor_weights) && length(object$robust_donor_weights) > 0L) {
    out$mean_robust_donor_weight <- mean(object$robust_donor_weights)
    out$min_robust_donor_weight <- min(object$robust_donor_weights)
    out$downweighted_donors <- names(object$robust_donor_weights)[object$robust_donor_weights < 0.5]
  }
  class(out) <- "summary.robustcause_sc"
  out
}

print.summary.robustcause_sc <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat(sprintf("RobustCause synthetic control | Method: %s | Treated unit: %s\n\n", x$method, x$treated_unit))
  cat(sprintf("Pre-RMSPE: %s\n", format(signif(x$pre_rmspe, digits))))
  if (is.finite(x$predictor_rmse)) {
    cat(sprintf("Predictor RMSE: %s\n", format(signif(x$predictor_rmse, digits))))
  }
  cat(sprintf("Post mean gap: %s\n", format(signif(x$post_mean_gap, digits))))
  cat(sprintf("Post mean abs gap: %s\n", format(signif(x$post_mean_abs_gap, digits))))
  cat(sprintf("Effective donors: %s\n", format(signif(x$effective_donors, digits))))
  cat(sprintf("Max donor weight: %s\n", format(signif(x$max_weight, digits))))
  if (isTRUE(x$robust_donors)) {
    cat(sprintf("Robust donor penalties: on | lambda: %s\n", format(signif(x$donor_penalty_lambda, digits))))
  }
  cat("\nWeights:\n")
  print(round(x$weights, digits))
  if (!is.null(x$robust_predictor_weights) && length(x$robust_predictor_weights) > 0L) {
    cat("\nRobust predictor weights:\n")
    print(round(x$robust_predictor_weights, digits))
  }
  if (!is.null(x$robust_donor_weights) && length(x$robust_donor_weights) > 0L) {
    cat("\nLowest robust donor weights:\n")
    donor_order <- order(x$robust_donor_weights)
    print(round(x$robust_donor_weights[head(donor_order, min(5L, length(donor_order)))], digits))
  }
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
