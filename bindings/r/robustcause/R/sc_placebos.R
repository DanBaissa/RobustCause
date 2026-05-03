sc_run_placebos <- function(prepared,
                            method = c("standard", "mm"),
                            maxit = 1000L,
                            mm_max_iter = 25L,
                            tukey_c = 4.685,
                            min_time_weight = 1e-8) {
  method <- match.arg(method)
  outcomes <- prepared$outcomes
  n_units <- ncol(outcomes)
  unit_names <- colnames(outcomes)

  rows <- vector("list", n_units)
  post_gap_matrix <- matrix(NA_real_, nrow = n_units, ncol = length(prepared$post_periods))
  rownames(post_gap_matrix) <- unit_names
  colnames(post_gap_matrix) <- as.character(prepared$post_periods)

  for (j in seq_len(n_units)) {
    donor_cols <- setdiff(seq_len(n_units), j)
    placebo_predictors <- if (length(prepared$predictor_names) == 0L) NULL else prepared$predictors
    placebo_prepared <- sc_prepare_matrix(
      outcomes = outcomes,
      treated_unit = j,
      treatment_start = prepared$treatment_start,
      donors = donor_cols,
      predictors = placebo_predictors,
      predictor_weights = prepared$predictor_weights,
      predictor_lambda = prepared$predictor_lambda,
      predictor_scale = "none",
      na.action = "fail"
    )
    fit <- sc_fit_one_from_prepared(
      prepared = placebo_prepared,
      method = method,
      maxit = maxit,
      mm_max_iter = mm_max_iter,
      tukey_c = tukey_c,
      min_time_weight = min_time_weight
    )
    post_gap_matrix[j, ] <- as.numeric(fit$post_gaps)
    post_rmspe <- sc_rmspe(fit$post_gaps)
    pre_rmspe <- max(sc_rmspe(fit$pre_residuals), .Machine$double.eps)
    rows[[j]] <- data.frame(
      unit = unit_names[[j]],
      is_treated = identical(unit_names[[j]], prepared$treated_unit),
      pre_rmspe = pre_rmspe,
      post_rmspe = post_rmspe,
      rmspe_ratio = post_rmspe / pre_rmspe,
      avg_abs_post_gap = mean(abs(fit$post_gaps)),
      predictor_rmse = fit$predictor_rmse,
      stringsAsFactors = FALSE
    )
  }

  table <- do.call(rbind, rows)
  treated_idx <- match(prepared$treated_unit, table$unit)
  treated_ratio <- table$rmspe_ratio[[treated_idx]]
  treated_avg_abs_gap <- table$avg_abs_post_gap[[treated_idx]]

  list(
    table = table,
    post_gap_matrix = post_gap_matrix,
    treated_rmspe_ratio = treated_ratio,
    placebo_p_value_rmspe_ratio = mean(table$rmspe_ratio >= treated_ratio),
    treated_avg_abs_post_gap = treated_avg_abs_gap,
    placebo_p_value_avg_gap = mean(table$avg_abs_post_gap >= treated_avg_abs_gap)
  )
}
