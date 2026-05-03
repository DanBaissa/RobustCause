summarize_sc_benchmark <- function(results,
                                   group_vars = c("contamination", "n_donors", "n_pre", "n_predictors", "treatment_effect", "method")) {
  if (!is.data.frame(results)) stop("`results` must be a data.frame.", call. = FALSE)
  group_vars <- intersect(group_vars, names(results))
  if (length(group_vars) == 0L) group_vars <- "method"

  split_key <- interaction(results[group_vars], drop = TRUE, sep = " | ")
  groups <- split(results, split_key)

  safe_mean <- function(x) {
    if (is.null(x)) return(NA_real_)
    mean(x, na.rm = TRUE)
  }

  out <- lapply(groups, function(x) {
    key <- x[1L, group_vars, drop = FALSE]
    data.frame(
      key,
      n = nrow(x),
      ok_rate = mean(x$status == "ok"),
      gap_rmse = safe_mean(x$gap_rmse),
      gap_bias = safe_mean(x$gap_bias),
      gap_mae = safe_mean(x$gap_mae),
      weight_rmse = safe_mean(x$weight_rmse),
      weight_l1 = safe_mean(x$weight_l1),
      pre_rmspe = safe_mean(x$pre_rmspe),
      predictor_rmse = safe_mean(x$predictor_rmse),
      mean_robust_predictor_weight = safe_mean(x$mean_robust_predictor_weight),
      min_robust_predictor_weight = safe_mean(x$min_robust_predictor_weight),
      avg_effect_hat = safe_mean(x$avg_estimated_effect),
      avg_effect_true = safe_mean(x$avg_true_effect),
      placebo_p_rmspe_ratio = safe_mean(x$placebo_p_value_rmspe_ratio),
      placebo_p_avg_gap = safe_mean(x$placebo_p_value_avg_gap),
      row.names = NULL,
      check.names = FALSE
    )
  })

  ans <- do.call(rbind, out)
  rownames(ans) <- NULL
  ans
}
