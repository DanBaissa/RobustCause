summarize_sc_benchmark <- function(results,
                                   group_vars = c("contamination", "n_donors", "n_pre", "treatment_effect", "method")) {
  if (!is.data.frame(results)) stop("`results` must be a data.frame.", call. = FALSE)
  group_vars <- intersect(group_vars, names(results))
  if (length(group_vars) == 0L) group_vars <- "method"

  split_key <- interaction(results[group_vars], drop = TRUE, sep = " | ")
  groups <- split(results, split_key)

  out <- lapply(groups, function(x) {
    key <- x[1L, group_vars, drop = FALSE]
    data.frame(
      key,
      n = nrow(x),
      ok_rate = mean(x$status == "ok"),
      gap_rmse = mean(x$gap_rmse, na.rm = TRUE),
      gap_bias = mean(x$gap_bias, na.rm = TRUE),
      gap_mae = mean(x$gap_mae, na.rm = TRUE),
      weight_rmse = mean(x$weight_rmse, na.rm = TRUE),
      weight_l1 = mean(x$weight_l1, na.rm = TRUE),
      pre_rmspe = mean(x$pre_rmspe, na.rm = TRUE),
      avg_effect_hat = mean(x$avg_estimated_effect, na.rm = TRUE),
      avg_effect_true = mean(x$avg_true_effect, na.rm = TRUE),
      placebo_p_rmspe_ratio = mean(x$placebo_p_value_rmspe_ratio, na.rm = TRUE),
      placebo_p_avg_gap = mean(x$placebo_p_value_avg_gap, na.rm = TRUE),
      row.names = NULL,
      check.names = FALSE
    )
  })

  ans <- do.call(rbind, out)
  rownames(ans) <- NULL
  ans
}
