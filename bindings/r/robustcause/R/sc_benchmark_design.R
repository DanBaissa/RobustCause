make_sc_benchmark_design <- function(contamination = c("none", "treated_pre_spike", "donor_pre_spike", "treated_and_donor_pre_spike", "heavy_tails", "donor_pool_mismatch", "predictor_spike", "mixed_outcome_predictor", "bad_donor_pre_path", "bad_donor_predictors", "bad_donor_mixed"),
                                     n_donors = c(10L, 25L),
                                     n_pre = c(15L, 30L, 60L),
                                     n_post = 10L,
                                     n_factors = 2L,
                                     true_sparsity = 5L,
                                     treatment_effect = c(0, 1),
                                     effect_shape = c("constant", "ramp"),
                                     contamination_rate = c(0, 0.05, 0.10),
                                     contamination_magnitude = c(4, 8),
                                     donor_mismatch_strength = c(0, 1),
                                     noise_sd = 1,
                                     n_predictors = c(0L, 8L),
                                     predictor_noise_sd = 0.25,
                                     predictor_contamination_rate = c(0, 0.10),
                                     predictor_contamination_magnitude = c(4, 8),
                                     bad_donor_rate = c(0.10, 0.20)) {
  design <- expand.grid(
    contamination = contamination,
    n_donors = n_donors,
    n_pre = n_pre,
    n_post = n_post,
    n_factors = n_factors,
    true_sparsity = true_sparsity,
    treatment_effect = treatment_effect,
    effect_shape = effect_shape,
    contamination_rate = contamination_rate,
    contamination_magnitude = contamination_magnitude,
    donor_mismatch_strength = donor_mismatch_strength,
    noise_sd = noise_sd,
    n_predictors = n_predictors,
    predictor_noise_sd = predictor_noise_sd,
    predictor_contamination_rate = predictor_contamination_rate,
    predictor_contamination_magnitude = predictor_contamination_magnitude,
    bad_donor_rate = bad_donor_rate,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  min_contam_mag <- contamination_magnitude[[1]]
  min_pred_contam_mag <- predictor_contamination_magnitude[[1]]
  min_bad_donor_rate <- bad_donor_rate[[1]]

  keep_none <- design$contamination == "none" &
    design$contamination_rate == 0 &
    design$predictor_contamination_rate == 0 &
    design$contamination_magnitude == min_contam_mag &
    design$predictor_contamination_magnitude == min_pred_contam_mag &
    design$bad_donor_rate == min_bad_donor_rate &
    design$donor_mismatch_strength == 0

  keep_mismatch <- design$contamination == "donor_pool_mismatch" &
    design$contamination_rate == 0 &
    design$predictor_contamination_rate == 0 &
    design$contamination_magnitude == min_contam_mag &
    design$predictor_contamination_magnitude == min_pred_contam_mag &
    design$bad_donor_rate == min_bad_donor_rate &
    design$donor_mismatch_strength > 0

  keep_outcome <- design$contamination %in% c("treated_pre_spike", "donor_pre_spike", "treated_and_donor_pre_spike", "heavy_tails") &
    design$contamination_rate > 0 &
    design$predictor_contamination_rate == 0 &
    design$bad_donor_rate == min_bad_donor_rate &
    design$donor_mismatch_strength == 0

  keep_predictor <- design$contamination == "predictor_spike" &
    design$n_predictors > 0 &
    design$contamination_rate == 0 &
    design$predictor_contamination_rate > 0 &
    design$contamination_magnitude == min_contam_mag &
    design$donor_mismatch_strength == 0

  keep_mixed <- design$contamination == "mixed_outcome_predictor" &
    design$n_predictors > 0 &
    design$contamination_rate > 0 &
    design$predictor_contamination_rate > 0 &
    design$donor_mismatch_strength == 0

  keep_bad_donor_pre <- design$contamination == "bad_donor_pre_path" &
    design$contamination_rate == 0 &
    design$predictor_contamination_rate == 0 &
    design$predictor_contamination_magnitude == min_pred_contam_mag &
    design$donor_mismatch_strength == 0

  keep_bad_donor_predictors <- design$contamination == "bad_donor_predictors" &
    design$n_predictors > 0 &
    design$contamination_rate == 0 &
    design$predictor_contamination_rate == 0 &
    design$contamination_magnitude == min_contam_mag &
    design$donor_mismatch_strength == 0

  keep_bad_donor_mixed <- design$contamination == "bad_donor_mixed" &
    design$n_predictors > 0 &
    design$contamination_rate == 0 &
    design$predictor_contamination_rate == 0 &
    design$donor_mismatch_strength == 0

  design <- design[
    keep_none |
      keep_mismatch |
      keep_outcome |
      keep_predictor |
      keep_mixed |
      keep_bad_donor_pre |
      keep_bad_donor_predictors |
      keep_bad_donor_mixed,
    ,
    drop = FALSE
  ]
  rownames(design) <- NULL
  design$design_id <- seq_len(nrow(design))
  design
}
