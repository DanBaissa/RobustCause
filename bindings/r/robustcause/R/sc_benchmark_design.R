make_sc_benchmark_design <- function(contamination = c("none", "treated_pre_spike", "donor_pre_spike", "treated_and_donor_pre_spike", "heavy_tails", "donor_pool_mismatch"),
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
                                     noise_sd = 1) {
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
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  keep_none <- design$contamination == "none" &
    design$contamination_rate == 0 &
    design$contamination_magnitude == contamination_magnitude[[1]] &
    design$donor_mismatch_strength == 0

  keep_mismatch <- design$contamination == "donor_pool_mismatch" &
    design$contamination_rate == 0 &
    design$contamination_magnitude == contamination_magnitude[[1]] &
    design$donor_mismatch_strength > 0

  keep_other <- !(design$contamination %in% c("none", "donor_pool_mismatch")) &
    design$contamination_rate > 0 &
    design$donor_mismatch_strength == 0

  design <- design[keep_none | keep_mismatch | keep_other, , drop = FALSE]
  rownames(design) <- NULL
  design$design_id <- seq_len(nrow(design))
  design
}
