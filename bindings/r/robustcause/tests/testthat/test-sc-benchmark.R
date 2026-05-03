test_that("SC benchmark helpers run on a tiny design", {
  design <- make_sc_benchmark_design(
    contamination = c("none", "treated_pre_spike"),
    n_donors = 4,
    n_pre = 6,
    n_post = 3,
    n_factors = 1,
    true_sparsity = 2,
    treatment_effect = 1,
    effect_shape = "constant",
    contamination_rate = c(0, 0.2),
    contamination_magnitude = 4,
    donor_mismatch_strength = 0
  )

  results <- run_sc_benchmark_grid(
    design = design,
    reps = 1,
    methods = c("standard", "mm"),
    run_placebos = FALSE,
    n_cores = 1,
    seed = 123
  )

  summary <- summarize_sc_benchmark(results)

  expect_true(nrow(design) >= 2L)
  expect_true(nrow(results) >= 4L)
  expect_true(nrow(summary) >= 2L)
  expect_true("gap_rmse" %in% names(summary))
  expect_true(all(results$status == "ok"))
})
