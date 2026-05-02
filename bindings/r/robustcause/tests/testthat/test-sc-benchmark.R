test_that("SC simulation helper returns usable panel and truth", {
  sim <- simulate_sc_panel(
    n_donors = 5,
    n_pre = 8,
    n_post = 3,
    treatment_effect = 1,
    contamination = "treated_pre_spike",
    contamination_rate = 0.2,
    seed = 123
  )

  expect_s3_class(sim, "robustcause_sc_simulation")
  expect_equal(ncol(sim$outcomes), 6L)
  expect_equal(nrow(sim$outcomes), 11L)
  expect_equal(nrow(sim$panel), 66L)
  expect_equal(length(sim$truth$true_weights), 5L)
  expect_true(abs(sum(sim$truth$true_weights) - 1) < 1e-8)
  expect_true(nrow(sim$contamination_info) > 0L)
})

test_that("SC benchmark runner compares standard and MM fits", {
  result <- run_sc_benchmark_once(
    seed = 321,
    methods = c("standard", "mm"),
    run_placebos = FALSE,
    n_donors = 5,
    n_pre = 8,
    n_post = 3,
    n_factors = 1,
    true_sparsity = 2,
    treatment_effect = 1,
    effect_shape = "constant",
    contamination = "treated_pre_spike",
    contamination_rate = 0.2,
    contamination_magnitude = 4,
    maxit = 100,
    subproblem_maxit = 100
  )

  expect_equal(nrow(result), 2L)
  expect_equal(sort(result$method), c("mm", "standard"))
  expect_true(all(result$status == "ok"))
  expect_true(all(is.finite(result$gap_rmse)))
  expect_true(all(is.finite(result$weight_rmse)))
})
