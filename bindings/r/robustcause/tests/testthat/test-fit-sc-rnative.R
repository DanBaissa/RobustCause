test_that("R-native fit_sc fits standard and MM synthetic controls", {
  set.seed(123)
  n_time <- 12
  donor_1 <- rnorm(n_time)
  donor_2 <- rnorm(n_time)
  donor_3 <- rnorm(n_time)
  treated <- 0.6 * donor_1 + 0.4 * donor_2 + rnorm(n_time, sd = 0.05)
  treated[9:12] <- treated[9:12] + 1
  outcomes <- cbind(treated = treated, donor_1 = donor_1, donor_2 = donor_2, donor_3 = donor_3)

  fit_standard <- fit_sc(outcomes, treated_unit = "treated", treatment_start = 9, method = "standard")
  fit_mm <- fit_sc(outcomes, treated_unit = "treated", treatment_start = 9, method = "mm")

  expect_s3_class(fit_standard, "robustcause_sc")
  expect_s3_class(fit_mm, "robustcause_sc")
  expect_equal(sum(coef(fit_standard)), 1, tolerance = 1e-6)
  expect_equal(sum(coef(fit_mm)), 1, tolerance = 1e-6)
  expect_true(all(coef(fit_standard) >= 0))
  expect_true(all(coef(fit_mm) >= 0))
  expect_equal(length(fit_standard$post_gaps), 4L)
  expect_equal(nrow(as.data.frame(fit_standard)), n_time)
})

test_that("R-native fit_sc returns placebo diagnostics", {
  set.seed(456)
  n_time <- 14
  donor_1 <- rnorm(n_time)
  donor_2 <- rnorm(n_time)
  donor_3 <- rnorm(n_time)
  treated <- 0.5 * donor_1 + 0.5 * donor_2 + rnorm(n_time, sd = 0.05)
  treated[10:14] <- treated[10:14] + 1
  outcomes <- cbind(treated = treated, donor_1 = donor_1, donor_2 = donor_2, donor_3 = donor_3)

  fit <- fit_sc(
    outcomes,
    treated_unit = "treated",
    treatment_start = 10,
    method = "mm",
    run_placebos = TRUE
  )

  expect_s3_class(fit, "robustcause_sc")
  expect_true(is.list(fit$placebo))
  expect_true(is.data.frame(fit$placebo$table))
  expect_equal(nrow(fit$placebo$table), ncol(outcomes))
  expect_equal(nrow(fit$placebo$post_gap_matrix), ncol(outcomes))
  expect_true(is.finite(fit$inference$placebo_p_value_rmspe_ratio))
  expect_true(is.finite(fit$inference$placebo_p_value_avg_gap))
})
