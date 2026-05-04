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

test_that("fit_sc supports predictor-augmented SC and robust predictor weights", {
  set.seed(789)
  n_time <- 16
  donor_1 <- rnorm(n_time)
  donor_2 <- rnorm(n_time)
  donor_3 <- rnorm(n_time)
  treated <- 0.7 * donor_1 + 0.3 * donor_2 + rnorm(n_time, sd = 0.05)
  treated[11:16] <- treated[11:16] + 1
  outcomes <- cbind(treated = treated, donor_1 = donor_1, donor_2 = donor_2, donor_3 = donor_3)

  predictors <- rbind(
    good_feature = c(treated = 1.0, donor_1 = 0.9, donor_2 = 1.1, donor_3 = -0.2),
    weird_feature = c(treated = 10.0, donor_1 = -1.0, donor_2 = 0.5, donor_3 = 1.0)
  )

  fit_standard <- fit_sc(
    outcomes,
    treated_unit = "treated",
    treatment_start = 11,
    predictors = predictors,
    predictor_lambda = 1,
    method = "standard"
  )

  fit_mm <- fit_sc(
    outcomes,
    treated_unit = "treated",
    treatment_start = 11,
    predictors = predictors,
    predictor_lambda = 1,
    method = "mm"
  )

  expect_s3_class(fit_standard, "robustcause_sc")
  expect_s3_class(fit_mm, "robustcause_sc")
  expect_true(isTRUE(fit_standard$has_predictors))
  expect_true(isTRUE(fit_mm$has_predictors))
  expect_equal(length(fit_mm$predictor_residuals), nrow(predictors))
  expect_equal(length(fit_mm$robust_predictor_weights), nrow(predictors))
  expect_true(all(is.finite(fit_mm$robust_predictor_weights)))
  expect_equal(sum(coef(fit_mm)), 1, tolerance = 1e-6)
})

test_that("fit_sc carries predictors through placebo diagnostics", {
  set.seed(321)
  n_time <- 12
  donor_1 <- rnorm(n_time)
  donor_2 <- rnorm(n_time)
  donor_3 <- rnorm(n_time)
  treated <- 0.5 * donor_1 + 0.5 * donor_2 + rnorm(n_time, sd = 0.05)
  treated[9:12] <- treated[9:12] + 1
  outcomes <- cbind(treated = treated, donor_1 = donor_1, donor_2 = donor_2, donor_3 = donor_3)
  predictors <- rbind(
    x1 = c(treated = 1.0, donor_1 = 0.8, donor_2 = 1.2, donor_3 = -0.5),
    x2 = c(treated = 0.0, donor_1 = 0.1, donor_2 = -0.1, donor_3 = 0.5)
  )

  fit <- fit_sc(
    outcomes,
    treated_unit = "treated",
    treatment_start = 9,
    predictors = predictors,
    predictor_lambda = 0.5,
    method = "mm",
    run_placebos = TRUE
  )

  expect_s3_class(fit, "robustcause_sc")
  expect_true(is.data.frame(fit$placebo$table))
  expect_true("predictor_rmse" %in% names(fit$placebo$table))
  expect_true(all(is.finite(fit$placebo$table$predictor_rmse)))
})

test_that("fit_sc reports donor diagnostics and supports donor penalties", {
  set.seed(654)
  n_time <- 18
  donor_1 <- rnorm(n_time)
  donor_2 <- rnorm(n_time)
  donor_3 <- rnorm(n_time)
  donor_4 <- rnorm(n_time)
  donor_4[1:10] <- donor_4[1:10] + 8
  treated <- 0.6 * donor_1 + 0.4 * donor_2 + rnorm(n_time, sd = 0.05)
  treated[12:18] <- treated[12:18] + 1
  outcomes <- cbind(treated = treated, donor_1 = donor_1, donor_2 = donor_2, donor_3 = donor_3, donor_4 = donor_4)

  fit_plain <- fit_sc(
    outcomes,
    treated_unit = "treated",
    treatment_start = 12,
    method = "mm"
  )

  fit_donor <- fit_sc(
    outcomes,
    treated_unit = "treated",
    treatment_start = 12,
    method = "mm",
    robust_donors = TRUE,
    donor_penalty_lambda = 1
  )

  expect_s3_class(fit_donor, "robustcause_sc")
  expect_true(is.data.frame(fit_donor$donor_diagnostics))
  expect_equal(nrow(fit_donor$donor_diagnostics), 4L)
  expect_equal(length(fit_donor$robust_donor_weights), 4L)
  expect_equal(length(fit_donor$donor_penalties), 4L)
  expect_true(all(is.finite(fit_donor$robust_donor_weights)))
  expect_true(all(is.finite(fit_donor$donor_penalties)))
  expect_true(any(fit_donor$donor_penalties > 0))
  expect_true(isTRUE(fit_donor$robust_donors))
  expect_equal(sum(coef(fit_donor)), 1, tolerance = 1e-6)
  expect_true(all(coef(fit_donor) >= 0))
  expect_false(isTRUE(all.equal(as.numeric(coef(fit_plain)), as.numeric(coef(fit_donor)), tolerance = 1e-5)))
})

test_that("donor robustness carries through placebo diagnostics", {
  set.seed(987)
  n_time <- 14
  donor_1 <- rnorm(n_time)
  donor_2 <- rnorm(n_time)
  donor_3 <- rnorm(n_time)
  donor_3[1:8] <- donor_3[1:8] + 6
  treated <- 0.5 * donor_1 + 0.5 * donor_2 + rnorm(n_time, sd = 0.05)
  treated[10:14] <- treated[10:14] + 1
  outcomes <- cbind(treated = treated, donor_1 = donor_1, donor_2 = donor_2, donor_3 = donor_3)

  fit <- fit_sc(
    outcomes,
    treated_unit = "treated",
    treatment_start = 10,
    method = "mm",
    robust_donors = TRUE,
    donor_penalty_lambda = 1,
    run_placebos = TRUE
  )

  expect_s3_class(fit, "robustcause_sc")
  expect_true(is.data.frame(fit$placebo$table))
  expect_true("min_robust_donor_weight" %in% names(fit$placebo$table))
  expect_true(all(is.finite(fit$placebo$table$min_robust_donor_weight)))
})
