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
