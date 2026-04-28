test_that("standard synthetic control works on matrix input", {
  outcomes <- cbind(
    treated = c(1, 2, 3, 10, 11, 12),
    d1 = c(1, 2, 2.8, 3, 3.2, 3.3),
    d2 = c(0.9, 1.9, 3.1, 3.5, 3.9, 4.1),
    d3 = c(1.2, 2.1, 2.9, 3.2, 3.4, 3.6)
  )

  fit <- fit_sc(
    outcomes = outcomes,
    treated_unit = "treated",
    treatment_start = 4,
    method = "standard",
    run_placebos = FALSE
  )

  expect_s3_class(fit, "robustcause_sc")
  expect_identical(fit$method, "standard")
  expect_equal(length(fit$weights), 3L)
  expect_true(is.numeric(fit$pre_rmspe))
  expect_true(abs(sum(as.numeric(fit$weights)) - 1) < 1e-6)
  expect_equal(names(coef(fit)), c("d1", "d2", "d3"))
  expect_equal(nrow(fit$path), 6L)
  expect_equal(colnames(fit$path), c("time", "period", "treated", "synthetic", "gap"))
})

test_that("synthetic control supports panel input and predictors", {
  panel <- data.frame(
    unit = rep(c("treated", "d1", "d2", "d3"), each = 6),
    time = rep(1:6, times = 4),
    outcome = c(
      1, 2, 3, 10, 11, 12,
      1, 2, 2.8, 3, 3.2, 3.3,
      0.9, 1.9, 3.1, 3.5, 3.9, 4.1,
      1.2, 2.1, 2.9, 3.2, 3.4, 3.6
    ),
    predictor_a = c(
      5, 5, 5, 9, 9, 9,
      5, 5, 5, 6, 6, 6,
      4, 4, 4, 7, 7, 7,
      6, 6, 6, 8, 8, 8
    ),
    predictor_b = c(
      2, 2, 2, 4, 4, 4,
      2, 2, 2, 3, 3, 3,
      1, 1, 1, 3, 3, 3,
      2, 2, 2, 3, 3, 3
    )
  )

  fit <- fit_sc(
    data = panel,
    unit = "unit",
    time = "time",
    outcome = "outcome",
    treated_unit = "treated",
    treatment_start = 4,
    predictors = c("predictor_a", "predictor_b"),
    predictor_weights = c(2, 1),
    method = "standard",
    run_placebos = FALSE
  )

  expect_equal(fit$predictors, c("predictor_a", "predictor_b"))
  expect_equal(length(fit$donors), 3L)
  expect_named(fit$weights, c("d1", "d2", "d3"))
})

test_that("mm synthetic control exposes robust diagnostics and placebo inference", {
  outcomes <- cbind(
    treated = c(1, 2, 20, 10, 11, 12),
    d1 = c(1, 2, 2.8, 3, 3.2, 3.3),
    d2 = c(0.9, 1.9, 3.1, 3.5, 3.9, 4.1),
    d3 = c(1.2, 2.1, 2.9, 3.2, 3.4, 3.6)
  )

  fit <- fit_sc(
    outcomes = outcomes,
    treated_unit = "treated",
    treatment_start = 4,
    method = "mm",
    run_placebos = TRUE
  )

  expect_identical(fit$method, "mm")
  expect_true(is.numeric(fit$pre_rmspe))
  expect_true(is.numeric(fit$mm_rmspe))
  expect_equal(length(fit$robust_time_weights), 3L)
  expect_false(is.null(fit$placebo))
  expect_false(is.null(fit$placebos))
  expect_true(is.matrix(fit$placebo$post_gap_matrix))
  expect_true(is.list(fit$inference))
  expect_true(is.numeric(fit$inference$placebo_p_value_rmspe_ratio))
  expect_equal(nrow(fit$placebo$post_gap_matrix), 4L)
})

test_that("synthetic control methods return stable user-facing objects", {
  outcomes <- cbind(
    treated = c(1, 2, 3, 10, 11, 12),
    d1 = c(1, 2, 2.8, 3, 3.2, 3.3),
    d2 = c(0.9, 1.9, 3.1, 3.5, 3.9, 4.1),
    d3 = c(1.2, 2.1, 2.9, 3.2, 3.4, 3.6)
  )

  fit <- fit_sc(
    outcomes = outcomes,
    treated_unit = "treated",
    treatment_start = 4,
    method = "standard",
    run_placebos = TRUE
  )

  expect_equal(length(fitted(fit)), 6L)
  expect_equal(length(residuals(fit)), 6L)
  expect_equal(length(predict(fit, type = "gap")), 6L)
  expect_equal(length(predict(fit, type = "synthetic", period = "post")), 3L)
  expect_equal(nrow(as.data.frame(fit)), 6L)
  expect_invisible(plot(fit))
  expect_invisible(plot(fit, type = "weights"))
  expect_invisible(plot(fit, type = "placebo_gaps"))
})
