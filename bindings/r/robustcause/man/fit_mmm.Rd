\name{fit_mmm}
\alias{fit_mmm}
\title{Fit a Robust Marketing Mix Model}
\usage{
fit_mmm(
  y = NULL,
  unit_ids = NULL,
  channels = NULL,
  preclean_x = NULL,
  controls = NULL,
  data = NULL,
  outcome = NULL,
  unit = NULL,
  channel_cols = NULL,
  preclean_cols = NULL,
  control_cols = NULL,
  fit_method = c("huber", "ols", "s", "mm"),
  rho = 0.85,
  increment_method = "adaptive_clip",
  hill_lambda = 0.10,
  preclean = TRUE,
  preclean_method = "mm",
  nonnegative = TRUE,
  residual_clip_multiplier = 2.5,
  clip_c = 1.5,
  adaptive_upper = 0,
  apply_hill = TRUE,
  huber_c = 1.345,
  huber_max_iter = 75L,
  s_starts = 20L,
  s_max_iter = 50L,
  mm_max_iter = 50L,
  mm_tukey_c = 4.685,
  seed = 1234L
)
}
\description{
Fits a robust multi-channel marketing mix model with optional first-step channel
cleaning, adstock carryover, Hill saturation, and a robust final regression.
}
\value{
An object of class \code{robustcause_mmm} containing the fitted coefficients,
per-channel cleaned signals, stocks, transformed regressors, and the final
design matrix used in the regression.
}
