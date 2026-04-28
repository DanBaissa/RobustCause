# robustcause

`robustcause` is an R interface to the RobustCause C++ core for robust causal
estimation and marketing measurement.

Current package coverage includes:

- robust linear regression with HC sandwich inference
- Tukey-bisquare S-estimation
- standard synthetic control
- robust MM synthetic control with placebo summaries
- robust MM-DML
- adstock construction and robust MMM helpers

## Synthetic control quick start

```r
library(robustcause)

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
  method = "mm",
  run_placebos = TRUE
)

summary(fit)
coef(fit)
head(as.data.frame(fit))
plot(fit)
plot(fit, type = "weights")
plot(fit, type = "placebo_gaps")
```

The MM synthetic-control estimator robustifies pre-treatment fitting with a
simplex-constrained startup, a MAD residual scale, and Tukey-bisquare
refinement. It reduces sensitivity to contaminated pre-periods, but it does
not solve donor-pool mismatch or lack of overlap.

## Current status

The synthetic-control frontend is usable on real data today for standard
synthetic control and robust MM synthetic control, using either matrix input or
panel-data input, with optional predictor balancing and placebo runs.

Current limits:

- no augmented or ridge synthetic control yet
- no matrix-completion variants
- no richer formal inference beyond placebo-style summaries
- no native solver-side missing-data support after preprocessing
- no dedicated QP backend; the simplex subproblem still uses projected gradient

## Package notes

The package compiles vendored copies of the current core sources so it can be
built as a self-contained R package inside this repository.
