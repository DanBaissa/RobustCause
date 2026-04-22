# RobustCause

RobustCause is a C++ core library for robust causal inference, robust marketing-mix modeling, and robust regression, with an R package already built on top of it and a Python binding layer intended next.

The project is organized around a shared numerical backend:

- robust linear regression with M and MM estimation
- Tukey bisquare S-estimation
- heteroskedasticity-consistent covariance estimators for robust and classical linear models
- MM double machine learning for partially linear treatment effects
- robust adstock construction with optional S/MM precleaning
- robust marketing mix modeling with per-channel cleaning, carryover, Hill saturation, and robust final regression

## What The API Does

Today, the most complete user-facing API is the R package in [bindings/r/robustcause](/A:/RobustCause/bindings/r/robustcause:1).

That package exposes a consistent set of fit objects with `print()` and `summary()` methods for:

- `fit_rlm()` for robust linear models
- `fit_s_estimator()` for standalone S-estimation
- `fit_mm_dml()` for robust double machine learning
- `build_adstock()` for robust carryover stock construction
- `fit_mmm()` for multi-channel marketing mix models

The underlying C++ core lives in:

- [include/robust](/A:/RobustCause/include/robust:1)
- [include/robustcause](/A:/RobustCause/include/robustcause:1)
- [src](/A:/RobustCause/src:1)

That core is intended to be the shared engine used by both R and future Python bindings.

## Current Status

Implemented now:

- installable C++ core library with CMake packaging
- shared C++ modules for robust regression, S-estimation, adstock, and MMM
- R package wrapper with compiled native code
- robust inference helpers for both RobustCause fits and base R `lm` objects
- smoke-test/showcase notebook at [bindings/r/robustcause/robustcause-smoke-test.Rmd](/A:/RobustCause/bindings/r/robustcause/robustcause-smoke-test.Rmd:1)

Not implemented yet:

- Python bindings
- a polished end-user Python package
- a unified high-level frontend across R and Python

## Repository Layout

- [include/robust](/A:/RobustCause/include/robust:1): core C++ headers
- [include/robustcause](/A:/RobustCause/include/robustcause:1): umbrella headers and C ABI
- [src](/A:/RobustCause/src:1): core C++ implementations
- [examples](/A:/RobustCause/examples:1): C++ example programs
- [tests](/A:/RobustCause/tests:1): smoke tests for the C++ core
- [bindings/r/robustcause](/A:/RobustCause/bindings/r/robustcause:1): R package
- [docs](/A:/RobustCause/docs:1): architecture notes

## R Package Installation

Install from GitHub:

```r
install.packages(c("remotes", "Rcpp", "RcppArmadillo", "RcppEigen"))
remotes::install_github("DanBaissa/RobustCause", subdir = "bindings/r/robustcause")
```

Then load it:

```r
library(robustcause)
```

## Core R Functions

### `fit_rlm()`

Fits a robust linear model using either:

- `method = "m"` for standard M-estimation
- `method = "mm"` for MM-estimation initialized from the S-estimator

Key behavior:

- accepts either matrix input or formula/data input
- defaults to Huber psi for M-estimation
- defaults to Tukey bisquare psi for MM-estimation
- includes native support for:
  - observation `weights`
  - `wt.method = "case"` and `"inv.var"`
  - built-in psi families `"huber"`, `"tukey_bisquare"`, and `"hampel"`
  - `init = "ls"`, `init = "lts"`, and user-supplied numeric starts
  - `scale.est = "MAD"`, `"Huber"`, and `"proposal 2"`
  - `k2`
  - `acc`
  - `test.vec = "coef"`, `"resid"`, or `"w"`
- defaults `tuning = 4.685` when `method = "mm"`
- handles formula-layer extras in the R frontend:
  - `subset`
  - `na.action`
  - `contrasts`
  - `offset`
- handles custom function-valued psi in the package wrapper without delegating estimation to `MASS::rlm()`
- robust covariance is handled separately through `vcov()` and `confint()`

Example:

```r
set.seed(42)
n <- 200
x1 <- rnorm(n)
x2 <- rnorm(n)
y <- 1 + 2 * x1 - 1.5 * x2 + rnorm(n)

fit <- fit_rlm(y ~ x1 + x2, data = data.frame(y, x1, x2), method = "mm")
summary(fit)
confint(fit, hc_type = "HC3")

w <- runif(n, 0.8, 1.2)
fit_hampel <- fit_rlm(
  y ~ x1 + x2,
  data = data.frame(y, x1, x2, w = w),
  method = "m",
  psi = "hampel",
  weights = w,
  wt.method = "case",
  init = "lts",
  scale.est = "Huber",
  test.vec = "coef"
)
summary(fit_hampel)
```

### `fit_s_estimator()`

Fits a standalone S-estimator for robust linear regression.

Use this when you want:

- a high-breakdown robust regression fit directly
- access to S-estimation without the MM refinement stage
- a robust starting point conceptually separated from MM fitting

Example:

```r
X <- cbind(x1 = x1, x2 = x2)
s_fit <- fit_s_estimator(X, y)
summary(s_fit)
```

### `vcov_robust()` and `confint_robust()`

These provide HC-style sandwich inference for classical `lm` objects.

Supported covariance types:

- `HC0`
- `HC1`
- `HC2`
- `HC3`
- `HC4`
- `HC5`

Example:

```r
ols_fit <- lm(y ~ x1 + x2)
vcov_robust(ols_fit, hc_type = "HC3")
confint_robust(ols_fit, hc_type = "HC3")
```

This split is deliberate:

- `fit_rlm()` chooses the estimator
- `vcov()` / `confint()` choose the inference method

## Robust DML

### `fit_mm_dml()`

Fits a robust partially linear double machine learning model with an MM second stage.

Current wrapper supports:

- matrix interface: `x`, `d`, `y`
- data-frame interface: `outcome`, `treatment`, `controls`, `data`
- nuisance learners:
  - `"lasso"`
  - `"elastic_net"`
  - `"random_forest"`
  - `"hist_gradient_boosting"`
- `se_type = "analytic"` or `se_type = "bootstrap"`
- bootstrap repetition control and multicore bootstrap execution
- `summary()` output with coefficient table, learner, fold count, backend, and convergence info

Example:

```r
set.seed(42)
n <- 250
p <- 6
x <- matrix(rnorm(n * p), nrow = n, ncol = p)
d <- 0.9 * x[, 1] - 0.4 * x[, 2] + rnorm(n, sd = 0.35)
y <- 1.35 * d + 0.5 * x[, 1] - 0.2 * x[, 3] + rnorm(n, sd = 0.30)
toy_df <- data.frame(y = y, d = d, x)
names(toy_df) <- c("y", "d", paste0("x", 1:6))

fit <- fit_mm_dml(
  outcome = "y",
  treatment = "d",
  controls = paste0("x", 1:6),
  data = toy_df,
  learner = "lasso",
  se_type = "bootstrap",
  bootstrap_replications = 50,
  n_cores = 2L,
  folds = 3,
  seed = 101
)

summary(fit)
```

## Robust Adstock

### `build_adstock()`

Builds a carryover stock from a raw signal, with optional robust precleaning before the recursive adstock step.

Supported increment rules:

- `"plain"`
- `"huber"`
- `"tanh"`
- `"softsign"`
- `"adaptive_clip"`

Supported preclean methods:

- `"s"`
- `"mm"`

The function supports:

- direct vector/matrix input
- data-frame input
- summary output showing signal, cleaned signal, stock summaries, and head comparisons

Example:

```r
signal <- c(10, 12, 60, 9, 8, 40)
unit_ids <- c(1, 1, 1, 2, 2, 2)
x_preclean <- cbind("(Intercept)" = 1, trend = c(0, 1, 2, 0, 1, 2))

ad_fit <- build_adstock(
  signal = signal,
  unit_ids = unit_ids,
  x_preclean = x_preclean,
  preclean = TRUE,
  preclean_method = "mm",
  increment_method = "huber",
  rho = 0.8
)

summary(ad_fit)
```

## Robust MMM

### `fit_mmm()`

Fits a multi-channel marketing mix model with the following pipeline for each channel:

1. optional robust precleaning
2. carryover/adstock construction
3. optional Hill saturation
4. final regression using OLS, Huber, S, or MM

Important inputs:

- `channel_cols`: media channels
- `preclean_cols`: covariates used in the first-step cleaning stage
- `control_cols`: final regression controls
- `rho`: per-channel carryover parameters
- `increment_method`: per-channel adstock rule
- `hill_lambda`: per-channel Hill saturation parameter
- `fit_method`: one of `"ols"`, `"huber"`, `"s"`, `"mm"`

The fit object includes:

- full coefficient vector
- channel-only coefficients
- cleaned channel signals
- carryover stocks
- Hill-transformed channel regressors
- final design matrix
- structured `summary()` output with channel diagnostics and coefficient table

Example:

```r
mmm_fit <- fit_mmm(
  data = mmm_dat,
  outcome = "y",
  unit = "unit",
  channel_cols = c("channel_1", "channel_2", "channel_3"),
  preclean_cols = c("x1", "x2"),
  control_cols = c("x1", "x2"),
  fit_method = "huber",
  rho = c(0.8, 0.75, 0.7),
  increment_method = c("adaptive_clip", "huber", "tanh"),
  hill_lambda = c(0.12, 0.10, 0.08),
  preclean = c(TRUE, TRUE, TRUE),
  preclean_method = c("mm", "mm", "s")
)

summary(mmm_fit)
```

## C++ Core

The C++ library is the shared backend that the R package compiles against.

Main public headers:

- [include/robustcause/robustcause.hpp](/A:/RobustCause/include/robustcause/robustcause.hpp:1)
- [include/robustcause/capi.h](/A:/RobustCause/include/robustcause/capi.h:1)

The design is matrix-first and binding-friendly:

- column-major layout aligns naturally with Armadillo and R
- the C ABI is intended to be stable enough for future Python and R frontends
- higher-level ergonomics such as formulas and data frames belong mostly in the frontend layer

For `fit_rlm()`, this means:

- named robust-regression controls live in the native solver
- formula parsing and model-matrix construction live in the R frontend
- arbitrary user-defined psi callbacks are currently an R-wrapper capability rather than part of the language-agnostic C++ ABI

## Building The C++ Library

The core library depends on [Armadillo](https://arma.sourceforge.net/).

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build
```

Downstream CMake usage after install:

```cmake
find_package(RobustCause REQUIRED)
target_link_libraries(my_target PRIVATE RobustCause::robustcause)
```

## Notebook Showcase

The repository includes an end-to-end notebook:

- source: [bindings/r/robustcause/robustcause-smoke-test.Rmd](/A:/RobustCause/bindings/r/robustcause/robustcause-smoke-test.Rmd:1)
- rendered output: [bindings/r/robustcause/robustcause-smoke-test.html](/A:/RobustCause/bindings/r/robustcause/robustcause-smoke-test.html:1)

The notebook acts as both:

- a feature showcase
- a smoke test for package installation and the main estimators

It currently covers:

- OLS
- M-estimation
- MM-estimation
- S-estimation
- HC0-HC5 inference for `lm`
- MM-DML
- adstock
- MMM

## Roadmap

The next high-value steps are:

- add Python bindings under `bindings/python/`
- reduce duplication between the root C++ core and vendored R package sources
- add stronger regression tests against reference implementations
- expand the causal-estimation catalog on top of the current robust core
