# RobustCause

> Status: experimental / research software.

RobustCause is a C++-backed library for robust regression, robust causal inference, and robust marketing-mix modeling. It is built around a shared numerical core intended to support multiple frontends, with an R package available now and Python bindings planned next.

The project combines high-breakdown robust estimation, practical model-fitting workflows, and a reusable backend that can serve both causal inference and media-response applications. Today, the most complete public API is the R package under `bindings/r/robustcause`.

## Table of Contents

- [Introduction](#introduction)
- [Current Status](#current-status)
- [Installation](#installation)
  - [R Package Installation](#r-package-installation)
  - [Building the C++ Library](#building-the-c-library)
- [Quick Start](#quick-start)
- [What the API Does](#what-the-api-does)
- [Core R Functions](#core-r-functions)
  - [`fit_rlm()`](#fit_rlm)
  - [`fit_s_estimator()`](#fit_s_estimator)
  - [`vcov_robust()` and `confint_robust()`](#vcov_robust-and-confint_robust)
- [Robust DML](#robust-dml)
  - [`fit_mm_dml()`](#fit_mm_dml)
- [Robust Adstock](#robust-adstock)
  - [`build_adstock()`](#build_adstock)
  - [`fit_adstock_model()`](#fit_adstock_model)
- [Robust MMM](#robust-mmm)
  - [`fit_mmm()`](#fit_mmm)
- [C++ Core](#c-core)
- [Repository Layout](#repository-layout)
- [Notebook Showcase](#notebook-showcase)
- [Roadmap](#roadmap)

## Introduction

RobustCause is designed for users who want a shared robust estimation stack rather than a collection of disconnected functions. The library centers on a common C++ backend that supports robust linear regression, robust double machine learning, robust adstock construction, and robust marketing-mix modeling.

At the user level, the project currently exposes its most complete interface through the R package. That package provides a consistent API for fitting robust linear models, standalone S-estimators, robust partially linear double machine learning models, robust adstock transformations, end-to-end adstock response models, and multi-channel marketing-mix models.

Under the hood, these pieces are built on a shared numerical core intended to serve as the long-run engine for both R and Python frontends. The goal is to provide a coherent robust modeling framework that can be reused across regression, causal inference, and media-response workflows.

## Current Status

RobustCause currently includes:

- an installable C++ core library with CMake packaging
- shared C++ modules for robust regression, S-estimation, adstock, and MMM
- an R package wrapper with compiled native code
- native C++ MM-DML nuisance learners with tunable lasso and random-forest backends
- robust inference helpers for both RobustCause fits and base R `lm` objects
- a smoke-test and showcase notebook in `bindings/r/robustcause/robustcause-smoke-test.Rmd`

Not yet implemented:

- Python bindings
- a polished user-facing Python package
- a unified high-level frontend across R and Python

If you want to use RobustCause today, the R package is the primary public API.

## Installation

### R Package Installation

Install from GitHub:

```r
install.packages(c("remotes", "Rcpp", "RcppArmadillo", "RcppEigen"))
remotes::install_github("DanBaissa/RobustCause", subdir = "bindings/r/robustcause")
```

Then load it:

```r
library(robustcause)
```

### Building the C++ Library

The core library depends on Armadillo.

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build
```

For downstream CMake usage after installation:

```cmake
find_package(RobustCause REQUIRED)
target_link_libraries(my_target PRIVATE RobustCause::robustcause)
```

## Quick Start

Install the R package and fit a robust linear model:

```r
set.seed(42)
n <- 200
x1 <- rnorm(n)
x2 <- rnorm(n)
y <- 1 + 2 * x1 - 1.5 * x2 + rnorm(n)

fit <- fit_rlm(y ~ x1 + x2, data = data.frame(y, x1, x2), method = "mm")
summary(fit)
confint(fit, hc_type = "HC3")
```

From there, you can move to:

- `fit_mm_dml()` for robust partially linear double machine learning
- `build_adstock()` for robust carryover stock construction
- `fit_adstock_model()` for end-to-end adstock response modeling
- `fit_mmm()` for multi-channel marketing-mix models

## What the API Does

Today, the most complete user-facing API lives in the R package under `bindings/r/robustcause`.

That package exposes a consistent family of fit objects with `print()` and `summary()` methods for:

- `fit_rlm()` for robust linear models
- `fit_s_estimator()` for standalone S-estimation
- `fit_mm_dml()` for robust double machine learning
- `build_adstock()` for robust carryover stock construction
- `fit_adstock_model()` for end-to-end single-signal adstock response modeling
- `fit_mmm()` for multi-channel marketing-mix models

The underlying C++ core lives in:

- `include/robust`
- `include/robustcause`
- `src`

That core is intended to be the shared engine for both the current R package and future Python bindings.

## Core R Functions

### `fit_rlm()`

`fit_rlm()` fits a robust linear model using either:

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
- separates estimation from inference, with covariance and intervals handled through `vcov()` and `confint()`

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

`fit_s_estimator()` fits a standalone S-estimator for robust linear regression.

Use it when you want:

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

These functions provide HC-style sandwich inference for classical `lm` objects.

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
- `vcov()` and `confint()` choose the inference method

## Robust DML

### `fit_mm_dml()`

`fit_mm_dml()` fits a robust partially linear double machine learning model with an MM second stage.

Current wrapper supports:

- matrix interface: `x`, `d`, `y`
- data-frame interface: `outcome`, `treatment`, `controls`, `data`
- nuisance learners:
  - `"lasso"`
  - `"elastic_net"`
  - `"random_forest"`
  - `"hist_gradient_boosting"`
- `learner_params` for learner-specific tuning
- `se_type = "analytic"` or `se_type = "bootstrap"`
- bootstrap repetition control and multicore bootstrap execution
- `summary()` output with coefficient table, learner, fold count, backend, convergence info, and native learner diagnostics when available

For the native compiled single-treatment random-fold path:

- `"lasso"` now exposes:
  - `cv_folds`
  - `n_lambda`
  - `lambda_min_ratio`
  - optional `lambda_grid`
  - `use_lambda_1se`
  - `standardize`
  - `max_iter`
  - `tolerance`
- `"random_forest"` now exposes:
  - `n_estimators`
  - `max_depth`
  - `min_samples_split`
  - `min_samples_leaf`
  - `max_features`
  - `max_features_fraction`
  - `sample_fraction`
  - `bootstrap`
  - `replacement`
  - `split_candidates`
  - `compute_oob`
  - `compute_importance`

Native fit objects now carry learner diagnostics in `fit$learner_details`. For lasso that includes selected lambda, `lambda.min`, `lambda.1se`, CV error, and nonzero count. For random forest that includes OOB error, tree count, feature importance, and feature use counts.

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
  learner_params = list(
    cv_folds = 4L,
    n_lambda = 80L,
    use_lambda_1se = TRUE
  ),
  se_type = "bootstrap",
  bootstrap_replications = 50,
  n_cores = 2L,
  folds = 3,
  seed = 101
)

summary(fit)
```

Native random-forest tuning works the same way:

```r
fit_rf <- fit_mm_dml(
  outcome = "y",
  treatment = "d",
  controls = paste0("x", 1:6),
  data = toy_df,
  learner = "random_forest",
  learner_params = list(
    n_estimators = 200L,
    max_depth = 7L,
    min_samples_split = 8L,
    min_samples_leaf = 4L,
    sample_fraction = 0.7
  ),
  folds = 3,
  seed = 101
)

summary(fit_rf)
fit_rf$learner_details
```

Use `make_custom_learner()` when you want to plug in a fully custom nuisance model such as your own neural net or GBM pipeline. Use named learners plus `learner_params` when you want to stay on the compiled native C++ path.

## Robust Adstock

### `build_adstock()`

`build_adstock()` builds a carryover stock from a raw signal, with optional robust precleaning before the recursive adstock step.

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

- direct vector or matrix input
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

Use `build_adstock()` when you want the transformed stock itself. If you want the full adstock-response model with parameter search and a downstream regression, use `fit_adstock_model()`.

### `fit_adstock_model()`

`fit_adstock_model()` fits a complete single-signal adstock response model end to end.

What it does:

- searches over `rho` and Hill-saturation candidates
- builds the stock and transformed stock
- fits the downstream response model
- supports full robustness levels:
  - `"classical"`
  - `"huber"`
  - `"s"`
  - `"mm"`
- returns the selected model along with the full search table

Example:

```r
fit_ad <- fit_adstock_model(
  data = dat,
  outcome = "y",
  signal_col = "signal",
  unit_col = "unit",
  control_cols = "x1",
  robustness = "classical"
)

summary(fit_ad)
```

## Robust MMM

### `fit_mmm()`

`fit_mmm()` fits a multi-channel marketing-mix model with the following per-channel pipeline:

1. optional robust precleaning
2. carryover or adstock construction
3. optional Hill saturation
4. final regression using OLS, Huber, S, or MM

It supports two workflows:

- fixed-specification MMM, where you provide `rho` and `hill_lambda` directly
- searched MMM, where you provide `rho_grid` and or `hill_grid` and rank candidate models by `BIC`, `AIC`, or `MSE`

It also supports high-level robustness profiles:

- `robustness = "classical"`
- `robustness = "huber"`
- `robustness = "s"`
- `robustness = "mm"`

Important inputs:

- `channel_cols`: media channels
- `preclean_cols`: covariates used in the first-step cleaning stage
- `control_cols`: final regression controls
- `rho`: per-channel fixed carryover parameters
- `rho_grid`: optional shared or per-channel carryover search grid
- `increment_method`: per-channel adstock rule
- `hill_lambda`: per-channel fixed Hill saturation parameter
- `hill_grid`: optional shared or per-channel Hill search grid
- `fit_method`: one of `"ols"`, `"huber"`, `"s"`, `"mm"`

The fit object includes:

- full coefficient vector
- channel-only coefficients
- cleaned channel signals
- carryover stocks
- Hill-transformed channel regressors
- final design matrix
- optional `grid_results` table when search is used
- structured `summary()` output with channel diagnostics, selected channel configuration, and top search candidates

Example:

```r
mmm_fit <- fit_mmm(
  data = mmm_dat,
  outcome = "y",
  unit = "unit",
  channel_cols = c("channel_1", "channel_2", "channel_3"),
  preclean_cols = c("x1", "x2"),
  control_cols = c("x1", "x2"),
  robustness = "mm",
  rho_grid = list(c(0.6, 0.8), c(0.5, 0.75), c(0.4, 0.7)),
  hill_grid = list(c(0.05, 0.12), c(0.05, 0.10), c(0.04, 0.08)),
  choose_by = "bic"
)

summary(mmm_fit)
```

## C++ Core

The C++ library is the shared backend that the R package compiles against.

Main public headers:

- `include/robustcause/robustcause.hpp`
- `include/robustcause/capi.h`

The design is matrix-first and binding-friendly:

- column-major layout aligns naturally with Armadillo and R
- the C ABI is intended to be stable enough for future Python and R frontends
- higher-level ergonomics such as formulas and data frames belong mostly in the frontend layer

For `fit_rlm()`, that means:

- named robust-regression controls live in the native solver
- formula parsing and model-matrix construction live in the R frontend
- arbitrary user-defined psi callbacks are currently an R-wrapper capability rather than part of the language-agnostic C++ ABI

## Repository Layout

- `include/robust`: core C++ headers
- `include/robustcause`: umbrella headers and C ABI
- `src`: core C++ implementations
- `examples`: C++ example programs
- `tests`: smoke tests for the C++ core
- `bindings/r/robustcause`: R package
- `docs`: architecture notes

## Notebook Showcase

The repository includes an end-to-end notebook:

- source: `bindings/r/robustcause/robustcause-smoke-test.Rmd`
- rendered output: `bindings/r/robustcause/robustcause-smoke-test.html`

The notebook acts as both:

- a feature showcase
- a smoke test for package installation and the main estimators

It currently covers:

- OLS
- M-estimation
- MM-estimation
- S-estimation
- HC0 through HC5 inference for `lm`
- MM-DML
- adstock
- end-to-end adstock modeling
- MMM

## Roadmap

The next high-value steps are:

- add Python bindings under `bindings/python/`
- reduce duplication between the root C++ core and vendored R package sources
- add stronger regression tests against reference implementations
- expand the causal-estimation catalog on top of the current robust core
