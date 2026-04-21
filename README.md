# RobustCause

RobustCause is a C++ core library for robust causal-inference tooling, intended to sit underneath R and Python frontends.

The current focus is the numerical kernel:
- robust adstock construction with optional S/MM precleaning
- robust MMM fitting with per-channel cleaning, carryover, Hill saturation, and robust regression
- robust linear regression with M and MM estimation
- Tukey bisquare S-estimation for linear models
- robust MM-DML for partially linear treatment-effect estimation
- HC0-HC5 heteroskedasticity-consistent covariance estimators
- a lightweight formula/data interface on the C++ side
- a small C ABI that can be wrapped from `pybind11`, `Rcpp`, or direct FFI

## Current status

The repository now includes:
- a C++ core library with installable CMake packaging
- a robust adstock module in the shared C++ core
- a robust MMM module in the shared C++ core
- an R package at `bindings/r/robustcause`
- M-estimation and MM-estimation in `fit_rlm()`
- MM-DML in `fit_mm_dml()`
- a separate robust covariance layer so HC inference is not tied to fitting
- robust covariance helpers for both RobustCause fits and base R `lm` objects
- a GitHub-installing R Markdown notebook that acts as both a showcase and smoke test

In the R package, the intended split is:
- fitting: `fit_rlm(..., method = "m" | "mm")`, `fit_s_estimator(...)`, and `fit_mm_dml(...)`
- inference: `vcov(..., hc_type = ...)`, `confint(..., hc_type = ...)`, `vcov_robust()`, and `confint_robust()`

## Design goals

- Keep the numerical core in C++ for speed and shared behavior.
- Expose matrix-first APIs that are easy to bind from Python and R.
- Treat formula parsing and dataframe ergonomics as frontend concerns when possible.
- Ship a stable, installable CMake target for downstream packages.

## Repository layout

- `include/robust/`: primary C++ headers for the current core API
- `include/robustcause/`: umbrella headers and the C ABI boundary
- `src/`: implementation files
- `examples/`: example consumers of the C++ and C APIs
- `tests/`: smoke tests
- `docs/`: architecture and binding notes
- `bindings/r/robustcause/`: self-contained R package wrapper around the C++ core
- `bindings/`: language bindings and related package scaffolding

## Build

The library depends on [Armadillo](https://arma.sourceforge.net/).

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build
```

## Consuming from CMake

After install, downstream CMake projects can use:

```cmake
find_package(RobustCause REQUIRED)
target_link_libraries(my_target PRIVATE RobustCause::robustcause)
```

## Binding strategy

Two layers are available:
- C++ API via `#include <robustcause/robustcause.hpp>`
- C ABI via `#include <robustcause/capi.h>`

The C ABI is deliberately matrix-first and column-major so it aligns well with Armadillo and R, while Python wrappers can opt into Fortran-order arrays or transpose on entry.

## R package

The R package currently provides:
- `fit_rlm()` with `method = "m"` and `method = "mm"`
- `fit_s_estimator()`
- `fit_mm_dml()` for robust orthogonalized treatment-effect estimation
- `build_adstock()` for robust carryover construction
- `fit_mmm()` for multi-channel MMM fitting
- `vcov()` and `confint()` methods for robust fits with selectable HC type
- `vcov_robust()` and `confint_robust()` for base `lm` objects
- a smoke-test/showcase notebook at `bindings/r/robustcause/robustcause-smoke-test.Rmd`

The notebook installs the package from GitHub before running, so it doubles as an end-to-end installation check.

## Adstock API

The shared C++ core now also exposes:
- `build_adstock(...)`
- `build_robust_adstock(...)`
- `preclean_signal_s(...)`
- `preclean_signal_mm(...)`
- `fit_mmm(...)`

The adstock module supports plain, Huber, tanh, softsign, and adaptive-clip increment rules, plus optional first-step S or MM signal cleaning before stock accumulation.

The MMM module builds on top of those adstock primitives and adds:
- per-channel precleaning
- per-channel Hill saturation
- final OLS, Huber, S, or MM regression

## Next recommended steps

- add `pybind11` bindings under `bindings/python/`
- grow the estimator catalog with causal-model-specific routines
- reduce source duplication between the root C++ core and the vendored R package sources
- formalize tests against reference implementations in R and Python
