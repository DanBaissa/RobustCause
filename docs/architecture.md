# Architecture

## Core principle

This repository is the shared numerical backend for robust causal-inference tooling. Language packages should be thin wrappers over a common C++ implementation so that estimation logic, defaults, and numerical edge-case handling remain aligned.

## Layers

### 1. C++ estimator layer

Implemented today in `include/robust/` and `src/`.

- `fit_rlm(...)`
- `vcov_hc(...)`
- `confint_normal(...)`
- `fit_s_estimator(...)`
- `prepare_sc_data(...)`
- `fit_sc(...)`
- `fit_mm_sc(...)`
- `fit_sc_placebos(...)`
- `fit_mm_sc_placebos(...)`
- `build_model_frame(...)`

This layer is the richest API and is the best place to evolve new algorithms first.

### 2. C ABI layer

Implemented in `include/robustcause/capi.h` and `src/capi.cpp`.

This layer exists to make language bindings easier to write and easier to keep stable. It currently focuses on matrix-based entry points and plain buffers:

- no STL types across the boundary
- no exceptions across the boundary
- caller-owned output buffers
- column-major matrix layout

That layout matches Armadillo and R naturally. Python wrappers can either require Fortran-order arrays or perform a conversion at the boundary.

### 3. Language frontends

Not implemented yet, but the expected structure is:

- `bindings/python/`: `pybind11` module exposing NumPy-first APIs
- `bindings/r/`: `Rcpp` package exposing R matrices, formulas, and data frames

Frontend packages should:
- translate user-facing data structures into the matrix-first core
- handle formula parsing in the host language when that is more idiomatic
- keep option objects aligned with the core defaults

## API direction

The matrix API should be treated as the stable contract.

Formula parsing is useful in C++, but for R and Python it is usually better as a frontend convenience layer. That keeps the core small and makes language-specific semantics explicit instead of baking them into the shared backend.
