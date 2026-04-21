# RobustCause

RobustCause is a C++ core library for robust causal-inference tooling, intended to sit underneath R and Python frontends.

The current focus is the numerical kernel:
- robust linear regression with HC sandwich inference
- Tukey bisquare S-estimation for linear models
- a lightweight formula/data interface on the C++ side
- a small C ABI that can be wrapped from `pybind11`, `Rcpp`, or direct FFI

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
- `bindings/`: placeholder directory for future R and Python packages

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

## Next recommended steps

- add `pybind11` bindings under `bindings/python/`
- add `Rcpp` bindings under `bindings/r/`
- grow the estimator catalog with causal-model-specific routines
- formalize tests against reference implementations in R and Python
