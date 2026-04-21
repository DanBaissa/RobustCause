# robustcause

`robustcause` is an R package wrapper around the shared RobustCause C++ core.

Current exported entry points:
- `fit_rlm()`
- `fit_s_estimator()`
- `vcov()`, `confint()`, and `plot_weights()` methods for fitted objects

The package compiles vendored copies of the current core sources so it can be
built as a self-contained R package inside this repository.
