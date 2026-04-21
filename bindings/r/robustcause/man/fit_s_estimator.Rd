\name{fit_s_estimator}
\alias{fit_s_estimator}
\title{Fit a Robust S-Estimator}
\usage{
fit_s_estimator(
  x,
  y = NULL,
  data = NULL,
  n_starts = 500L,
  n_best_starts = 25L,
  max_refine = 100L,
  max_scale_iter = 100L,
  tol = 1e-08,
  scale_tol = 1e-10,
  c = 1.54764,
  b = 0.5,
  add_intercept = TRUE,
  ridge = 1e-10,
  min_weight = 1e-12,
  seed = 123456789,
  use_fast_s = TRUE,
  include_ols_start = TRUE,
  fast_s_screen_subsets = 500L,
  fast_s_screen_iters = 2L,
  fast_s_keep = 25L
)
}
\arguments{
  \item{x}{A numeric design matrix or a formula.}
  \item{y}{Numeric response vector when \code{x} is not a formula.}
  \item{data}{A data frame used with formula input.}
  \item{n_starts}{Number of random starts.}
  \item{n_best_starts}{Number of best starts to keep before full refinement.}
  \item{max_refine}{Maximum refinement iterations.}
  \item{max_scale_iter}{Maximum scale-solver iterations.}
  \item{tol}{Coefficient convergence tolerance.}
  \item{scale_tol}{Scale-solver tolerance.}
  \item{c}{Tukey bisquare constant.}
  \item{b}{Target rho expectation.}
  \item{add_intercept}{Whether to add an intercept for matrix input.}
  \item{ridge}{Small ridge term used in linear solves.}
  \item{min_weight}{Lower bound for weights in the solver.}
  \item{seed}{Random seed passed to the C++ core.}
  \item{use_fast_s}{Whether to use Fast-S-style screening.}
  \item{include_ols_start}{Whether to include an OLS starting point.}
  \item{fast_s_screen_subsets}{Number of screening subsets.}
  \item{fast_s_screen_iters}{Screening refinement iterations.}
  \item{fast_s_keep}{Number of screening candidates to keep.}
}
\value{
An object of class \code{robustcause_s}.
}
\description{
Fits the RobustCause Tukey bisquare S-estimator through the compiled C++ core.
}
