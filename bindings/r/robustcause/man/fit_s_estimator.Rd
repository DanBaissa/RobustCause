\name{fit_s_estimator}
\alias{fit_s_estimator}
\alias{print.robustcause_s}
\alias{summary.robustcause_s}
\alias{print.summary.robustcause_s}
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
\details{
\strong{When to use this model:}

The S-estimator is the high-breakdown robust regression option in the package.
Use it when you expect substantial contamination, leverage points, or a small
fraction of observations that may be badly corrupted. It is also the conceptual
starting point for MM regression.

\strong{How it differs from \code{fit_rlm(method = "mm")}:}

\itemize{
\item \code{fit_s_estimator()} returns the pure S-estimator.
\item \code{fit_rlm(method = "mm")} uses an S-estimator for initialization and then refines to an MM solution.
\item In practice, use S when you want the higher-breakdown fit itself. Use MM when you want a more efficient final estimator with S-style protection at the start.
}

\strong{Tuning and starts:}

\itemize{
\item Increase \code{n_starts} when the problem is difficult or when you suspect many local minima.
\item Leave \code{use_fast_s = TRUE} on unless you are debugging or doing very small exact comparisons.
\item The \code{fast_s_*} controls govern the screening stage. Larger values are slower but can improve stability on harder problems.
}

\strong{What to inspect after fitting:}

\itemize{
\item \code{summary(fit)} for the fitted coefficients and scale
\item \code{fit\$weights} and \code{plot_weights(fit)} to see how strongly observations are being downweighted
}
}
\examples{
set.seed(1)
n <- 150
x1 <- rnorm(n)
x2 <- rnorm(n)
y <- 1 + 1.8 * x1 - 0.9 * x2 + rnorm(n, sd = 0.25)
y[c(5, 25, 90)] <- y[c(5, 25, 90)] + c(7, -8, 6)

fit <- fit_s_estimator(y ~ x1 + x2, data = data.frame(y, x1, x2))
summary(fit)
}
