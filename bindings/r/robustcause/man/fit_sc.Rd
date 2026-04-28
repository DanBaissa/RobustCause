\name{fit_sc}
\alias{fit_sc}
\alias{summary.robustcause_sc}
\alias{print.robustcause_sc}
\alias{print.summary.robustcause_sc}
\alias{plot.robustcause_sc}
\alias{coef.robustcause_sc}
\alias{fitted.robustcause_sc}
\alias{residuals.robustcause_sc}
\alias{predict.robustcause_sc}
\alias{as.data.frame.robustcause_sc}
\title{Fit Standard or Robust Synthetic Control}
\usage{
fit_sc(
  outcomes = NULL,
  treated_unit = NULL,
  treatment_start = NULL,
  donors = NULL,
  unit = NULL,
  time = NULL,
  outcome = NULL,
  data = NULL,
  predictors = NULL,
  predictor_agg = c("mean", "last_pre"),
  predictor_weights = NULL,
  method = c("standard", "mm"),
  run_placebos = TRUE,
  na.action = c("omit", "fail"),
  maxit = 500L,
  tol = 1e-08,
  ridge = 1e-10,
  startup_maxit = 75L,
  subproblem_maxit = 500L,
  subproblem_tol = 1e-08,
  l1_smoothing = 1e-06,
  tukey_c = 4.685,
  min_scale = 1e-08,
  min_time_weight = 1e-08
)

\method{coef}{robustcause_sc}(object, ...)

\method{fitted}{robustcause_sc}(object, ...)

\method{residuals}{robustcause_sc}(object, ...)

\method{predict}{robustcause_sc}(
  object,
  newdata = NULL,
  type = c("synthetic", "gap"),
  period = c("all", "pre", "post"),
  ...
)

\method{as.data.frame}{robustcause_sc}(x, ...)
}
\arguments{
  \item{outcomes}{Numeric matrix with rows as time periods and columns as units.}
  \item{treated_unit}{Treated unit identifier. For matrix input this can be a column name or index. For data-frame input it must match values in \code{unit}.}
  \item{treatment_start}{First treated period. For matrix input use a row index. For data-frame input use a value from the time column.}
  \item{donors}{Optional donor-unit identifiers or matrix column indices. Defaults to all non-treated units.}
  \item{unit}{Unit column name for panel data input.}
  \item{time}{Time column name for panel data input.}
  \item{outcome}{Outcome column name for panel data input.}
  \item{data}{Optional panel data frame.}
  \item{predictors}{Optional predictor specification. For matrix input, a numeric matrix with one column per unit. For data-frame input, one or more column names to aggregate over the pre-treatment period.}
  \item{predictor_agg}{How to aggregate panel predictors over pre-treatment periods.}
  \item{predictor_weights}{Optional positive weights for predictor rows.}
  \item{method}{Either standard synthetic control or robust MM synthetic control.}
  \item{run_placebos}{Whether to run placebo fits over donor units.}
  \item{na.action}{Missing-data handling. \code{"omit"} drops incomplete rows or periods before fitting.}
  \item{maxit}{Maximum iterations for the main solver.}
  \item{tol}{Convergence tolerance.}
  \item{ridge}{Small ridge stabilization used in the simplex subproblem.}
  \item{startup_maxit}{Maximum startup iterations for MM synthetic control.}
  \item{subproblem_maxit}{Maximum projected-gradient iterations for the simplex subproblem.}
  \item{subproblem_tol}{Tolerance for the simplex subproblem.}
  \item{l1_smoothing}{Smoothing parameter for the startup L1 approximation.}
  \item{tukey_c}{Tukey bisquare tuning constant for MM synthetic control.}
  \item{min_scale}{Lower bound for the robust residual scale.}
  \item{min_time_weight}{Lower bound for MM time weights.}
  \item{object}{A fitted \code{robustcause_sc} object.}
  \item{x}{A fitted \code{robustcause_sc} object.}
  \item{newdata}{Not used. Synthetic-control predictions are tied to the fitted panel.}
  \item{type}{For \code{predict()}, return either the synthetic path or the gap path.}
  \item{period}{For \code{predict()}, subset the returned path to all, pre-, or post-treatment periods.}
  \item{\dots}{Additional arguments passed to methods.}
}
\value{
For \code{fit_sc()}, an object of class \code{robustcause_sc}. Core elements
include \code{weights}, \code{synthetic_pre}, \code{synthetic_post},
\code{pre_residuals}, \code{post_gaps}, \code{path}, and placebo/inference
outputs when requested.

For \code{coef()}, the donor-weight vector.

For \code{fitted()}, the full synthetic path over pre- and post-treatment
periods.

For \code{residuals()}, the full gap path over pre- and post-treatment
periods.

For \code{predict()}, the requested synthetic or gap path.

For \code{as.data.frame()}, a data frame with columns \code{time},
\code{period}, \code{treated}, \code{synthetic}, and \code{gap}.
}
\description{
Fits standard or MM-style robust synthetic control through the compiled C++ backend.
}
\details{
\strong{Interfaces:}

\itemize{
\item Matrix interface: pass \code{outcomes}, \code{treated_unit}, and \code{treatment_start}.
\item Panel-data interface: pass \code{data}, \code{unit}, \code{time}, \code{outcome}, \code{treated_unit}, and \code{treatment_start}.
}

\strong{Robust synthetic control:}

With \code{method = "mm"}, the estimator uses a simplex-constrained robust startup, a MAD residual scale, and Tukey-bisquare refinement. This downweights contaminated pre-treatment periods but does not solve donor-pool mismatch or lack of overlap.

\strong{Inference:}

If \code{run_placebos = TRUE}, the wrapper runs donor-placebo fits and reports RMSPE-ratio and average-post-gap placebo p-values.
}
\examples{
outcomes <- cbind(
  treated = c(1.0, 2.1, 9.0, 4.0, 5.1, 6.0, 7.0),
  donor_a = c(1.1, 2.0, 3.0, 4.1, 5.0, 6.1, 7.2),
  donor_b = c(0.8, 1.9, 3.2, 3.9, 5.2, 5.9, 6.9),
  donor_c = c(1.3, 2.2, 2.8, 4.2, 4.9, 6.0, 7.1)
)

fit <- fit_sc(
  outcomes = outcomes,
  treated_unit = "treated",
  treatment_start = 6,
  method = "mm",
  run_placebos = TRUE
)

summary(fit)
coef(fit)
fitted(fit)
predict(fit, type = "gap", period = "post")
head(as.data.frame(fit))
plot(fit)
plot(fit, type = "weights")
plot(fit, type = "placebo_gaps")
}
