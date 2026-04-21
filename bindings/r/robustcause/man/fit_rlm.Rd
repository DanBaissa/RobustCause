\name{fit_rlm}
\alias{fit_rlm}
\title{Fit a Robust Linear Model}
\usage{
fit_rlm(
  x,
  y = NULL,
  data = NULL,
  psi = c("huber", "tukey_bisquare"),
  tuning = 1.345,
  maxit = 100L,
  tol = 1e-08,
  add_intercept = TRUE,
  ridge = 1e-10,
  min_weight = 1e-12,
  hc_type = "HC3",
  zcrit = qnorm(0.975)
)
}
\arguments{
  \item{x}{A numeric design matrix or a formula.}
  \item{y}{Numeric response vector when \code{x} is not a formula.}
  \item{data}{A data frame used with formula input.}
  \item{psi}{Psi function name.}
  \item{tuning}{Psi tuning constant.}
  \item{maxit}{Maximum IRLS iterations.}
  \item{tol}{Convergence tolerance.}
  \item{add_intercept}{Whether to add an intercept for matrix input.}
  \item{ridge}{Small ridge term used in linear solves.}
  \item{min_weight}{Lower bound for weights in the solver.}
  \item{hc_type}{HC covariance estimator name.}
  \item{zcrit}{Critical value used for normal confidence intervals.}
}
\value{
An object of class \code{robustcause_rlm}.
}
\description{
Fits the RobustCause robust linear model and computes HC sandwich inference.
}
