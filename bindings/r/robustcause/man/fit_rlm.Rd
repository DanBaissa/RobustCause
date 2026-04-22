\name{fit_rlm}
\alias{fit_rlm}
\title{Fit a Robust Linear Model}
\usage{
fit_rlm(
  x,
  y = NULL,
  data = NULL,
  method = c("m", "mm"),
  psi = NULL,
  tuning = NULL,
  maxit = 100L,
  tol = 1e-08,
  add_intercept = TRUE,
  ridge = 1e-10,
  min_weight = 1e-12,
  mm_s_control = NULL
)
}
\arguments{
  \item{x}{A numeric design matrix or a formula.}
  \item{y}{Numeric response vector when \code{x} is not a formula.}
  \item{data}{A data frame used with formula input.}
  \item{method}{Either \code{"m"} or \code{"mm"}.}
  \item{psi}{Psi function name. If \code{NULL}, the default is Huber for M and Tukey bisquare for MM.}
  \item{tuning}{Psi tuning constant. If \code{NULL}, the default is \code{1.345} for M and \code{4.685} for MM.}
  \item{maxit}{Maximum IRLS iterations.}
  \item{tol}{Convergence tolerance.}
  \item{add_intercept}{Whether to add an intercept for matrix input.}
  \item{ridge}{Small ridge term used in linear solves.}
  \item{min_weight}{Lower bound for weights in the solver.}
  \item{mm_s_control}{Optional named list of S-estimator controls used to initialize MM fits.}
}
\value{
An object of class \code{robustcause_rlm}.
}
\description{
Fits the RobustCause robust linear model using either a plain M-estimator or an MM-estimator seeded by the S-estimator.
}
