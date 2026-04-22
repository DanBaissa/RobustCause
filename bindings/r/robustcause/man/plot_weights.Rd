\name{plot_weights}
\alias{plot_weights}
\title{Plot Robust Weights}
\usage{
plot_weights(object, ...)
}
\arguments{
  \item{object}{A fitted \code{robustcause_rlm} or \code{robustcause_s} object.}
  \item{...}{Additional graphical parameters passed to \code{plot}.}
}
\description{
Plots per-observation robustness weights from a fitted model.
}
\details{
This is a quick diagnostic for how aggressively a robust fit downweighted
individual observations. Large clusters of weights near zero suggest the model
is treating many rows as highly incompatible with the fitted relationship.

Typical workflow:
\enumerate{
\item fit a robust model with \code{fit_rlm()} or \code{fit_s_estimator()}
\item inspect \code{summary(fit)}
\item use \code{plot_weights(fit)} to see which observations received little weight
}
}
\seealso{
\code{\link{vcov_robust}}, \code{\link{confint_robust}}
}
