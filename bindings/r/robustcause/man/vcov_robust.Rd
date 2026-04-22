\name{vcov_robust}
\alias{vcov_robust}
\alias{confint_robust}
\title{Robust Covariance For Fitted Models}
\usage{
vcov_robust(object, hc_type = "HC3")

confint_robust(object, hc_type = "HC3", level = 0.95)
}
\arguments{
  \item{object}{A \code{robustcause_rlm} object or a base \code{lm} object.}
  \item{hc_type}{HC covariance estimator name.}
  \item{level}{Confidence level for \code{confint_robust}.}
}
\description{
Computes heteroskedasticity-consistent covariance estimates and confidence
intervals for fitted RobustCause models or ordinary least squares \code{lm}
objects.
}
\details{
\strong{When to use these helpers:}

Use \code{vcov_robust()} and \code{confint_robust()} when you want robust
inference but do not want that choice to be baked into model fitting. This is
especially useful when:

\itemize{
\item you fit \code{fit_rlm()} and want to compare different HC covariance estimators
\item you fit a base \code{lm} object and want sandwich-style standard errors
}

\strong{Practical default:}

\code{HC3} is a good default in many applied settings because it behaves well in
moderate samples and under heteroskedasticity.
}
\examples{
set.seed(1)
n <- 100
x <- rnorm(n)
y <- 1 + 2 * x + rnorm(n, sd = abs(x) + 0.2)

ols_fit <- lm(y ~ x)
vcov_robust(ols_fit, hc_type = "HC3")
confint_robust(ols_fit, hc_type = "HC3")
}
