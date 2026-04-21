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
