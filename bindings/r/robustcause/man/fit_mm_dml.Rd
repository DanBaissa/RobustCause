\name{fit_mm_dml}
\alias{fit_mm_dml}
\alias{print.robustcause_mm_dml}
\alias{summary.robustcause_mm_dml}
\title{Fit MM Double Machine Learning}
\usage{
fit_mm_dml(
  x = NULL,
  d = NULL,
  y = NULL,
  outcome = NULL,
  treatment = NULL,
  controls = NULL,
  data = NULL,
  learner = c("lasso", "elastic_net", "random_forest", "hist_gradient_boosting"),
  se_type = c("analytic", "bootstrap"),
  bootstrap_replications = 200L,
  bootstrap_seed = NULL,
  progress = interactive(),
  n_cores = 1L,
  folds = 2L,
  fold_type = c("random"),
  seed = 123L,
  mm_c = 4.685,
  max_iter = 100L,
  tolerance = 1e-07,
  ci_level = 0.95
)
}
\arguments{
  \item{x}{Numeric control matrix.}
  \item{d}{Numeric treatment vector.}
  \item{y}{Numeric outcome vector.}
  \item{data}{Optional data frame interface.}
  \item{outcome}{Outcome column name in \code{data}.}
  \item{treatment}{Treatment column name in \code{data}.}
  \item{controls}{Control column names in \code{data}.}
  \item{learner}{Nuisance learner type.}
  \item{se_type}{Either \code{"analytic"} or \code{"bootstrap"}.}
  \item{bootstrap_replications}{Number of bootstrap resamples when \code{se_type = "bootstrap"}.}
  \item{bootstrap_seed}{Optional bootstrap seed override.}
  \item{progress}{Whether to show bootstrap progress.}
  \item{n_cores}{Requested bootstrap worker count.}
  \item{folds}{Cross-fitting fold count.}
  \item{fold_type}{Currently only \code{"random"} is supported.}
  \item{seed}{Random seed.}
  \item{mm_c}{Tukey bisquare constant for the MM second stage.}
  \item{max_iter}{Maximum MM stage iterations.}
  \item{tolerance}{Convergence tolerance for the MM stage.}
  \item{ci_level}{Normal-approximation confidence interval level.}
}
\value{
An object of class \code{robustcause_mm_dml}.
}
\description{
Fits a matrix-first MM double machine learning estimator using the native C++
backend adapted from the Robust_ML project, with an R-facing wrapper that adds
bootstrap standard errors and a data-frame interface.
}
