\name{fit_mm_dml}
\alias{fit_mm_dml}
\alias{print.robustcause_mm_dml}
\alias{summary.robustcause_mm_dml}
\alias{print.summary.robustcause_mm_dml}
\title{Fit MM Double Machine Learning}
\usage{
fit_mm_dml(
  x = NULL,
  d = NULL,
  y = NULL,
  outcome = NULL,
  treatment = NULL,
  controls = NULL,
  second_stage_controls = NULL,
  data = NULL,
  learner = c("lasso", "elastic_net", "random_forest", "hist_gradient_boosting"),
  learner_params = NULL,
  se_type = c("analytic", "bootstrap"),
  bootstrap_replications = 200L,
  bootstrap_seed = NULL,
  progress = interactive(),
  n_cores = 1L,
  folds = 2L,
  fold_type = c("random", "spatial_block"),
  block_id = NULL,
  spatial_k = NULL,
  spatial_by = NULL,
  seed = 123L,
  mm_c = 4.685,
  max_iter = 100L,
  tolerance = 1e-07,
  ci_level = 0.95
)
}
\arguments{
  \item{x}{Numeric control matrix for the matrix interface.}
  \item{d}{Numeric treatment vector or matrix for the matrix interface.}
  \item{y}{Numeric outcome vector for the matrix interface.}
  \item{outcome}{Outcome column name in \code{data}.}
  \item{treatment}{One or more treatment column names in \code{data}.}
  \item{controls}{Control column names in \code{data}.}
  \item{second_stage_controls}{Optional columns to include directly in the final MM second stage after nuisance residualization.}
  \item{data}{Optional data frame interface.}
  \item{learner}{A named nuisance learner or a custom learner object created with \code{make_custom_learner()}.}
  \item{learner_params}{Optional named list of learner-specific tuning controls. Native C++ \code{"lasso"} and \code{"random_forest"} now expose substantial tuning knobs through this argument.}
  \item{se_type}{Either \code{"analytic"} or \code{"bootstrap"}.}
  \item{bootstrap_replications}{Number of bootstrap resamples when \code{se_type = "bootstrap"}.}
  \item{bootstrap_seed}{Optional bootstrap seed override.}
  \item{progress}{Whether to show bootstrap progress.}
  \item{n_cores}{Requested bootstrap worker count.}
  \item{folds}{Cross-fitting fold count.}
  \item{fold_type}{Fold assignment strategy. Use \code{"random"} or \code{"spatial_block"}.}
  \item{block_id}{When \code{fold_type = "spatial_block"}, a column name in \code{data} identifying the block for each row.}
  \item{spatial_k}{Reserved spatial learner tuning parameter for API compatibility with the \code{Robust_ML} interface.}
  \item{spatial_by}{Reserved spatial learner modifier for API compatibility with the \code{Robust_ML} interface.}
  \item{seed}{Random seed.}
  \item{mm_c}{Tukey bisquare constant for the MM second stage.}
  \item{max_iter}{Maximum MM stage iterations.}
  \item{tolerance}{Convergence tolerance for the MM stage.}
  \item{ci_level}{Confidence interval level.}
}
\value{
An object of class \code{robustcause_mm_dml}.
}
\description{
Fits a robust MM double machine learning estimator. For the simple single-treatment random-fold case with built-in compiled learners, the package uses the native C++ backend. For richer use cases such as multiple treatments, second-stage controls, spatial block folds, or user-supplied nuisance learners, the package falls back to the higher-level joint residualization path in R.
}
\details{
\strong{When to use this model:}

Use \code{fit_mm_dml()} when you want orthogonalized treatment-effect estimation
but do not want the final stage to be ordinary least squares. The MM second stage
is especially useful when the residualized treatment-outcome relationship may still
contain outliers or heavy tails.

\strong{Choosing the interface:}

\itemize{
\item Use the matrix interface for simple low-level work or when you are binding to the core API from another system.
\item Use the data-frame interface for practical applications, especially when you need multiple treatments, second-stage controls, or spatial block folds.
}

\strong{Choosing learners:}

\itemize{
\item Built-in names are convenient when you want a standard workflow.
\item Use \code{make_custom_learner()} when you want to plug in your own nuisance model, such as a neural net, custom boosted model, or another pipeline.
\item The custom learner hook is the preferred extension path rather than adding a new hard-coded learner name for every model family.
}

\strong{Choosing uncertainty estimates:}

\itemize{
\item \code{se_type = "analytic"} is faster and fine for quick diagnostics.
\item \code{se_type = "bootstrap"} is the better empirical default because it propagates nuisance fitting, cross-fitting, and the final MM stage together.
}

\strong{What to inspect after fitting:}

\itemize{
\item \code{summary(fit)} for coefficients, uncertainty, fold setup, backend, and convergence
\item \code{fit\$bootstrap_estimates} when using bootstrap inference
\item \code{plot_marginal_effect()} when your second stage contains interactions
}

Built-in learner names currently include \code{"lasso"}, \code{"elastic_net"}, \code{"random_forest"}, \code{"hist_gradient_boosting"}, \code{"ols"}, \code{"xgboost"}, \code{"gam"}, and \code{"spatial_gam"}. Some of these rely on optional packages that are loaded only when needed.

Native C++ \code{"lasso"} supports richer controls through \code{learner_params}, including:
\itemize{
\item \code{cv_folds}
\item \code{n_lambda}
\item \code{lambda_min_ratio}
\item \code{lambda_grid}
\item \code{use_lambda_1se}
\item \code{standardize}
\item \code{max_iter}
\item \code{tolerance}
}

Native C++ \code{"random_forest"} supports:
\itemize{
\item \code{n_estimators}
\item \code{max_depth}
\item \code{min_samples_split}
\item \code{min_samples_leaf}
\item \code{max_features}
\item \code{max_features_fraction}
\item \code{sample_fraction}
\item \code{bootstrap}
\item \code{replacement}
\item \code{split_candidates}
\item \code{compute_oob}
\item \code{compute_importance}
}

For native fits, \code{summary(fit)} now prints learner diagnostics. The raw values are also available in \code{fit$learner_details}, including selected lambda / nonzero counts for lasso and OOB error / feature importance for random forest.

For custom nuisance models, use \code{make_custom_learner()} and provide \code{model_y} and \code{model_t} objects as either:
\enumerate{
\item prediction functions of the form \code{function(x_train, y_train, x_test, seed, role) ...}
\item lists with \code{fit} and \code{predict} functions
}

The higher-level data-frame interface supports:
\itemize{
\item multiple treatment columns
\item direct second-stage controls
\item spatial block fold assignment
\item marginal-effect plotting with \code{plot_marginal_effect()}
}
}
\examples{
set.seed(42)
n <- 250
p <- 6
x <- matrix(rnorm(n * p), nrow = n, ncol = p)
d <- 0.9 * x[, 1] - 0.4 * x[, 2] + rnorm(n, sd = 0.35)
y <- 1.35 * d + 0.5 * x[, 1] - 0.2 * x[, 3] + rnorm(n, sd = 0.30)

toy_df <- data.frame(y = y, d = d, x)
names(toy_df) <- c("y", "d", paste0("x", 1:6))

fit <- fit_mm_dml(
  outcome = "y",
  treatment = "d",
  controls = paste0("x", 1:6),
  data = toy_df,
  learner = "lasso",
  se_type = "bootstrap",
  bootstrap_replications = 50,
  n_cores = 2L,
  folds = 3,
  seed = 101
)

summary(fit)

custom_lm <- make_custom_learner(
  model_y = list(
    fit = function(x_train, y_train, seed, role) stats::lm(y_train ~ ., data = x_train),
    predict = function(model, x_test, seed, role) stats::predict(model, newdata = x_test)
  ),
  model_t = list(
    fit = function(x_train, y_train, seed, role) stats::lm(y_train ~ ., data = x_train),
    predict = function(model, x_test, seed, role) stats::predict(model, newdata = x_test)
  )
)

fit_custom <- fit_mm_dml(
  outcome = "y",
  treatment = "d",
  controls = paste0("x", 1:6),
  data = toy_df,
  learner = custom_lm,
  folds = 3,
  seed = 123
)

summary(fit_custom)
}
