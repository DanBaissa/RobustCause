\name{fit_adstock_model}
\alias{fit_adstock_model}
\alias{print.robustcause_adstock_model}
\alias{summary.robustcause_adstock_model}
\alias{print.summary.robustcause_adstock_model}
\title{Fit an End-to-End Adstock Response Model}
\usage{
fit_adstock_model(
  y = NULL,
  signal = NULL,
  unit_ids = NULL,
  controls = NULL,
  data = NULL,
  outcome = NULL,
  signal_col = NULL,
  unit_col = NULL,
  control_cols = NULL,
  robustness = c("classical", "huber", "s", "mm"),
  preclean_method = NULL,
  increment_method = NULL,
  fit_method = NULL,
  rho_grid = seq(0.1, 0.9, by = 0.05),
  hill_grid = c(0, 0.05, 0.1, 0.2),
  choose_by = c("bic", "aic", "mse"),
  nonnegative = FALSE,
  residual_clip_multiplier = 2.5,
  clip_c = 1.5,
  adaptive_upper = 0,
  s_starts = 20L,
  s_max_iter = 50L,
  mm_max_iter = 50L,
  mm_tukey_c = 4.685,
  huber_c = 1.345,
  huber_max_iter = 75L,
  seed = 1234L
)
}
\arguments{
  \item{y}{Numeric outcome vector for the matrix interface.}
  \item{signal}{Numeric exposure or media signal to be transformed into an adstock regressor.}
  \item{unit_ids}{Unit or panel identifiers aligned with \code{signal}.}
  \item{controls}{Optional numeric control matrix.}
  \item{data}{Optional data-frame interface.}
  \item{outcome}{Outcome column name in \code{data}.}
  \item{signal_col}{Signal column name in \code{data}.}
  \item{unit_col}{Unit-id column name in \code{data}.}
  \item{control_cols}{Optional control columns in \code{data}.}
  \item{robustness}{One of \code{"classical"}, \code{"huber"}, \code{"s"}, or \code{"mm"}. This selects a full default adstock-modeling level.}
  \item{preclean_method}{Optional override for the signal-cleaning stage. Use \code{"none"}, \code{"s"}, or \code{"mm"}.}
  \item{increment_method}{Optional override for the adstock increment rule.}
  \item{fit_method}{Optional override for the final regression stage.}
  \item{rho_grid}{Candidate carryover values to search over.}
  \item{hill_grid}{Candidate Hill-saturation values to search over. Include \code{0} to allow no saturation.}
  \item{choose_by}{Selection criterion used to choose the best candidate from the search grid.}
  \item{nonnegative}{Whether to truncate the cleaned signal at zero.}
  \item{residual_clip_multiplier}{Residual clipping multiplier used in robust precleaning.}
  \item{clip_c}{Bounding constant for robust increment rules.}
  \item{adaptive_upper}{Optional threshold for \code{"adaptive_clip"}.}
  \item{s_starts}{Number of S-estimator starts used in robust stages.}
  \item{s_max_iter}{Maximum local S iterations.}
  \item{mm_max_iter}{Maximum MM iterations.}
  \item{mm_tukey_c}{Tukey bisquare constant used in MM stages.}
  \item{huber_c}{Huber tuning constant.}
  \item{huber_max_iter}{Maximum Huber IRLS iterations.}
  \item{seed}{Random seed for robust stages.}
}
\value{
An object of class \code{robustcause_adstock_model} containing the selected
adstock parameters, the chosen stock and transformed stock, the final fitted
regression object, and the full grid-search table.
}
\description{
Fits a complete adstock response model end to end. The function searches over a
grid of adstock carryover and Hill-saturation parameters, builds the stock
series, fits the downstream regression, and returns the best specification.
}
\details{
\strong{What this function adds beyond \code{build_adstock()}:}

\itemize{
\item \code{build_adstock()} only transforms a signal into a carryover stock.
\item \code{fit_adstock_model()} adds model fitting and parameter selection on top of that transformation.
\item This is the standalone entry point when you want a full single-signal adstock model before moving up to MMM.
}

\strong{Robustness levels:}

\itemize{
\item \code{"classical"} uses no precleaning, \code{"plain"} increments, and an OLS final regression.
\item \code{"huber"} uses Huber-style defaults for the transformation and final regression.
\item \code{"s"} uses S-based cleaning and an S final regression.
\item \code{"mm"} uses MM-based cleaning and an MM final regression.
}

You can override the defaults by supplying \code{preclean_method},
\code{increment_method}, or \code{fit_method} directly.

\strong{Model search:}

The function evaluates each \code{rho} and \code{hill_lambda} combination,
fits the final regression, and ranks candidates by \code{BIC}, \code{AIC}, or
\code{MSE}. The top grid candidates are included in \code{summary(fit)}.

\strong{What to inspect after fitting:}

\itemize{
\item \code{summary(fit)} for the selected specification and top candidates
\item \code{fit\$stock} for the selected carryover series
\item \code{fit\$transformed_stock} for the final regressor entering the model
\item \code{fit\$grid_results} to inspect the entire search path
\item \code{fit\$final_fit} for the underlying regression object
}
}
\examples{
set.seed(11)
n <- 180
unit_ids <- rep(1:36, each = 5)
x1 <- rnorm(n)
signal <- exp(0.4 + 0.3 * x1 + rnorm(n, sd = 0.2))
y <- 1 + 0.08 * signal + 0.5 * x1 + rnorm(n, sd = 0.3)

dat <- data.frame(
  y = y,
  signal = signal,
  unit = unit_ids,
  x1 = x1
)

fit <- fit_adstock_model(
  data = dat,
  outcome = "y",
  signal_col = "signal",
  unit_col = "unit",
  control_cols = "x1",
  robustness = "classical"
)

summary(fit)
}
