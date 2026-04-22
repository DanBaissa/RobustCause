\name{fit_mmm}
\alias{fit_mmm}
\title{Fit a Robust Marketing Mix Model}
\usage{
fit_mmm(
  y = NULL,
  unit_ids = NULL,
  channels = NULL,
  preclean_x = NULL,
  controls = NULL,
  data = NULL,
  outcome = NULL,
  unit = NULL,
  channel_cols = NULL,
  preclean_cols = NULL,
  control_cols = NULL,
  robustness = NULL,
  fit_method = c("huber", "ols", "s", "mm"),
  rho = 0.85,
  increment_method = "adaptive_clip",
  hill_lambda = 0.10,
  preclean = TRUE,
  preclean_method = "mm",
  nonnegative = TRUE,
  residual_clip_multiplier = 2.5,
  clip_c = 1.5,
  adaptive_upper = 0,
  apply_hill = TRUE,
  rho_grid = NULL,
  hill_grid = NULL,
  choose_by = c("bic", "aic", "mse"),
  max_grid_size = 2500L,
  huber_c = 1.345,
  huber_max_iter = 75L,
  s_starts = 20L,
  s_max_iter = 50L,
  mm_max_iter = 50L,
  mm_tukey_c = 4.685,
  seed = 1234L
)
}
\description{
Fits a multi-channel marketing mix model with optional first-step channel
cleaning, adstock carryover, Hill saturation, robust final regression, and
optional grid search over channel-level carryover and saturation parameters.
}
\value{
An object of class \code{robustcause_mmm} containing the fitted coefficients,
per-channel cleaned signals, stocks, transformed regressors, the final design
matrix used in the regression, and optionally a model-search table when grid
search is used.
}
\arguments{
  \item{y}{Numeric outcome vector for the matrix interface.}
  \item{unit_ids}{Unit or panel identifiers for the matrix interface.}
  \item{channels}{Numeric channel matrix, one column per media input.}
  \item{preclean_x}{Numeric design matrix used in the first-step cleaning stage.}
  \item{controls}{Numeric control matrix for the final MMM regression.}
  \item{data}{Optional data frame interface.}
  \item{outcome}{Outcome column name in \code{data}.}
  \item{unit}{Unit or panel identifier column name in \code{data}.}
  \item{channel_cols}{Media channel columns in \code{data}.}
  \item{preclean_cols}{Columns used in the channel-cleaning stage.}
  \item{control_cols}{Columns used in the final regression stage.}
  \item{robustness}{Optional high-level profile: \code{"classical"}, \code{"huber"}, \code{"s"}, \code{"mm"}, or \code{"custom"}. When supplied, this sets coherent defaults for the cleaning stage, adstock increment rule, and final regression.}
  \item{fit_method}{Final regression type: \code{"ols"}, \code{"huber"}, \code{"s"}, or \code{"mm"}.}
  \item{rho}{Carryover parameter, scalar or one value per channel.}
  \item{increment_method}{Adstock increment rule, scalar or one value per channel.}
  \item{hill_lambda}{Hill saturation parameter, scalar or one value per channel.}
  \item{preclean}{Whether to clean each channel before adstock accumulation, scalar or one value per channel.}
  \item{preclean_method}{Either \code{"mm"} or \code{"s"}, scalar or one value per channel.}
  \item{nonnegative}{Whether to truncate cleaned channels at zero.}
  \item{residual_clip_multiplier}{Residual clipping multiplier used in channel precleaning.}
  \item{clip_c}{Bounding constant for robust increment rules.}
  \item{adaptive_upper}{Optional threshold for \code{"adaptive_clip"}.}
  \item{apply_hill}{Whether to apply Hill saturation, scalar or one value per channel.}
  \item{rho_grid}{Optional shared or per-channel grid for searching over carryover parameters. Supply either a numeric vector shared across channels or a list with one numeric vector per channel.}
  \item{hill_grid}{Optional shared or per-channel grid for searching over Hill saturation parameters. Supply either a numeric vector shared across channels or a list with one numeric vector per channel.}
  \item{choose_by}{Criterion used to rank search candidates.}
  \item{max_grid_size}{Upper bound on the total number of search candidates.}
  \item{huber_c}{Huber tuning constant for the final Huber regression.}
  \item{huber_max_iter}{Maximum Huber iterations.}
  \item{s_starts}{Number of S-estimator starts.}
  \item{s_max_iter}{Maximum S iterations.}
  \item{mm_max_iter}{Maximum MM iterations.}
  \item{mm_tukey_c}{Tukey bisquare constant for MM stages.}
  \item{seed}{Random seed.}
}
\details{
\strong{What this model does:}

For each media channel, \code{fit_mmm()} can:

\enumerate{
\item robustly preclean the raw signal
\item build a carryover stock
\item optionally apply Hill saturation
\item fit a final regression of the outcome on transformed media and controls
}

\strong{When to use this model:}

Use \code{fit_mmm()} when you want a full media-response pipeline rather than
just a stock series. It is the right entry point when you care about channel
coefficients, transformed regressors, and the final design actually used in the
MMM.

\strong{Fixed specification vs. model search:}

\itemize{
\item If you supply fixed values for \code{rho} and \code{hill_lambda}, \code{fit_mmm()} fits that exact MMM specification.
\item If you supply \code{rho_grid} and/or \code{hill_grid}, the function evaluates the candidate MMM specifications and chooses the best one by \code{BIC}, \code{AIC}, or \code{MSE}.
\item Search grids can be shared across channels or specified per channel.
}

\strong{Practical defaults:}

\itemize{
\item Use \code{robustness = "classical"} for a plain baseline MMM.
\item Start with \code{fit_method = "huber"} or \code{"mm"}.
\item Use \code{preclean = TRUE} when raw channel spikes are plausible.
\item Use \code{increment_method = "adaptive_clip"} or \code{"huber"} as robust first choices.
\item Keep \code{apply_hill = TRUE} if you want diminishing-returns behavior.
}

\strong{What to inspect after fitting:}

\itemize{
\item \code{summary(fit)} for channel-level diagnostics and the final coefficient table
\item \code{fit\$cleaned_signals} to see the precleaned channel series
\item \code{fit\$stocks} for carryover stocks
\item \code{fit\$transformed_channels} for the regressors actually entering the final regression
\item \code{fit\$design_matrix} for the full final design
\item \code{fit\$grid_results} when using model search
}
}
\examples{
set.seed(7)
n <- 120
unit_ids <- rep(1:24, each = 5)
x1 <- rnorm(n)
x2 <- rnorm(n)
channel_1 <- exp(0.5 + 0.4 * x1 + rnorm(n, sd = 0.15))
channel_2 <- exp(0.3 + 0.3 * x2 + rnorm(n, sd = 0.15))
y <- 1 + 0.05 * channel_1 + 0.03 * channel_2 + 0.4 * x1 - 0.2 * x2 + rnorm(n, sd = 0.4)

mmm_dat <- data.frame(
  y = y,
  unit = unit_ids,
  x1 = x1,
  x2 = x2,
  channel_1 = channel_1,
  channel_2 = channel_2
)

fit <- fit_mmm(
  data = mmm_dat,
  outcome = "y",
  unit = "unit",
  channel_cols = c("channel_1", "channel_2"),
  preclean_cols = c("x1", "x2"),
  control_cols = c("x1", "x2"),
  robustness = "mm",
  rho_grid = list(c(0.6, 0.8), c(0.5, 0.75)),
  hill_grid = list(c(0.05, 0.12), c(0.05, 0.10)),
  choose_by = "bic"
)

summary(fit)
}
