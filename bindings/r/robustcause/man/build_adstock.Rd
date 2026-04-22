\name{build_adstock}
\alias{build_adstock}
\title{Build Robust Adstock Series}
\usage{
build_adstock(
  signal = NULL,
  unit_ids = NULL,
  x_preclean = NULL,
  data = NULL,
  signal_col = NULL,
  unit_col = NULL,
  controls = NULL,
  preclean = FALSE,
  preclean_method = c("mm", "s"),
  nonnegative = FALSE,
  residual_clip_multiplier = 2.5,
  s_starts = 20L,
  s_max_iter = 50L,
  mm_max_iter = 50L,
  mm_tukey_c = 4.685,
  seed = 1234L,
  rho = 0.85,
  increment_method = c("plain", "huber", "tanh", "softsign", "adaptive_clip"),
  clip_c = 1.5,
  adaptive_upper = 0
)
}
\arguments{
  \item{signal}{Numeric signal vector.}
  \item{unit_ids}{Unit or panel identifiers aligned with \code{signal}.}
  \item{x_preclean}{Numeric design matrix for the preclean step.}
  \item{data}{Optional data frame interface.}
  \item{signal_col}{Signal column name in \code{data}.}
  \item{unit_col}{Unit-id column name in \code{data}.}
  \item{controls}{Optional control columns used for precleaning.}
  \item{preclean}{Whether to run first-step robust cleaning before adstock accumulation.}
  \item{preclean_method}{Either \code{"mm"} or \code{"s"}.}
  \item{nonnegative}{Whether to truncate the cleaned signal at zero.}
  \item{residual_clip_multiplier}{Residual clipping multiplier applied around the fitted preclean mean.}
  \item{s_starts}{Number of random S-estimator starts.}
  \item{s_max_iter}{Maximum local S-refinement iterations.}
  \item{mm_max_iter}{Maximum MM-refinement iterations.}
  \item{mm_tukey_c}{Tukey bisquare tuning constant for the MM step.}
  \item{seed}{Random seed used inside the robust preclean routine.}
  \item{rho}{Adstock carryover parameter.}
  \item{increment_method}{One of \code{"plain"}, \code{"huber"}, \code{"tanh"}, \code{"softsign"}, or \code{"adaptive_clip"}.}
  \item{clip_c}{Clipping or bounding constant for robust increment rules.}
  \item{adaptive_upper}{Optional upper threshold for \code{"adaptive_clip"}; if non-positive it is estimated from the signal.}
}
\value{
An object of class \code{robustcause_adstock} containing the final stock series,
the cleaned signal used to build it, and metadata about the chosen preclean and
increment rules.
}
\description{
Builds a plain or robust adstock series with optional first-step S- or MM-based
signal cleaning. The interface supports both direct vector input and a data-frame
workflow.
}
\details{
\strong{When to use this model:}

Use \code{build_adstock()} when you want a carryover stock but do not yet want to
fit a full marketing mix model. It is useful for:

\itemize{
\item inspecting how different carryover rules transform a raw media series
\item pre-building a stock variable for downstream causal or regression models
\item stress-testing whether outliers in the raw signal are driving your carryover measure
}

\strong{Precleaning:}

\itemize{
\item Set \code{preclean = TRUE} when the raw signal may contain spikes, recording errors, or one-off bursts you do not want to feed directly into recursive stock accumulation.
\item Use \code{preclean_method = "mm"} as the usual robust default.
\item Use \code{preclean_method = "s"} when you want the more aggressively high-breakdown first step itself.
}

\strong{Increment rules:}

\itemize{
\item \code{"plain"} gives the usual non-robust recursive adstock.
\item \code{"huber"} and \code{"adaptive_clip"} are good starting points when you want robustification but still want the series to behave like a familiar carryover stock.
\item \code{"tanh"} and \code{"softsign"} give smoother bounded updates.
}

\strong{What to inspect after fitting:}

\itemize{
\item \code{summary(fit)} for signal, cleaned-signal, and stock summaries
\item \code{fit\$cleaned_signal} to see how much the preclean stage moved the raw series
\item \code{fit\$stock} as the final carryover series
}
}
\examples{
signal <- c(10, 12, 60, 9, 8, 40)
unit_ids <- c(1, 1, 1, 2, 2, 2)
x_preclean <- cbind(1, c(0, 1, 2, 0, 1, 2))

fit <- build_adstock(
  signal = signal,
  unit_ids = unit_ids,
  x_preclean = x_preclean,
  preclean = TRUE,
  preclean_method = "mm",
  increment_method = "huber"
)

fit
summary(fit)
}
