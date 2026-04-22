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
  psi_args = NULL,
  maxit = 100L,
  tol = 1e-08,
  acc = NULL,
  add_intercept = TRUE,
  ridge = 1e-10,
  min_weight = 1e-12,
  weights = NULL,
  wt.method = c("inv.var", "case"),
  init = NULL,
  scale.est = NULL,
  k2 = NULL,
  test.vec = NULL,
  subset = NULL,
  na.action = stats::na.omit,
  contrasts = NULL,
  offset = NULL,
  mm_s_control = NULL
)
}
\arguments{
  \item{x}{A numeric design matrix or a formula.}
  \item{y}{Numeric response vector when \code{x} is not a formula.}
  \item{data}{A data frame used with formula input.}
  \item{method}{Either \code{"m"} or \code{"mm"}.}
  \item{psi}{Psi function name or a psi function. Built-in names are \code{"huber"}, \code{"tukey_bisquare"}, and \code{"hampel"}. If \code{NULL}, the default is Huber for M and Tukey bisquare for MM. Named psi families are handled by the native C++ solver; function-valued psi specifications are handled in the package wrapper without delegating to MASS.}
  \item{tuning}{Psi tuning constant. If \code{NULL}, the default is \code{1.345} for M and \code{4.685} for MM.}
  \item{psi_args}{Optional named list of extra psi-function arguments for custom psi functions.}
  \item{maxit}{Maximum IRLS iterations.}
  \item{tol}{Convergence tolerance.}
  \item{acc}{Optional convergence threshold used for the IRLS stopping rule.}
  \item{add_intercept}{Whether to add an intercept for matrix input.}
  \item{ridge}{Small ridge term used in linear solves.}
  \item{min_weight}{Lower bound for weights in the solver.}
  \item{weights}{Optional observation weights.}
  \item{wt.method}{Weight interpretation for weighted fits: \code{"inv.var"} for inverse-variance style weights or \code{"case"} for case weights.}
  \item{init}{Optional starting specification for M-estimation. Supported values are \code{"ls"}, \code{"lts"}, a numeric coefficient vector, or a list containing \code{$coef}.}
  \item{scale.est}{Optional scale estimator for M-estimation, such as \code{"MAD"}, \code{"Huber"}, or \code{"proposal 2"}.}
  \item{k2}{Optional tuning constant used by some MASS scale estimators.}
  \item{test.vec}{Optional convergence target for the IRLS iterations. Common values are \code{"resid"}, \code{"coef"}, and \code{"w"}.}
  \item{subset}{Optional subset expression for formula fits.}
  \item{na.action}{Missing-data handling for formula fits.}
  \item{contrasts}{Optional contrasts list passed to \code{model.matrix()} for formula fits.}
  \item{offset}{Optional offset for formula fits.}
  \item{mm_s_control}{Optional named list of S-estimator controls used to initialize MM fits.}
}
\value{
An object of class \code{robustcause_rlm}.
}
\description{
Fits the RobustCause robust linear model using either a plain M-estimator or an MM-estimator seeded by the S-estimator. Formula parsing stays in R, but the named robust-regression controls are handled natively by RobustCause rather than delegated to \code{MASS::rlm()}.
}
\details{
\strong{When to use this model:}

\itemize{
\item Use \code{method = "m"} when you want a robust regression that downweights moderate outliers but stays fairly close to ordinary robust IRLS workflows.
\item Use \code{method = "mm"} when you want higher breakdown protection. The MM path starts from an S-estimator and is the safer default when contamination may be severe.
}

\strong{Choosing \code{psi} and \code{tuning}:}

\itemize{
\item For \code{method = "m"}, the default is Huber with \code{tuning = 1.345}. This is a good general-purpose choice for mild to moderate contamination.
\item For \code{method = "mm"}, the default is Tukey bisquare with \code{tuning = 4.685}. This is the standard high-efficiency MM default and is usually the right starting point.
\item Hampel psi is available through the MASS-compatible backend when you want a redescending alternative with explicit shoulder and rejection points.
\item If you supply a psi function directly, use \code{psi_args} for any additional constants that function needs.
}

\strong{Implementation split:}

\itemize{
\item Formula handling features such as \code{subset}, \code{na.action}, \code{contrasts}, and \code{offset} are resolved in the R wrapper before the model matrix is handed to the solver.
\item Named psi families, weighted fits, initialization choices, scale-estimation choices, and convergence controls are handled by the RobustCause implementation.
\item Function-valued psi inputs are supported through the package wrapper without switching to \code{MASS::rlm()}.
}

\strong{Inference:}

Fitting and inference are intentionally separate. After estimating the model with
\code{fit_rlm()}, use \code{vcov()}, \code{confint()}, \code{vcov_robust()}, or
\code{confint_robust()} to choose the heteroskedasticity-consistent covariance
estimator. In practice, \code{HC3} is a good default when sample size is not large.

\strong{What to inspect after fitting:}

\itemize{
\item \code{summary(fit)} for coefficients, robust standard errors, and convergence
\item \code{fit\$weights} to see which observations were downweighted
\item \code{plot_weights(fit)} for a quick weight diagnostic
}
}
\examples{
set.seed(42)
n <- 200
x1 <- rnorm(n)
x2 <- rnorm(n)
y <- 1 + 2 * x1 - 1.5 * x2 + rnorm(n, sd = 0.3)
y[c(10, 40, 150)] <- y[c(10, 40, 150)] + c(6, -7, 5)

fit_m <- fit_rlm(y ~ x1 + x2, data = data.frame(y, x1, x2), method = "m")
summary(fit_m)

fit_mm <- fit_rlm(y ~ x1 + x2, data = data.frame(y, x1, x2), method = "mm")
summary(fit_mm)
confint(fit_mm, hc_type = "HC3")

w <- runif(n, 0.8, 1.2)
off <- 0.25 * x2
keep <- x1 > -1.5
fit_mass <- fit_rlm(
  y ~ x1 + x2,
  data = data.frame(y, x1, x2, w, off, keep),
  method = "m",
  psi = "hampel",
  weights = w,
  wt.method = "case",
  subset = keep,
  offset = off,
  scale.est = "Huber",
  test.vec = "coef"
)
summary(fit_mass)
}
