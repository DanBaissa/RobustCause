\name{make_custom_learner}
\alias{make_custom_learner}
\title{Create a Custom MM-DML Learner}
\usage{
make_custom_learner(model_y, model_t = model_y, label = "custom")
}
\arguments{
  \item{model_y}{Outcome nuisance learner specification.}
  \item{model_t}{Treatment nuisance learner specification. Defaults to \code{model_y}.}
  \item{label}{Short label shown in summaries.}
}
\value{
A custom learner object for use in \code{fit_mm_dml()}.
}
\description{
Builds a user-supplied nuisance learner specification for \code{fit_mm_dml()}.
Each learner can be either a direct prediction function or a list with \code{fit}
and \code{predict} functions.
}
\details{
\strong{When to use this helper:}

Use \code{make_custom_learner()} when you want \code{fit_mm_dml()} to use your
own nuisance model rather than one of the package's named learners. This is the
recommended route for custom neural networks, tuned gradient boosting pipelines,
or project-specific regression stacks.

Supported specification styles:
\enumerate{
\item \code{function(x_train, y_train, x_test, seed, role) ...} returning predictions for \code{x_test}
\item \code{list(fit = function(...), predict = function(model, x_test, ...))}
}

This lets users plug in arbitrary model pipelines without waiting for the package
to add a named learner wrapper.

Use this helper when you want a completely custom nuisance model. If you want to
stay on the native compiled backend for the built-in learners, prefer
\code{learner = "lasso"} or \code{learner = "random_forest"} together with
\code{learner_params = list(...)} in \code{fit_mm_dml()}.
}
\examples{
custom_lm <- make_custom_learner(
  model_y = list(
    fit = function(x_train, y_train, seed, role) stats::lm(y_train ~ ., data = x_train),
    predict = function(model, x_test, seed, role) stats::predict(model, newdata = x_test)
  ),
  model_t = list(
    fit = function(x_train, y_train, seed, role) stats::lm(y_train ~ ., data = x_train),
    predict = function(model, x_test, seed, role) stats::predict(model, newdata = x_test)
  ),
  label = "custom_lm"
)
}
