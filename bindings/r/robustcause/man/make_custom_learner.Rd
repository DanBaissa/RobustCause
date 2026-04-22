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
Supported specification styles:
\enumerate{
\item \code{function(x_train, y_train, x_test, seed, role) ...} returning predictions for \code{x_test}
\item \code{list(fit = function(...), predict = function(model, x_test, ...))}
}

This lets users plug in arbitrary model pipelines without waiting for the package
to add a named learner wrapper.
}
