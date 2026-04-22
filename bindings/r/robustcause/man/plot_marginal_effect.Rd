\name{plot_marginal_effect}
\alias{plot_marginal_effect}
\title{Plot MM-DML Marginal Effects}
\usage{
plot_marginal_effect(
  fit,
  focal,
  interaction,
  moderator_values = c(0, 1),
  moderator_labels = NULL,
  ci_level = 0.95,
  xlab = NULL,
  ylab = NULL,
  main = NULL,
  pch = 19,
  col = "black",
  ci_col = "black",
  ci_lty = 1,
  ci_lwd = 1.2,
  ...
)
}
\arguments{
  \item{fit}{An object returned by \code{fit_mm_dml()}.}
  \item{focal}{Coefficient name for the focal treatment effect.}
  \item{interaction}{Coefficient name for the interaction term.}
  \item{moderator_values}{Moderator values at which to evaluate the marginal effect.}
  \item{moderator_labels}{Optional plot labels for the moderator values.}
  \item{ci_level}{Confidence interval level.}
  \item{xlab}{Optional x-axis label.}
  \item{ylab}{Optional y-axis label.}
  \item{main}{Optional plot title.}
  \item{pch}{Point character.}
  \item{col}{Point color.}
  \item{ci_col}{Confidence interval color.}
  \item{ci_lty}{Confidence interval line type.}
  \item{ci_lwd}{Confidence interval line width.}
  \item{...}{Additional graphical arguments passed to \code{plot()}.}
}
\value{
An invisible data frame containing the marginal effect estimates and intervals.
}
\description{
Plots the marginal effect implied by a focal coefficient and an interaction
coefficient from a fitted \code{fit_mm_dml()} model.
}
