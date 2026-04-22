fit_s_estimator <- function(x,
                            y = NULL,
                            data = NULL,
                            n_starts = 500L,
                            n_best_starts = 25L,
                            max_refine = 100L,
                            max_scale_iter = 100L,
                            tol = 1e-8,
                            scale_tol = 1e-10,
                            c = 1.54764,
                            b = 0.5,
                            add_intercept = TRUE,
                            ridge = 1e-10,
                            min_weight = 1e-12,
                            seed = 123456789,
                            use_fast_s = TRUE,
                            include_ols_start = TRUE,
                            fast_s_screen_subsets = 500L,
                            fast_s_screen_iters = 2L,
                            fast_s_keep = 25L) {
  inputs <- .as_matrix_inputs(x, y, data, add_intercept = add_intercept)

  fit <- .Call(
    rc_r_fit_s_estimator,
    inputs$X,
    inputs$y,
    inputs$core_add_intercept,
    as.integer(n_starts),
    as.integer(n_best_starts),
    as.integer(max_refine),
    as.integer(max_scale_iter),
    as.numeric(tol),
    as.numeric(scale_tol),
    as.numeric(c),
    as.numeric(b),
    as.numeric(ridge),
    as.numeric(min_weight),
    as.numeric(seed),
    isTRUE(use_fast_s),
    isTRUE(include_ols_start),
    as.integer(fast_s_screen_subsets),
    as.integer(fast_s_screen_iters),
    as.integer(fast_s_keep)
  )

  fit$coef_names <- .finalize_coef_names(inputs$coef_names, inputs$core_add_intercept, length(fit$coef))
  fit$formula <- inputs$formula
  fit$call <- match.call()
  class(fit) <- "robustcause_s"
  fit
}

print.robustcause_s <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

summary.robustcause_s <- function(object, ...) {
  coef_table <- cbind(Estimate = object$coef)
  rownames(coef_table) <- object$coef_names

  structure(
    list(
      call = object$call,
      formula = object$formula,
      converged = object$converged,
      iterations = object$iterations,
      scale = object$scale,
      objective = object$objective,
      starts_tried = object$starts_tried,
      starts_used = object$starts_used,
      message = object$message,
      coefficients = coef_table
    ),
    class = "summary.robustcause_s"
  )
}

print.summary.robustcause_s <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("RobustCause S-estimator fit\n")
  if (!is.null(x$formula)) {
    cat("Formula:", x$formula, "\n")
  }
  cat("Converged:", x$converged, "\n")
  cat("Iterations:", x$iterations, "\n")
  cat("Scale:", format(signif(x$scale, 6)), "\n")
  cat("Starts tried:", x$starts_tried, "\n")
  cat("Starts used:", x$starts_used, "\n")
  if (!is.null(x$message) && nzchar(x$message)) {
    cat("Message:", x$message, "\n")
  }
  cat("\nCoefficients:\n")
  print(round(x$coefficients, digits))
  invisible(x)
}

plot_weights <- function(object, ...) {
  weights <- object$weights
  if (is.null(weights)) {
    stop("`object` does not contain weights.", call. = FALSE)
  }
  graphics::plot(
    seq_along(weights),
    weights,
    xlab = "Observation",
    ylab = "Weight",
    main = "Robust observation weights",
    ...
  )
}
