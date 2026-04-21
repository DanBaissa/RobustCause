.as_matrix_inputs <- function(x, y = NULL, data = NULL, add_intercept = TRUE) {
  if (inherits(x, "formula")) {
    if (is.null(data)) {
      stop("`data` is required when `x` is a formula.", call. = FALSE)
    }
    mf <- model.frame(x, data = data, na.action = na.omit)
    response <- model.response(mf)
    design <- model.matrix(attr(mf, "terms"), data = mf)
    if (!add_intercept && "(Intercept)" %in% colnames(design)) {
      design <- design[, colnames(design) != "(Intercept)", drop = FALSE]
    }
    return(list(
      X = unname(as.matrix(design)),
      y = as.numeric(response),
      coef_names = colnames(design),
      formula = deparse(x),
      core_add_intercept = FALSE
    ))
  }

  if (is.null(y)) {
    stop("`y` is required when `x` is not a formula.", call. = FALSE)
  }

  X <- as.matrix(x)
  if (!is.numeric(X)) {
    stop("`x` must coerce to a numeric matrix.", call. = FALSE)
  }

  list(
    X = unname(X),
    y = as.numeric(y),
    coef_names = colnames(X),
    formula = NULL,
    core_add_intercept = isTRUE(add_intercept)
  )
}

.psi_name <- function(psi) {
  psi <- match.arg(psi, c("huber", "tukey_bisquare"))
  psi
}

.hc_name <- function(hc_type) {
  match.arg(hc_type, c("HC0", "HC1", "HC2", "HC3", "HC4", "HC4m", "HC5"))
}

.finalize_coef_names <- function(coef_names, add_intercept, n_coef) {
  if (is.null(coef_names)) {
    coef_names <- character()
  }

  coef_names <- as.character(coef_names)
  if (isTRUE(add_intercept) && length(coef_names) == n_coef - 1L) {
    coef_names <- c("(Intercept)", coef_names)
  }
  if (length(coef_names) != n_coef) {
    coef_names <- paste0("b", seq_len(n_coef))
  }
  coef_names
}

fit_rlm <- function(x,
                    y = NULL,
                    data = NULL,
                    psi = c("huber", "tukey_bisquare"),
                    tuning = 1.345,
                    maxit = 100L,
                    tol = 1e-8,
                    add_intercept = TRUE,
                    ridge = 1e-10,
                    min_weight = 1e-12,
                    hc_type = "HC3",
                    zcrit = qnorm(0.975)) {
  psi <- .psi_name(psi)
  hc_type <- .hc_name(hc_type)
  inputs <- .as_matrix_inputs(x, y, data, add_intercept = add_intercept)

  fit <- .Call(
    rc_r_fit_rlm,
    inputs$X,
    inputs$y,
    psi,
    as.numeric(tuning),
    as.integer(maxit),
    as.numeric(tol),
    inputs$core_add_intercept,
    as.numeric(ridge),
    as.numeric(min_weight)
  )

  inference <- .Call(
    rc_r_confint_rlm,
    inputs$X,
    inputs$y,
    psi,
    as.numeric(tuning),
    as.integer(maxit),
    as.numeric(tol),
    inputs$core_add_intercept,
    as.numeric(ridge),
    as.numeric(min_weight),
    hc_type,
    as.numeric(zcrit)
  )

  fit$vcov <- inference$vcov
  fit$se <- inference$se
  fit$ci <- inference$ci
  fit$coef_names <- .finalize_coef_names(inputs$coef_names, inputs$core_add_intercept, length(fit$coef))
  fit$formula <- inputs$formula
  fit$call <- match.call()
  fit$psi <- psi
  fit$hc_type <- hc_type
  class(fit) <- "robustcause_rlm"
  fit
}

print.robustcause_rlm <- function(x, ...) {
  cat("RobustCause RLM fit\n")
  if (!is.null(x$formula)) {
    cat("Formula:", x$formula, "\n")
  }
  cat("Psi:", x$psi, "\n")
  cat("HC type:", x$hc_type, "\n")
  cat("Converged:", x$converged, "\n")
  cat("Iterations:", x$iterations, "\n")
  cat("Scale:", format(signif(x$scale, 6)), "\n")
  cat("\nCoefficients:\n")
  print(stats::setNames(x$coef, x$coef_names))
  invisible(x)
}

vcov.robustcause_rlm <- function(object, ...) {
  object$vcov
}

confint.robustcause_rlm <- function(object, ...) {
  object$ci
}

`%||%` <- function(lhs, rhs) {
  if (is.null(lhs)) rhs else lhs
}
