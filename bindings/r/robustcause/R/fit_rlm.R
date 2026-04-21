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

.rlm_method_name <- function(method) {
  match.arg(method, c("m", "mm"))
}

.psi_name <- function(psi, method) {
  if (is.null(psi)) {
    return(if (identical(method, "mm")) "tukey_bisquare" else "huber")
  }
  match.arg(psi, c("huber", "tukey_bisquare"))
}

.hc_name <- function(hc_type) {
  match.arg(hc_type, c("HC0", "HC1", "HC2", "HC3", "HC4", "HC4m", "HC5"))
}

.default_mm_s_control <- function() {
  list(
    n_starts = 500L,
    n_best_starts = 25L,
    max_refine = 100L,
    max_scale_iter = 100L,
    tol = 1e-8,
    scale_tol = 1e-10,
    c = 1.54764,
    b = 0.5,
    ridge = 1e-10,
    min_weight = 1e-12,
    seed = 123456789,
    use_fast_s = TRUE,
    include_ols_start = TRUE,
    fast_s_screen_subsets = 500L,
    fast_s_screen_iters = 2L,
    fast_s_keep = 25L
  )
}

.coerce_mm_s_control <- function(mm_s_control) {
  defaults <- .default_mm_s_control()
  if (is.null(mm_s_control)) {
    return(defaults)
  }

  unknown <- setdiff(names(mm_s_control), names(defaults))
  if (length(unknown) > 0L) {
    stop("Unknown `mm_s_control` fields: ", paste(unknown, collapse = ", "), call. = FALSE)
  }

  modifyList(defaults, mm_s_control)
}

.make_rlm_control <- function(method,
                              psi,
                              tuning,
                              maxit,
                              tol,
                              add_intercept,
                              ridge,
                              min_weight,
                              mm_s_control = NULL) {
  list(
    method = .rlm_method_name(method),
    psi = .psi_name(psi, method),
    tuning = as.numeric(tuning),
    maxit = as.integer(maxit),
    tol = as.numeric(tol),
    add_intercept = isTRUE(add_intercept),
    ridge = as.numeric(ridge),
    min_weight = as.numeric(min_weight),
    mm_s_control = .coerce_mm_s_control(mm_s_control)
  )
}

.fit_rlm_from_inputs <- function(inputs, control) {
  fit <- .Call(rc_r_fit_rlm, inputs$X, inputs$y, control)
  fit$coef_names <- .finalize_coef_names(inputs$coef_names, inputs$core_add_intercept, length(fit$coef))
  fit$formula <- inputs$formula
  fit$control <- control
  fit$X_input <- inputs$X
  fit$y_input <- inputs$y
  class(fit) <- "robustcause_rlm"
  fit
}

fit_rlm <- function(x,
                    y = NULL,
                    data = NULL,
                    method = c("m", "mm"),
                    psi = NULL,
                    tuning = 1.345,
                    maxit = 100L,
                    tol = 1e-8,
                    add_intercept = TRUE,
                    ridge = 1e-10,
                    min_weight = 1e-12,
                    mm_s_control = NULL) {
  method <- .rlm_method_name(method)
  inputs <- .as_matrix_inputs(x, y, data, add_intercept = add_intercept)
  control <- .make_rlm_control(
    method = method,
    psi = psi,
    tuning = tuning,
    maxit = maxit,
    tol = tol,
    add_intercept = inputs$core_add_intercept,
    ridge = ridge,
    min_weight = min_weight,
    mm_s_control = mm_s_control
  )

  fit <- .fit_rlm_from_inputs(inputs, control)
  fit$call <- match.call()
  fit
}

print.robustcause_rlm <- function(x, ...) {
  cat("RobustCause RLM fit\n")
  if (!is.null(x$formula)) {
    cat("Formula:", x$formula, "\n")
  }
  cat("Method:", x$method, "\n")
  cat("Psi:", x$psi, "\n")
  cat("Converged:", x$converged, "\n")
  cat("Iterations:", x$iterations, "\n")
  cat("Scale:", format(signif(x$scale, 6)), "\n")
  cat("\nCoefficients:\n")
  print(stats::setNames(x$coef, x$coef_names))
  invisible(x)
}

vcov.robustcause_rlm <- function(object, hc_type = "HC3", ...) {
  .Call(
    rc_r_vcov_rlm,
    object$X_input,
    object$y_input,
    object$control,
    .hc_name(hc_type)
  )
}

confint.robustcause_rlm <- function(object,
                                    parm,
                                    level = 0.95,
                                    hc_type = "HC3",
                                    ...) {
  if (!missing(parm)) {
    warning("`parm` is ignored for robustcause fits; returning all coefficients.", call. = FALSE)
  }

  zcrit <- stats::qnorm((1 + level) / 2)
  ci <- .Call(
    rc_r_confint_rlm,
    object$X_input,
    object$y_input,
    object$control,
    .hc_name(hc_type),
    as.numeric(zcrit)
  )$ci
  rownames(ci) <- object$coef_names
  colnames(ci) <- c("lower", "upper")
  ci
}

.lm_hc_vcov <- function(object, hc_type = "HC3") {
  hc_type <- .hc_name(hc_type)
  X <- stats::model.matrix(object)
  e <- stats::residuals(object)
  n <- nrow(X)
  p <- ncol(X)
  h <- stats::hatvalues(object)
  one_minus_h <- pmax(1 - h, 1e-12)

  omega <- switch(
    hc_type,
    HC0 = e^2,
    HC1 = e^2 * (n / (n - p)),
    HC2 = e^2 / one_minus_h,
    HC3 = e^2 / (one_minus_h^2),
    HC4 = {
      delta <- pmin(4, n * h / p)
      e^2 / pmax((one_minus_h^delta), 1e-12)
    },
    HC4m = {
      nhp <- n * h / p
      delta <- pmin(1, nhp) + pmin(1.5, nhp)
      e^2 / pmax((one_minus_h^delta), 1e-12)
    },
    HC5 = {
      k <- 0.7
      cap <- max(4, n * k * max(h) / p)
      delta <- pmin(cap, n * h / p)
      e^2 / pmax((one_minus_h^(delta / 2)), 1e-12)
    }
  )

  bread <- solve(crossprod(X))
  meat <- crossprod(X, X * omega)
  bread %*% meat %*% bread
}

vcov_robust <- function(object, hc_type = "HC3") {
  if (inherits(object, "robustcause_rlm")) {
    return(vcov(object, hc_type = hc_type))
  }
  if (inherits(object, "lm")) {
    return(.lm_hc_vcov(object, hc_type = hc_type))
  }
  stop("`object` must inherit from `robustcause_rlm` or `lm`.", call. = FALSE)
}

confint_robust <- function(object, hc_type = "HC3", level = 0.95) {
  vc <- vcov_robust(object, hc_type = hc_type)
  zcrit <- stats::qnorm((1 + level) / 2)
  if (inherits(object, "robustcause_rlm")) {
    beta <- object$coef
    names(beta) <- object$coef_names
  } else {
    beta <- stats::coef(object)
  }
  se <- sqrt(pmax(diag(vc), 0))
  out <- cbind(lower = beta - zcrit * se, upper = beta + zcrit * se)
  rownames(out) <- names(beta)
  out
}

`%||%` <- function(lhs, rhs) {
  if (is.null(lhs)) rhs else lhs
}
