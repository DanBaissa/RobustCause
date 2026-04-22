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

.build_formula_rlm_frame <- function(formula,
                                     data,
                                     weights = NULL,
                                     weights_expr = NULL,
                                     subset = NULL,
                                     subset_expr = NULL,
                                     na.action = stats::na.omit,
                                     contrasts = NULL,
                                     offset = NULL,
                                     offset_expr = NULL) {
  if (is.null(data)) {
    stop("`data` is required when `x` is a formula.", call. = FALSE)
  }
  data_eval <- data
  if (!is.null(offset_expr) && !identical(offset_expr, quote(NULL)) && !identical(offset_expr, quote(expr = ))) {
    data_eval$.robustcause_offset <- eval(offset_expr, envir = data_eval, enclos = parent.frame())
    formula <- stats::update.formula(formula, . ~ . + offset(.robustcause_offset))
  }
  mf_call <- as.call(list(
    quote(stats::model.frame),
    formula = formula,
    data = data_eval,
    subset = subset_expr %||% subset,
    na.action = na.action,
    weights = weights_expr %||% weights
  ))
  mf <- eval(mf_call, envir = parent.frame())
  mt <- attr(mf, "terms")
  y <- model.response(mf)
  model_weights <- model.weights(mf)
  if (!length(model_weights)) {
    model_weights <- NULL
  }
  off <- model.offset(mf)
  if (!is.null(off)) {
    y <- y - off
  }
  X <- model.matrix(mt, mf, contrasts.arg = contrasts)
  list(
    X = X,
    y = as.numeric(y),
    model_weights = model_weights,
    formula = deparse(formula),
    mf = mf,
    terms = mt,
    offset = off,
    contrasts = attr(X, "contrasts"),
    xlevels = .getXlevels(mt, mf),
    na.action = attr(mf, "na.action")
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
  match.arg(tolower(method), c("m", "mm"))
}

.psi_name <- function(psi, method) {
  if (is.null(psi)) {
    return(if (identical(method, "mm")) "tukey_bisquare" else "huber")
  }
  if (is.function(psi)) {
    return("custom")
  }
  match.arg(psi, c("huber", "tukey_bisquare", "hampel"))
}

.scale_est_name <- function(scale.est) {
  if (is.null(scale.est)) {
    return("MAD")
  }
  match.arg(scale.est, c("MAD", "Huber", "proposal 2"))
}

.test_vec_name <- function(test.vec) {
  if (is.null(test.vec)) {
    return("coef")
  }
  match.arg(test.vec, c("coef", "resid", "w"))
}

.init_spec <- function(init, p) {
  if (is.null(init)) {
    return(list(method = "ls", coef = numeric()))
  }
  if (is.character(init) && length(init) == 1L) {
    return(list(method = match.arg(init, c("ls", "lts")), coef = numeric()))
  }
  if (is.numeric(init)) {
    init <- as.numeric(init)
    if (length(init) != p) {
      stop("Numeric `init` must have length equal to the number of coefficients.", call. = FALSE)
    }
    return(list(method = "user", coef = init))
  }
  if (is.list(init) && !is.null(init$coef)) {
    coef <- as.numeric(init$coef)
    if (length(coef) != p) {
      stop("`init$coef` must have length equal to the number of coefficients.", call. = FALSE)
    }
    return(list(method = "user", coef = coef))
  }
  stop("`init` must be NULL, 'ls', 'lts', a numeric coefficient vector, or a list with `$coef`.", call. = FALSE)
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
                              acc,
                              maxit,
                              tol,
                              add_intercept,
                              ridge,
                              min_weight,
                              prior_weights = numeric(),
                              wt.method = "case",
                              init = NULL,
                              scale.est = NULL,
                              k2 = NULL,
                              test.vec = NULL,
                              mm_s_control = NULL) {
  if (is.null(tuning)) {
    tuning <- if (identical(.rlm_method_name(method), "mm")) 4.685 else 1.345
  }
  if (is.function(psi)) {
    stop("Function-valued `psi` should be handled by the custom R wrapper path.", call. = FALSE)
  }
  psi_name <- .psi_name(psi, method)
  init_info <- .init_spec(init, if (length(prior_weights) > 0L) length(prior_weights) else 0L)
  list(
    method = .rlm_method_name(method),
    psi = psi_name,
    tuning = as.numeric(tuning),
    maxit = as.integer(maxit),
    tol = as.numeric(tol),
    acc = as.numeric(acc %||% tol),
    add_intercept = isTRUE(add_intercept),
    ridge = as.numeric(ridge),
    min_weight = as.numeric(min_weight),
    prior_weights = as.numeric(prior_weights),
    wt_method = match.arg(wt.method, c("inv.var", "case")),
    init_method = init_info$method,
    init_coef = as.numeric(init_info$coef),
    scale_est = .scale_est_name(scale.est),
    k2 = as.numeric(k2 %||% 1.345),
    test_vec = .test_vec_name(test.vec),
    mm_s_control = .coerce_mm_s_control(mm_s_control)
  )
}

.weighted_hatvalues <- function(X, w, ridge = 1e-10, min_weight = 1e-12) {
  active <- which(w > min_weight)
  h <- rep(0, nrow(X))
  if (length(active) == 0L) {
    return(h)
  }
  Xa <- X[active, , drop = FALSE]
  sw <- sqrt(w[active])
  Z <- Xa * sw
  XtWX <- crossprod(Z)
  diag(XtWX) <- diag(XtWX) + ridge
  XtWX_inv <- tryCatch(solve(XtWX), error = function(e) NULL)
  if (is.null(XtWX_inv)) {
    return(h)
  }
  for (j in seq_along(active)) {
    i <- active[[j]]
    zi <- X[i, , drop = FALSE] * sqrt(w[i])
    h[i] <- as.numeric(zi %*% XtWX_inv %*% t(zi))
  }
  h
}

.omega_hc <- function(e_working, h, n, p, hc_type) {
  one_minus_h <- pmax(1 - h, 1e-12)
  switch(
    .hc_name(hc_type),
    HC0 = e_working^2,
    HC1 = e_working^2 * (n / (n - p)),
    HC2 = e_working^2 / one_minus_h,
    HC3 = e_working^2 / (one_minus_h^2),
    HC4 = {
      delta <- pmin(4, n * h / p)
      e_working^2 / pmax(one_minus_h^delta, 1e-12)
    },
    HC4m = {
      nhp <- n * h / p
      delta <- pmin(1, nhp) + pmin(1.5, nhp)
      e_working^2 / pmax(one_minus_h^delta, 1e-12)
    },
    HC5 = {
      k <- 0.7
      cap <- max(4, n * k * max(h) / p)
      delta <- pmin(cap, n * h / p)
      e_working^2 / pmax(one_minus_h^(delta / 2), 1e-12)
    }
  )
}

.prepare_hc_storage <- function(X, resid, robust_w, wt.method = "case", prior_weights = NULL, ridge = 1e-10, min_weight = 1e-12) {
  wt.method <- match.arg(wt.method, c("inv.var", "case"))
  prior_weights <- if (is.null(prior_weights)) rep(1, nrow(X)) else as.numeric(prior_weights)
  if (wt.method == "inv.var") {
    fac <- sqrt(prior_weights)
    X_hc <- X * fac
    e_hc <- resid * fac
    w_hc <- robust_w
  } else {
    X_hc <- X
    e_hc <- resid
    w_hc <- robust_w * prior_weights
  }
  list(
    X = X_hc,
    resid = e_hc,
    w = w_hc,
    hat = .weighted_hatvalues(X_hc, w_hc, ridge = ridge, min_weight = min_weight)
  )
}

.fit_rlm_from_inputs <- function(inputs, control, prior_weights = NULL, wt.method = "case", metadata = list()) {
  control$prior_weights <- if (is.null(prior_weights)) numeric() else as.numeric(prior_weights)
  fit <- .Call(rc_r_fit_rlm, inputs$X, inputs$y, control)
  fit$coef_names <- .finalize_coef_names(inputs$coef_names, inputs$core_add_intercept, length(fit$coef))
  fit$formula <- inputs$formula
  fit$control <- control
  fit$X_input <- if (isTRUE(inputs$core_add_intercept)) cbind("(Intercept)" = 1, inputs$X) else inputs$X
  fit$y_input <- inputs$y
  fit$backend <- "cpp"
  hc_store <- .prepare_hc_storage(
    X = fit$X_input,
    resid = fit$resid,
    robust_w = fit$weights,
    wt.method = wt.method,
    prior_weights = prior_weights,
    ridge = control$ridge,
    min_weight = control$min_weight
  )
  fit$.hc_X <- hc_store$X
  fit$.hc_resid <- hc_store$resid
  fit$.hc_w <- hc_store$w
  fit$hat <- hc_store$hat
  fit$wt.method <- wt.method
  fit$prior_weights <- if (is.null(prior_weights)) NULL else as.numeric(prior_weights)
  fit$offset <- metadata$offset %||% NULL
  fit$terms <- metadata$terms %||% NULL
  fit$contrasts <- metadata$contrasts %||% NULL
  fit$xlevels <- metadata$xlevels %||% NULL
  fit$na.action <- metadata$na.action %||% NULL
  fit$model <- metadata$model %||% NULL
  class(fit) <- "robustcause_rlm"
  fit
}

.proposal2_gamma <- function(k) {
  phi <- dnorm(k)
  Phi <- pnorm(k)
  max(2 * (Phi - 0.5 - k * phi + k^2 * (1 - Phi)), 1e-8)
}

.robust_scale_r <- function(resid, scale.est = "MAD", k2 = 1.345) {
  scale.est <- .scale_est_name(scale.est)
  scale <- stats::mad(resid, constant = 1.4826)
  if (!is.finite(scale) || scale <= 0) {
    scale <- sqrt(max(sum(resid^2) / max(1, length(resid) - 1), 1e-16))
  }
  if (identical(scale.est, "MAD")) {
    return(scale)
  }
  gamma <- .proposal2_gamma(k2)
  for (i in seq_len(100L)) {
    cap <- k2 * scale
    updated <- sqrt(max(mean(pmin(resid^2, cap^2)) / gamma, 1e-16))
    if (!is.finite(updated) || updated <= 0) {
      break
    }
    if (abs(updated - scale) <= 1e-8 * max(1, scale)) {
      scale <- updated
      break
    }
    scale <- updated
  }
  scale
}

.solve_weighted_ls_r <- function(X, y, w, ridge = 1e-10, min_weight = 1e-12) {
  active <- which(w > min_weight)
  if (length(active) < ncol(X)) {
    stop("Too few active observations for weighted solve.", call. = FALSE)
  }
  Xa <- X[active, , drop = FALSE]
  ya <- y[active]
  sw <- sqrt(w[active])
  Xw <- Xa * sw
  yw <- ya * sw
  XtWX <- crossprod(Xw)
  diag(XtWX) <- diag(XtWX) + ridge
  XtWy <- crossprod(Xw, yw)
  as.numeric(solve(XtWX, XtWy))
}

.custom_psi_weights <- function(u, psi, psi_args = NULL, min_weight = 1e-12) {
  w <- do.call(psi, c(list(u), psi_args %||% list()))
  w <- as.numeric(w)
  if (length(w) != length(u)) {
    stop("Custom `psi` must return one weight per observation.", call. = FALSE)
  }
  w[!is.finite(w)] <- min_weight
  pmax(w, min_weight)
}

.fit_rlm_custom_psi <- function(inputs,
                                psi,
                                psi_args = NULL,
                                method,
                                tuning,
                                maxit,
                                tol,
                                acc,
                                ridge,
                                min_weight,
                                prior_weights = NULL,
                                wt.method = "case",
                                init = NULL,
                                scale.est = NULL,
                                k2 = NULL,
                                test.vec = NULL,
                                mm_s_control = NULL,
                                metadata = list()) {
  X <- if (isTRUE(inputs$core_add_intercept)) cbind("(Intercept)" = 1, inputs$X) else inputs$X
  y <- inputs$y
  prior_weights <- if (is.null(prior_weights)) rep(1, nrow(X)) else as.numeric(prior_weights)
  if (wt.method == "inv.var") {
    sw <- sqrt(prior_weights)
    X_work <- X * sw
    y_work <- y * sw
    base_weights <- rep(1, nrow(X))
  } else {
    X_work <- X
    y_work <- y
    base_weights <- prior_weights
  }

  init_info <- .init_spec(init, ncol(X))
  fixed_scale <- NA_real_
  if (identical(method, "mm")) {
    s_fit <- fit_s_estimator(
      X,
      y,
      add_intercept = FALSE,
      n_starts = mm_s_control$n_starts,
      n_best_starts = mm_s_control$n_best_starts,
      max_refine = mm_s_control$max_refine,
      max_scale_iter = mm_s_control$max_scale_iter,
      tol = mm_s_control$tol,
      scale_tol = mm_s_control$scale_tol,
      c = mm_s_control$c,
      b = mm_s_control$b,
      ridge = mm_s_control$ridge,
      min_weight = mm_s_control$min_weight,
      seed = mm_s_control$seed,
      use_fast_s = mm_s_control$use_fast_s,
      include_ols_start = mm_s_control$include_ols_start,
      fast_s_screen_subsets = mm_s_control$fast_s_screen_subsets,
      fast_s_screen_iters = mm_s_control$fast_s_screen_iters,
      fast_s_keep = mm_s_control$fast_s_keep
    )
    beta <- as.numeric(s_fit$coef)
    fixed_scale <- s_fit$scale
  } else if (identical(init_info$method, "user")) {
    beta <- init_info$coef
  } else if (identical(init_info$method, "lts")) {
    s_fit <- fit_s_estimator(
      X,
      y,
      add_intercept = FALSE,
      n_starts = mm_s_control$n_starts,
      n_best_starts = mm_s_control$n_best_starts,
      max_refine = mm_s_control$max_refine,
      max_scale_iter = mm_s_control$max_scale_iter,
      tol = mm_s_control$tol,
      scale_tol = mm_s_control$scale_tol,
      c = mm_s_control$c,
      b = mm_s_control$b,
      ridge = mm_s_control$ridge,
      min_weight = mm_s_control$min_weight,
      seed = mm_s_control$seed,
      use_fast_s = mm_s_control$use_fast_s,
      include_ols_start = mm_s_control$include_ols_start,
      fast_s_screen_subsets = mm_s_control$fast_s_screen_subsets,
      fast_s_screen_iters = mm_s_control$fast_s_screen_iters,
      fast_s_keep = mm_s_control$fast_s_keep
    )
    beta <- as.numeric(s_fit$coef)
  } else {
    beta <- .solve_weighted_ls_r(X_work, y_work, base_weights, ridge = ridge, min_weight = min_weight)
  }

  test.vec <- .test_vec_name(test.vec)
  k2 <- as.numeric(k2 %||% 1.345)
  scale.est <- .scale_est_name(scale.est)
  converged <- FALSE
  iterations <- 0L
  psi_w <- rep(1, nrow(X))

  for (iter in seq_len(maxit)) {
    resid <- y - as.numeric(X %*% beta)
    resid_work <- y_work - as.numeric(X_work %*% beta)
    scale <- if (is.finite(fixed_scale) && fixed_scale > 0) fixed_scale else .robust_scale_r(resid_work, scale.est, k2)
    psi_w <- .custom_psi_weights(resid_work / scale, psi, psi_args, min_weight = min_weight)
    combined_w <- pmax(base_weights * psi_w, min_weight)
    beta_new <- .solve_weighted_ls_r(X_work, y_work, combined_w, ridge = ridge, min_weight = min_weight)
    resid_work_new <- y_work - as.numeric(X_work %*% beta_new)
    scale_new <- if (is.finite(fixed_scale) && fixed_scale > 0) fixed_scale else .robust_scale_r(resid_work_new, scale.est, k2)
    psi_w_new <- .custom_psi_weights(resid_work_new / scale_new, psi, psi_args, min_weight = min_weight)
    diff <- switch(
      test.vec,
      resid = max(abs(resid_work_new - resid_work)) / max(1, max(abs(resid_work))),
      w = max(abs(psi_w_new - psi_w)) / max(1, max(abs(psi_w))),
      coef = max(abs(beta_new - beta)) / max(1, max(abs(beta)))
    )
    beta <- beta_new
    iterations <- iter
    if (diff < (acc %||% tol)) {
      converged <- TRUE
      break
    }
  }

  fitted <- as.numeric(X %*% beta)
  resid <- y - fitted
  resid_hc <- if (wt.method == "inv.var") resid * sqrt(prior_weights) else resid
  hc_weights <- if (wt.method == "inv.var") psi_w else psi_w * prior_weights
  hc_store <- .prepare_hc_storage(
    X = X,
    resid = resid,
    robust_w = psi_w,
    wt.method = wt.method,
    prior_weights = prior_weights,
    ridge = ridge,
    min_weight = min_weight
  )

  fit <- list(
    coef = as.numeric(beta),
    fitted = fitted,
    resid = resid,
    weights = psi_w,
    hat = hc_store$hat,
    scale = if (is.finite(fixed_scale) && fixed_scale > 0) fixed_scale else .robust_scale_r(if (wt.method == "inv.var") resid_hc else resid, scale.est, k2),
    converged = converged,
    iterations = iterations,
    method = toupper(method),
    psi = "custom",
    coef_names = .finalize_coef_names(inputs$coef_names, inputs$core_add_intercept, length(beta)),
    formula = inputs$formula,
    control = list(
      method = method,
      psi = "custom",
      tuning = tuning,
      maxit = maxit,
      tol = tol,
      acc = acc %||% tol,
      ridge = ridge,
      min_weight = min_weight,
      wt.method = wt.method,
      scale.est = scale.est,
      k2 = k2,
      test.vec = test.vec,
      init = init
    ),
    X_input = X,
    y_input = y,
    backend = "r_custom_psi",
    prior_weights = prior_weights,
    wt.method = wt.method,
    offset = metadata$offset %||% NULL,
    terms = metadata$terms %||% NULL,
    contrasts = metadata$contrasts %||% NULL,
    xlevels = metadata$xlevels %||% NULL,
    na.action = metadata$na.action %||% NULL,
    model = metadata$model %||% NULL,
    .hc_X = hc_store$X,
    .hc_resid = hc_store$resid,
    .hc_w = hc_store$w
  )
  class(fit) <- "robustcause_rlm"
  fit
}

.needs_mass_backend <- function(method,
                                psi,
                                weights,
                                init,
                                scale.est,
                                k2,
                                test.vec,
                                subset,
                                na.action,
                                contrasts,
                                offset) {
  psi_name <- .psi_name(psi, method)
  is.function(psi) ||
    identical(psi_name, "hampel") ||
    !is.null(weights) ||
    !is.null(init) ||
    !is.null(scale.est) ||
    !is.null(k2) ||
    !is.null(test.vec) ||
    !is.null(subset) ||
    !identical(na.action, stats::na.omit) ||
    !is.null(contrasts) ||
    !is.null(offset)
}

.mass_psi <- function(psi, method, tuning = NULL, psi_args = NULL) {
  if (is.function(psi)) {
    return(list(psi = psi, psi_args = psi_args %||% list()))
  }
  psi_name <- .psi_name(psi, method)
  out <- switch(
    psi_name,
    huber = list(psi = MASS::psi.huber, psi_args = c(list(k = tuning %||% 1.345), psi_args %||% list())),
    tukey_bisquare = list(psi = MASS::psi.bisquare, psi_args = c(list(c = tuning %||% 4.685), psi_args %||% list())),
    hampel = {
      cval <- tuning %||% 8
      args <- c(list(a = cval / 4, b = cval / 2, c = cval), psi_args %||% list())
      list(psi = MASS::psi.hampel, psi_args = args)
    },
    stop("Unsupported psi for MASS backend.", call. = FALSE)
  )
  out
}

.fit_rlm_mass <- function(x,
                          y = NULL,
                          data = NULL,
                          method,
                          psi,
                          tuning,
                          psi_args = NULL,
                          maxit,
                          tol,
                          acc,
                          add_intercept,
                          ridge,
                          min_weight,
                          weights = NULL,
                          wt.method = c("inv.var", "case"),
                          init = NULL,
                          scale.est = NULL,
                          k2 = NULL,
                          test.vec = NULL,
                          subset = NULL,
                          subset_expr = NULL,
                          na.action = stats::na.omit,
                          contrasts = NULL,
                          offset = NULL,
                          offset_expr = NULL,
                          weights_expr = NULL,
                          mm_s_control = NULL,
                          call = NULL) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("The `MASS` package is required for these `fit_rlm()` options.", call. = FALSE)
  }
  mass_rlm_default <- getS3method("rlm", "default", envir = asNamespace("MASS"))
  wt.method <- match.arg(wt.method)
  if (identical(method, "mm") && !is.null(mm_s_control)) {
    warning("`mm_s_control` is ignored on the MASS-backed path.", call. = FALSE)
  }
  psi_spec <- .mass_psi(psi, method, tuning = tuning, psi_args = psi_args)
  psi_label <- if (is.function(psi)) {
    fn_name <- deparse(substitute(psi))
    if (identical(fn_name, "psi")) "custom" else fn_name
  } else {
    .psi_name(psi, method)
  }
  mass_method <- toupper(method)
  fit_call <- list(
    method = mass_method,
    psi = psi_spec$psi,
    maxit = as.integer(maxit),
    acc = as.numeric(acc %||% tol),
    wt.method = wt.method
  )
  if (!is.null(scale.est) && identical(method, "m")) {
    fit_call$scale.est <- scale.est
  }
  if (!is.null(k2) && identical(method, "m")) {
    fit_call$k2 <- as.numeric(k2)
  }
  if (!is.null(test.vec)) {
    fit_call$test.vec <- if (identical(test.vec, "NULL")) NULL else test.vec
  }
  if (!is.null(init)) {
    fit_call$init <- init
  }
  fit_call <- c(fit_call, psi_spec$psi_args)

  if (inherits(x, "formula")) {
    frame <- .build_formula_rlm_frame(
      formula = x,
      data = data,
      weights = weights,
      weights_expr = weights_expr,
      subset = subset,
      subset_expr = subset_expr,
      na.action = na.action,
      contrasts = contrasts,
      offset = offset,
      offset_expr = offset_expr
    )
    fit_args <- list(x = frame$X, y = frame$y)
    if (!is.null(frame$model_weights)) {
      fit_args$weights <- frame$model_weights
    }
    fit <- do.call(mass_rlm_default, c(fit_args, fit_call))
    coef_names <- colnames(frame$X)
    formula_label <- deparse(x)
    prior_weights <- frame$model_weights
    X_for_store <- frame$X
    y_for_store <- frame$y
    offset_used <- frame$offset
    terms_used <- frame$terms
    contrasts_used <- frame$contrasts
    xlevels_used <- frame$xlevels
    na_action_used <- frame$na.action
    model_used <- frame$mf
  } else {
    inputs <- .as_matrix_inputs(x, y, data = NULL, add_intercept = add_intercept)
    prior_weights <- weights
    fit_args <- list(x = inputs$X, y = inputs$y)
    if (!is.null(weights)) {
      fit_args$weights <- weights
    }
    fit <- do.call(mass_rlm_default, c(fit_args, fit_call))
    coef_names <- .finalize_coef_names(inputs$coef_names, inputs$core_add_intercept, length(fit$coefficients))
    formula_label <- NULL
    X_for_store <- if (isTRUE(add_intercept)) cbind("(Intercept)" = 1, inputs$X) else inputs$X
    y_for_store <- inputs$y
    offset_used <- NULL
    terms_used <- NULL
    contrasts_used <- NULL
    xlevels_used <- NULL
    na_action_used <- NULL
    model_used <- NULL
  }

  coef_vec <- as.numeric(fit$coefficients)
  names(coef_vec) <- coef_names
  robust_w <- fit$w
  if (is.null(robust_w)) {
    robust_w <- rep(1, nrow(X_for_store))
  }
  hc_store <- .prepare_hc_storage(
    X = X_for_store,
    resid = fit$residuals,
    robust_w = robust_w,
    wt.method = wt.method,
    prior_weights = prior_weights,
    ridge = ridge,
    min_weight = min_weight
  )

  out <- list(
    coef = coef_vec,
    fitted = fit$fitted.values,
    resid = fit$residuals,
    weights = robust_w,
    hat = hc_store$hat,
    scale = fit$s,
    converged = isTRUE(fit$converged),
    iterations = length(fit$conv),
    method = toupper(method),
    psi = psi_label,
    coef_names = coef_names,
    formula = formula_label,
    control = list(
      method = method,
      psi = .psi_name(psi, method),
      tuning = tuning,
      maxit = maxit,
      tol = tol,
      acc = acc %||% tol,
      ridge = ridge,
      min_weight = min_weight,
      wt.method = wt.method,
      scale.est = scale.est,
      k2 = k2,
      test.vec = test.vec,
      init = init
    ),
    X_input = X_for_store,
    y_input = y_for_store,
    backend = "mass",
    mass_fit = fit,
    prior_weights = prior_weights,
    wt.method = wt.method,
    offset = offset_used,
    terms = terms_used,
    contrasts = contrasts_used,
    xlevels = xlevels_used,
    na.action = na_action_used,
    model = model_used,
    .hc_X = hc_store$X,
    .hc_resid = hc_store$resid,
    .hc_w = hc_store$w
  )
  class(out) <- "robustcause_rlm"
  out$call <- call %||% match.call()
  out
}

fit_rlm <- function(x,
                    y = NULL,
                    data = NULL,
                    method = c("m", "mm"),
                    psi = NULL,
                    tuning = NULL,
                    psi_args = NULL,
                    maxit = 100L,
                    tol = 1e-8,
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
                    mm_s_control = NULL) {
  method <- .rlm_method_name(method)
  wt.method <- match.arg(wt.method)
  mm_s_control <- .coerce_mm_s_control(mm_s_control)

  if (inherits(x, "formula")) {
    frame <- .build_formula_rlm_frame(
      formula = x,
      data = data,
      weights = weights,
      weights_expr = substitute(weights),
      subset = subset,
      subset_expr = substitute(subset),
      na.action = na.action,
      contrasts = contrasts,
      offset = offset,
      offset_expr = substitute(offset)
    )
    inputs <- list(
      X = unname(as.matrix(frame$X)),
      y = as.numeric(frame$y),
      coef_names = colnames(frame$X),
      formula = deparse(x),
      core_add_intercept = FALSE
    )
    if (is.function(psi)) {
      fit <- .fit_rlm_custom_psi(
        inputs = inputs,
        psi = psi,
        psi_args = psi_args,
        method = method,
        tuning = tuning %||% if (identical(method, "mm")) 4.685 else 1.345,
        maxit = maxit,
        tol = tol,
        acc = acc,
        ridge = ridge,
        min_weight = min_weight,
        prior_weights = frame$model_weights,
        wt.method = wt.method,
        init = init,
        scale.est = scale.est,
        k2 = k2,
        test.vec = test.vec,
        mm_s_control = mm_s_control,
        metadata = list(
          offset = frame$offset,
          terms = frame$terms,
          contrasts = frame$contrasts,
          xlevels = frame$xlevels,
          na.action = frame$na.action,
          model = frame$mf
        )
      )
      fit$call <- match.call()
      return(fit)
    }
    init_info <- .init_spec(init, ncol(frame$X))
    control <- list(
      method = method,
      psi = .psi_name(psi, method),
      tuning = as.numeric(tuning %||% if (identical(method, "mm")) 4.685 else 1.345),
      maxit = as.integer(maxit),
      tol = as.numeric(tol),
      acc = as.numeric(acc %||% tol),
      add_intercept = FALSE,
      ridge = as.numeric(ridge),
      min_weight = as.numeric(min_weight),
      prior_weights = if (is.null(frame$model_weights)) numeric() else as.numeric(frame$model_weights),
      wt_method = wt.method,
      init_method = init_info$method,
      init_coef = as.numeric(init_info$coef),
      scale_est = .scale_est_name(scale.est),
      k2 = as.numeric(k2 %||% 1.345),
      test_vec = .test_vec_name(test.vec),
      mm_s_control = mm_s_control
    )
    fit <- .fit_rlm_from_inputs(
      inputs,
      control,
      prior_weights = frame$model_weights,
      wt.method = wt.method,
      metadata = list(
        offset = frame$offset,
        terms = frame$terms,
        contrasts = frame$contrasts,
        xlevels = frame$xlevels,
        na.action = frame$na.action,
        model = frame$mf
      )
    )
  } else {
    inputs <- .as_matrix_inputs(x, y, data, add_intercept = add_intercept)
    if (is.function(psi)) {
      fit <- .fit_rlm_custom_psi(
        inputs = inputs,
        psi = psi,
        psi_args = psi_args,
        method = method,
        tuning = tuning %||% if (identical(method, "mm")) 4.685 else 1.345,
        maxit = maxit,
        tol = tol,
        acc = acc,
        ridge = ridge,
        min_weight = min_weight,
        prior_weights = weights,
        wt.method = wt.method,
        init = init,
        scale.est = scale.est,
        k2 = k2,
        test.vec = test.vec,
        mm_s_control = mm_s_control
      )
      fit$call <- match.call()
      return(fit)
    }
    p <- ncol(if (isTRUE(inputs$core_add_intercept)) cbind(1, inputs$X) else inputs$X)
    init_info <- .init_spec(init, p)
    control <- list(
      method = method,
      psi = .psi_name(psi, method),
      tuning = as.numeric(tuning %||% if (identical(method, "mm")) 4.685 else 1.345),
      maxit = as.integer(maxit),
      tol = as.numeric(tol),
      acc = as.numeric(acc %||% tol),
      add_intercept = isTRUE(inputs$core_add_intercept),
      ridge = as.numeric(ridge),
      min_weight = as.numeric(min_weight),
      prior_weights = if (is.null(weights)) numeric() else as.numeric(weights),
      wt_method = wt.method,
      init_method = init_info$method,
      init_coef = as.numeric(init_info$coef),
      scale_est = .scale_est_name(scale.est),
      k2 = as.numeric(k2 %||% 1.345),
      test_vec = .test_vec_name(test.vec),
      mm_s_control = mm_s_control
    )
    fit <- .fit_rlm_from_inputs(inputs, control, prior_weights = weights, wt.method = wt.method)
  }

  fit$call <- match.call()
  fit
}

print.robustcause_rlm <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

summary.robustcause_rlm <- function(object, hc_type = "HC3", level = 0.95, ...) {
  vc <- vcov(object, hc_type = hc_type)
  se <- sqrt(pmax(diag(vc), 0))
  zcrit <- stats::qnorm((1 + level) / 2)
  coef_table <- cbind(
    Estimate = object$coef,
    `Std. Error` = se,
    `z value` = object$coef / se,
    `Pr(>|z|)` = 2 * stats::pnorm(abs(object$coef / se), lower.tail = FALSE),
    lower = object$coef - zcrit * se,
    upper = object$coef + zcrit * se
  )
  rownames(coef_table) <- object$coef_names
  colnames(coef_table) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)", "lower", "upper")

  structure(
    list(
      call = object$call,
      formula = object$formula,
      method = object$method,
      psi = object$psi,
      converged = object$converged,
      iterations = object$iterations,
      scale = object$scale,
      hc_type = hc_type,
      level = level,
      coefficients = coef_table
    ),
    class = "summary.robustcause_rlm"
  )
}

print.summary.robustcause_rlm <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("RobustCause RLM fit\n")
  if (!is.null(x$formula)) {
    cat("Formula:", x$formula, "\n")
  }
  cat("Method:", x$method, "\n")
  cat("Psi:", x$psi, "\n")
  cat("Converged:", x$converged, "\n")
  cat("Iterations:", x$iterations, "\n")
  cat("Scale:", format(signif(x$scale, 6)), "\n")
  cat("Inference:", x$hc_type, "\n")
  cat("\nCoefficients:\n")
  stats::printCoefmat(x$coefficients[, c("Estimate", "Std. Error", "z value", "Pr(>|z|)")], digits = digits, signif.stars = TRUE, na.print = "NA")
  invisible(x)
}

.rlm_hc_vcov <- function(object, hc_type = "HC3") {
  X <- object$.hc_X %||% object$X_input
  e <- object$.hc_resid %||% object$resid
  w <- object$.hc_w %||% object$weights
  h <- object$hat %||% .weighted_hatvalues(X, w, ridge = object$control$ridge %||% 1e-10, min_weight = object$control$min_weight %||% 1e-12)
  n <- nrow(X)
  p <- ncol(X)
  e_working <- w * e
  omega <- .omega_hc(e_working, h, n, p, hc_type)
  XtWX <- crossprod(X, sweep(X, 1L, w, `*`))
  diag(XtWX) <- diag(XtWX) + (object$control$ridge %||% 1e-10)
  bread <- solve(XtWX)
  meat <- crossprod(X, sweep(X, 1L, omega, `*`))
  bread %*% meat %*% bread
}

vcov.robustcause_rlm <- function(object, hc_type = "HC3", ...) {
  .rlm_hc_vcov(object, hc_type = hc_type)
}

confint.robustcause_rlm <- function(object,
                                    parm,
                                    level = 0.95,
                                    hc_type = "HC3",
                                    ...) {
  if (!missing(parm)) {
    warning("`parm` is ignored for robustcause fits; returning all coefficients.", call. = FALSE)
  }
  vc <- vcov(object, hc_type = hc_type)
  zcrit <- stats::qnorm((1 + level) / 2)
  se <- sqrt(pmax(diag(vc), 0))
  ci <- cbind(lower = object$coef - zcrit * se, upper = object$coef + zcrit * se)
  rownames(ci) <- object$coef_names
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
  meat <- crossprod(X, sweep(X, 1L, omega, `*`))
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
