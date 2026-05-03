sc_prepare_matrix <- function(outcomes,
                              treated_unit = NULL,
                              treatment_start,
                              donors = NULL,
                              predictors = NULL,
                              predictor_weights = NULL,
                              predictor_lambda = 1,
                              predictor_scale = c("mad", "sd", "none"),
                              na.action = c("omit", "fail")) {
  predictor_scale <- match.arg(predictor_scale)
  na.action <- match.arg(na.action)
  outcomes <- as.matrix(outcomes)

  if (!is.numeric(outcomes)) {
    stop("`outcomes` must be numeric or coercible to a numeric matrix.", call. = FALSE)
  }
  if (identical(na.action, "omit")) {
    outcomes <- outcomes[stats::complete.cases(outcomes), , drop = FALSE]
  } else if (any(!is.finite(outcomes))) {
    stop("`outcomes` contains missing or non-finite values.", call. = FALSE)
  }
  if (ncol(outcomes) < 2L) {
    stop("`outcomes` must contain one treated unit and at least one donor.", call. = FALSE)
  }

  unit_names <- colnames(outcomes)
  if (is.null(unit_names)) {
    unit_names <- as.character(seq_len(ncol(outcomes)))
    colnames(outcomes) <- unit_names
  }

  treated_col <- if (is.null(treated_unit)) {
    1L
  } else if (is.numeric(treated_unit)) {
    as.integer(treated_unit[[1]])
  } else {
    match(as.character(treated_unit[[1]]), unit_names)
  }
  if (is.na(treated_col) || treated_col < 1L || treated_col > ncol(outcomes)) {
    stop("`treated_unit` could not be resolved.", call. = FALSE)
  }

  donor_cols <- if (is.null(donors)) {
    setdiff(seq_len(ncol(outcomes)), treated_col)
  } else if (is.numeric(donors)) {
    as.integer(donors)
  } else {
    match(as.character(donors), unit_names)
  }
  donor_cols <- unique(setdiff(donor_cols, treated_col))
  if (length(donor_cols) < 1L || anyNA(donor_cols)) {
    stop("`donors` could not be resolved.", call. = FALSE)
  }

  treatment_start <- as.integer(treatment_start[[1]])
  if (is.na(treatment_start) || treatment_start <= 1L || treatment_start > nrow(outcomes)) {
    stop("`treatment_start` must leave pre- and post-treatment periods.", call. = FALSE)
  }

  pre_idx <- seq_len(treatment_start - 1L)
  post_idx <- seq.int(treatment_start, nrow(outcomes))
  ordered_cols <- c(treated_col, donor_cols)

  pred <- sc_prepare_predictors(
    predictors = predictors,
    unit_names = unit_names,
    ordered_cols = ordered_cols,
    predictor_weights = predictor_weights,
    predictor_lambda = predictor_lambda,
    predictor_scale = predictor_scale
  )

  list(
    outcomes = outcomes[, ordered_cols, drop = FALSE],
    treated_pre = as.numeric(outcomes[pre_idx, treated_col]),
    donors_pre = unname(outcomes[pre_idx, donor_cols, drop = FALSE]),
    treated_post = as.numeric(outcomes[post_idx, treated_col]),
    donors_post = unname(outcomes[post_idx, donor_cols, drop = FALSE]),
    predictors = pred$predictors,
    predictor_treated = pred$predictor_treated,
    predictor_donors = pred$predictor_donors,
    predictor_weights = pred$predictor_weights,
    predictor_lambda = pred$predictor_lambda,
    predictor_scale = pred$predictor_scale,
    predictor_names = pred$predictor_names,
    treated_unit = unit_names[[treated_col]],
    donors = unit_names[donor_cols],
    pre_periods = pre_idx,
    post_periods = post_idx,
    treatment_start = treatment_start
  )
}

sc_prepare_predictors <- function(predictors,
                                  unit_names,
                                  ordered_cols,
                                  predictor_weights,
                                  predictor_lambda,
                                  predictor_scale) {
  n_ordered <- length(ordered_cols)

  if (is.null(predictors)) {
    return(list(
      predictors = matrix(numeric(), nrow = 0L, ncol = n_ordered),
      predictor_treated = numeric(),
      predictor_donors = matrix(numeric(), nrow = 0L, ncol = n_ordered - 1L),
      predictor_weights = numeric(),
      predictor_lambda = as.numeric(predictor_lambda),
      predictor_scale = predictor_scale,
      predictor_names = character()
    ))
  }

  predictors <- as.matrix(predictors)
  if (!is.numeric(predictors)) {
    stop("`predictors` must be numeric or coercible to a numeric matrix.", call. = FALSE)
  }
  if (ncol(predictors) != length(unit_names)) {
    stop("`predictors` must have one column per unit in `outcomes`.", call. = FALSE)
  }
  if (any(!is.finite(predictors))) {
    stop("`predictors` contains missing or non-finite values.", call. = FALSE)
  }

  predictor_names <- rownames(predictors)
  if (is.null(predictor_names)) {
    predictor_names <- paste0("x", seq_len(nrow(predictors)))
    rownames(predictors) <- predictor_names
  }

  predictors <- predictors[, ordered_cols, drop = FALSE]
  predictors <- sc_scale_predictors(predictors, predictor_scale)

  weights <- sc_predictor_weights(predictor_weights, nrow(predictors))
  predictor_lambda <- as.numeric(predictor_lambda[[1]])
  if (!is.finite(predictor_lambda) || predictor_lambda < 0) {
    stop("`predictor_lambda` must be a nonnegative finite number.", call. = FALSE)
  }

  list(
    predictors = predictors,
    predictor_treated = as.numeric(predictors[, 1L]),
    predictor_donors = unname(predictors[, -1L, drop = FALSE]),
    predictor_weights = weights,
    predictor_lambda = predictor_lambda,
    predictor_scale = predictor_scale,
    predictor_names = predictor_names
  )
}

sc_scale_predictors <- function(predictors, predictor_scale) {
  if (nrow(predictors) == 0L || identical(predictor_scale, "none")) {
    return(predictors)
  }

  center <- apply(predictors, 1L, stats::median)
  spread <- switch(
    predictor_scale,
    mad = apply(predictors, 1L, stats::mad, constant = 1.4826),
    sd = apply(predictors, 1L, stats::sd),
    rep(1, nrow(predictors))
  )
  spread[!is.finite(spread) | spread <= .Machine$double.eps] <- 1

  sweep(sweep(predictors, 1L, center, "-"), 1L, spread, "/")
}

sc_predictor_weights <- function(predictor_weights, n_predictors) {
  if (n_predictors == 0L) return(numeric())

  if (is.null(predictor_weights)) {
    return(rep(1, n_predictors))
  }

  weights <- as.numeric(predictor_weights)
  if (length(weights) != n_predictors || any(!is.finite(weights)) || any(weights < 0)) {
    stop("`predictor_weights` must be nonnegative and match the number of predictor rows.", call. = FALSE)
  }
  weights
}
