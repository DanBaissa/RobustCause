sc_prepare_matrix <- function(outcomes,
                              treated_unit = NULL,
                              treatment_start,
                              donors = NULL,
                              na.action = c("omit", "fail")) {
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

  list(
    outcomes = outcomes[, c(treated_col, donor_cols), drop = FALSE],
    treated_pre = as.numeric(outcomes[pre_idx, treated_col]),
    donors_pre = unname(outcomes[pre_idx, donor_cols, drop = FALSE]),
    treated_post = as.numeric(outcomes[post_idx, treated_col]),
    donors_post = unname(outcomes[post_idx, donor_cols, drop = FALSE]),
    treated_unit = unit_names[[treated_col]],
    donors = unit_names[donor_cols],
    pre_periods = pre_idx,
    post_periods = post_idx,
    treatment_start = treatment_start
  )
}
