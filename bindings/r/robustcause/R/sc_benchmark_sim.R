simulate_sc_panel <- function(n_donors = 20L,
                              n_pre = 30L,
                              n_post = 10L,
                              n_factors = 2L,
                              true_sparsity = min(5L, n_donors),
                              factor_rho = 0.7,
                              noise_sd = 1,
                              treatment_effect = 1,
                              effect_shape = c("constant", "ramp", "delayed", "temporary", "none"),
                              contamination = c("none", "treated_pre_spike", "donor_pre_spike", "treated_and_donor_pre_spike", "heavy_tails", "donor_pool_mismatch", "predictor_spike", "mixed_outcome_predictor", "bad_donor_pre_path", "bad_donor_predictors", "bad_donor_mixed"),
                              contamination_rate = 0.05,
                              contamination_magnitude = 8,
                              donor_mismatch_strength = 0,
                              heavy_tail_df = 3,
                              n_predictors = 0L,
                              predictor_noise_sd = 0.25,
                              predictor_contamination_rate = contamination_rate,
                              predictor_contamination_magnitude = contamination_magnitude,
                              bad_donor_rate = 0.10,
                              seed = NULL) {
  effect_shape <- match.arg(effect_shape)
  contamination <- match.arg(contamination)
  if (!is.null(seed)) set.seed(seed)

  n_donors <- as.integer(n_donors)
  n_pre <- as.integer(n_pre)
  n_post <- as.integer(n_post)
  n_factors <- as.integer(n_factors)
  n_predictors <- as.integer(n_predictors)
  true_sparsity <- as.integer(min(true_sparsity, n_donors))
  if (n_donors < 2L) stop("`n_donors` must be at least 2.", call. = FALSE)
  if (n_pre < 2L || n_post < 1L) stop("Need at least two pre-periods and one post-period.", call. = FALSE)
  if (n_predictors < 0L) stop("`n_predictors` must be nonnegative.", call. = FALSE)

  n_time <- n_pre + n_post
  factors <- matrix(0, n_time, n_factors)
  factors[1, ] <- stats::rnorm(n_factors)
  for (tt in 2:n_time) {
    factors[tt, ] <- factor_rho * factors[tt - 1L, ] + stats::rnorm(n_factors)
  }

  loadings <- matrix(stats::rnorm(n_donors * n_factors), n_donors, n_factors)
  donors <- factors %*% t(loadings)
  donors <- sweep(donors, 2L, stats::rnorm(n_donors, sd = 0.5), "+")
  donors <- donors + matrix(stats::rnorm(n_time * n_donors, sd = noise_sd), n_time, n_donors)
  colnames(donors) <- paste0("donor_", seq_len(n_donors))

  active <- sort(sample(seq_len(n_donors), true_sparsity))
  raw_weights <- stats::rexp(true_sparsity)
  true_weights <- rep(0, n_donors)
  true_weights[active] <- raw_weights / sum(raw_weights)
  names(true_weights) <- colnames(donors)

  counterfactual <- as.numeric(donors %*% true_weights) + stats::rnorm(n_time, sd = noise_sd)
  if (identical(contamination, "donor_pool_mismatch") || donor_mismatch_strength > 0) {
    mismatch <- as.numeric(stats::filter(stats::rnorm(n_time), filter = factor_rho, method = "recursive"))
    mismatch[!is.finite(mismatch)] <- 0
    counterfactual <- counterfactual + donor_mismatch_strength * mismatch
  }

  post_effect <- switch(
    effect_shape,
    constant = rep(treatment_effect, n_post),
    ramp = seq(0, treatment_effect, length.out = n_post),
    delayed = c(rep(0, floor(n_post / 2)), rep(treatment_effect, n_post - floor(n_post / 2))),
    temporary = c(rep(treatment_effect, ceiling(n_post / 2)), rep(0, n_post - ceiling(n_post / 2))),
    none = rep(0, n_post)
  )
  treatment_path <- c(rep(0, n_pre), post_effect)
  outcomes <- cbind(treated = counterfactual + treatment_path, donors)
  outcomes_clean <- outcomes

  predictors <- simulate_sc_predictors(
    outcomes_pre_clean = outcomes_clean[seq_len(n_pre), , drop = FALSE],
    true_weights = true_weights,
    n_predictors = n_predictors,
    predictor_noise_sd = predictor_noise_sd
  )
  predictors_clean <- predictors

  add_spikes <- function(rows, cols) {
    cells <- expand.grid(row = rows, col = cols, KEEP.OUT.ATTRS = FALSE)
    n_pick <- min(nrow(cells), max(1L, floor(contamination_rate * nrow(cells))))
    picked <- cells[sample(seq_len(nrow(cells)), n_pick), , drop = FALSE]
    shock <- sample(c(-1, 1), n_pick, replace = TRUE) * contamination_magnitude
    for (ii in seq_len(n_pick)) outcomes[picked$row[ii], picked$col[ii]] <<- outcomes[picked$row[ii], picked$col[ii]] + shock[ii]
    picked$shock <- shock
    picked
  }

  add_predictor_spikes <- function() {
    if (is.null(predictors) || nrow(predictors) == 0L) {
      return(data.frame(row = integer(), col = integer(), shock = numeric()))
    }
    cells <- expand.grid(row = seq_len(nrow(predictors)), col = seq_len(ncol(predictors)), KEEP.OUT.ATTRS = FALSE)
    n_pick <- min(nrow(cells), max(1L, floor(predictor_contamination_rate * nrow(cells))))
    picked <- cells[sample(seq_len(nrow(cells)), n_pick), , drop = FALSE]
    shock <- sample(c(-1, 1), n_pick, replace = TRUE) * predictor_contamination_magnitude
    for (ii in seq_len(n_pick)) predictors[picked$row[ii], picked$col[ii]] <<- predictors[picked$row[ii], picked$col[ii]] + shock[ii]
    picked$shock <- shock
    picked
  }

  pick_bad_donors <- function() {
    n_bad <- min(n_donors, max(1L, ceiling(bad_donor_rate * n_donors)))
    sort(sample(seq_len(n_donors), n_bad))
  }

  contaminate_bad_donor_paths <- function(bad_donors) {
    bad_cols <- bad_donors + 1L
    pre_shift <- matrix(stats::rnorm(n_pre * length(bad_cols), sd = contamination_magnitude), nrow = n_pre)
    drift <- seq(0, contamination_magnitude, length.out = n_pre)
    for (jj in seq_along(bad_cols)) {
      outcomes[seq_len(n_pre), bad_cols[jj]] <<- outcomes[seq_len(n_pre), bad_cols[jj]] + pre_shift[, jj] + sample(c(-1, 1), 1) * drift
    }
    data.frame(row = rep(seq_len(n_pre), times = length(bad_cols)), col = rep(bad_cols, each = n_pre), shock = NA_real_)
  }

  contaminate_bad_donor_predictors <- function(bad_donors) {
    if (is.null(predictors) || nrow(predictors) == 0L) {
      return(data.frame(row = integer(), col = integer(), shock = numeric()))
    }
    bad_cols <- bad_donors + 1L
    rows <- seq_len(nrow(predictors))
    shock <- matrix(stats::rnorm(length(rows) * length(bad_cols), sd = predictor_contamination_magnitude), nrow = length(rows))
    for (jj in seq_along(bad_cols)) {
      predictors[rows, bad_cols[jj]] <<- predictors[rows, bad_cols[jj]] + shock[, jj]
    }
    data.frame(row = rep(rows, times = length(bad_cols)), col = rep(bad_cols, each = length(rows)), shock = as.numeric(shock))
  }

  contamination_info <- data.frame(row = integer(), col = integer(), shock = numeric())
  predictor_contamination_info <- data.frame(row = integer(), col = integer(), shock = numeric())
  bad_donor_info <- data.frame(donor = character(), donor_col = integer())
  if (identical(contamination, "treated_pre_spike")) contamination_info <- add_spikes(seq_len(n_pre), 1L)
  if (identical(contamination, "donor_pre_spike")) contamination_info <- add_spikes(seq_len(n_pre), seq.int(2L, n_donors + 1L))
  if (identical(contamination, "treated_and_donor_pre_spike")) contamination_info <- add_spikes(seq_len(n_pre), seq_len(n_donors + 1L))
  if (identical(contamination, "heavy_tails")) outcomes <- outcomes + matrix(stats::rt(n_time * (n_donors + 1L), df = heavy_tail_df), n_time, n_donors + 1L)
  if (identical(contamination, "predictor_spike")) predictor_contamination_info <- add_predictor_spikes()
  if (identical(contamination, "mixed_outcome_predictor")) {
    contamination_info <- add_spikes(seq_len(n_pre), seq_len(n_donors + 1L))
    predictor_contamination_info <- add_predictor_spikes()
  }
  if (contamination %in% c("bad_donor_pre_path", "bad_donor_predictors", "bad_donor_mixed")) {
    bad_donors <- pick_bad_donors()
    bad_donor_info <- data.frame(donor = colnames(donors)[bad_donors], donor_col = bad_donors + 1L)
    if (contamination %in% c("bad_donor_pre_path", "bad_donor_mixed")) {
      contamination_info <- contaminate_bad_donor_paths(bad_donors)
    }
    if (contamination %in% c("bad_donor_predictors", "bad_donor_mixed")) {
      predictor_contamination_info <- contaminate_bad_donor_predictors(bad_donors)
    }
  }

  panel <- data.frame(unit = rep(colnames(outcomes), each = n_time), time = rep(seq_len(n_time), times = ncol(outcomes)), outcome = as.numeric(outcomes))

  structure(list(
    panel = panel,
    outcomes = outcomes,
    outcomes_clean = outcomes_clean,
    predictors = predictors,
    predictors_clean = predictors_clean,
    contamination_info = contamination_info,
    predictor_contamination_info = predictor_contamination_info,
    bad_donor_info = bad_donor_info,
    truth = list(true_weights = true_weights, treatment_path = treatment_path, treatment_effect_post = post_effect, treatment_start = n_pre + 1L),
    design = list(n_donors = n_donors, n_pre = n_pre, n_post = n_post, n_factors = n_factors, true_sparsity = true_sparsity, noise_sd = noise_sd, treatment_effect = treatment_effect, effect_shape = effect_shape, contamination = contamination, contamination_rate = contamination_rate, contamination_magnitude = contamination_magnitude, donor_mismatch_strength = donor_mismatch_strength, n_predictors = n_predictors, predictor_noise_sd = predictor_noise_sd, predictor_contamination_rate = predictor_contamination_rate, predictor_contamination_magnitude = predictor_contamination_magnitude, bad_donor_rate = bad_donor_rate)
  ), class = "robustcause_sc_simulation")
}

simulate_sc_predictors <- function(outcomes_pre_clean,
                                   true_weights,
                                   n_predictors,
                                   predictor_noise_sd) {
  if (n_predictors <= 0L) return(NULL)

  unit_names <- colnames(outcomes_pre_clean)
  n_units <- ncol(outcomes_pre_clean)
  n_pre <- nrow(outcomes_pre_clean)
  donor_weights <- as.numeric(true_weights)

  base_summaries <- rbind(
    pre_mean = colMeans(outcomes_pre_clean),
    pre_sd = apply(outcomes_pre_clean, 2L, stats::sd),
    pre_first_half = colMeans(outcomes_pre_clean[seq_len(max(1L, floor(n_pre / 2L))), , drop = FALSE]),
    pre_second_half = colMeans(outcomes_pre_clean[seq.int(max(1L, floor(n_pre / 2L)), n_pre), , drop = FALSE])
  )

  predictors <- matrix(NA_real_, nrow = n_predictors, ncol = n_units)
  for (k in seq_len(n_predictors)) {
    coeff <- stats::rnorm(nrow(base_summaries))
    predictors[k, ] <- as.numeric(crossprod(coeff, base_summaries)) + stats::rnorm(n_units, sd = predictor_noise_sd)
  }

  if (n_units > 1L) {
    treated_signal <- as.numeric(predictors[, -1L, drop = FALSE] %*% donor_weights)
    predictors[, 1L] <- treated_signal + stats::rnorm(n_predictors, sd = predictor_noise_sd)
  }

  rownames(predictors) <- paste0("x", seq_len(n_predictors))
  colnames(predictors) <- unit_names
  predictors
}
