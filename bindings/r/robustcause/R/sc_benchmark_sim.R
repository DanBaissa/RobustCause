simulate_sc_panel <- function(n_donors = 20L,
                              n_pre = 30L,
                              n_post = 10L,
                              n_factors = 2L,
                              true_sparsity = min(5L, n_donors),
                              factor_rho = 0.7,
                              noise_sd = 1,
                              treatment_effect = 1,
                              effect_shape = c("constant", "ramp", "delayed", "temporary", "none"),
                              contamination = c("none", "treated_pre_spike", "donor_pre_spike", "treated_and_donor_pre_spike", "heavy_tails", "donor_pool_mismatch"),
                              contamination_rate = 0.05,
                              contamination_magnitude = 8,
                              donor_mismatch_strength = 0,
                              heavy_tail_df = 3,
                              seed = NULL) {
  effect_shape <- match.arg(effect_shape)
  contamination <- match.arg(contamination)
  if (!is.null(seed)) set.seed(seed)

  n_donors <- as.integer(n_donors)
  n_pre <- as.integer(n_pre)
  n_post <- as.integer(n_post)
  n_factors <- as.integer(n_factors)
  true_sparsity <- as.integer(min(true_sparsity, n_donors))
  if (n_donors < 2L) stop("`n_donors` must be at least 2.", call. = FALSE)
  if (n_pre < 2L || n_post < 1L) stop("Need at least two pre-periods and one post-period.", call. = FALSE)

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

  add_spikes <- function(rows, cols) {
    cells <- expand.grid(row = rows, col = cols, KEEP.OUT.ATTRS = FALSE)
    n_pick <- min(nrow(cells), max(1L, floor(contamination_rate * nrow(cells))))
    picked <- cells[sample(seq_len(nrow(cells)), n_pick), , drop = FALSE]
    shock <- sample(c(-1, 1), n_pick, replace = TRUE) * contamination_magnitude
    for (ii in seq_len(n_pick)) outcomes[picked$row[ii], picked$col[ii]] <<- outcomes[picked$row[ii], picked$col[ii]] + shock[ii]
    picked$shock <- shock
    picked
  }

  contamination_info <- data.frame(row = integer(), col = integer(), shock = numeric())
  if (identical(contamination, "treated_pre_spike")) contamination_info <- add_spikes(seq_len(n_pre), 1L)
  if (identical(contamination, "donor_pre_spike")) contamination_info <- add_spikes(seq_len(n_pre), seq.int(2L, n_donors + 1L))
  if (identical(contamination, "treated_and_donor_pre_spike")) contamination_info <- add_spikes(seq_len(n_pre), seq_len(n_donors + 1L))
  if (identical(contamination, "heavy_tails")) outcomes <- outcomes + matrix(stats::rt(n_time * (n_donors + 1L), df = heavy_tail_df), n_time, n_donors + 1L)

  panel <- data.frame(unit = rep(colnames(outcomes), each = n_time), time = rep(seq_len(n_time), times = ncol(outcomes)), outcome = as.numeric(outcomes))

  structure(list(
    panel = panel,
    outcomes = outcomes,
    outcomes_clean = outcomes_clean,
    contamination_info = contamination_info,
    truth = list(true_weights = true_weights, treatment_path = treatment_path, treatment_effect_post = post_effect, treatment_start = n_pre + 1L),
    design = list(n_donors = n_donors, n_pre = n_pre, n_post = n_post, n_factors = n_factors, true_sparsity = true_sparsity, noise_sd = noise_sd, treatment_effect = treatment_effect, effect_shape = effect_shape, contamination = contamination, contamination_rate = contamination_rate, contamination_magnitude = contamination_magnitude, donor_mismatch_strength = donor_mismatch_strength)
  ), class = "robustcause_sc_simulation")
}
