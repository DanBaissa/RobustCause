sc_null_coalesce <- function(x, y) {
  if (is.null(x)) y else x
}

sc_softmax <- function(z) {
  z <- as.numeric(z)
  z <- z - max(z)
  exp_z <- exp(z)
  exp_z / sum(exp_z)
}

sc_simplex_weights <- function(n_donors) {
  rep(1 / n_donors, n_donors)
}

sc_rmspe <- function(x) {
  sqrt(mean(as.numeric(x)^2))
}

sc_tukey_weights <- function(residuals, scale, c = 4.685, min_weight = 1e-8) {
  u <- as.numeric(residuals) / max(scale, .Machine$double.eps)
  w <- (1 - (u / c)^2)^2
  w[abs(u) >= c] <- 0
  pmax(w, min_weight)
}
