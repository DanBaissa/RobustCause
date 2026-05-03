sc_compute_donor_diagnostics <- function(y,
                                         x,
                                         donor_names,
                                         donor_tukey_c = 4.685,
                                         min_donor_weight = 1e-4,
                                         donor_score = c("balance_distance")) {
  donor_score <- match.arg(donor_score)
  y <- as.numeric(y)
  x <- as.matrix(x)

  if (length(y) != nrow(x)) {
    stop("Internal error: donor diagnostic outcome and matrix rows do not match.", call. = FALSE)
  }
  if (ncol(x) != length(donor_names)) {
    stop("Internal error: donor diagnostic names do not match donor columns.", call. = FALSE)
  }

  if (ncol(x) == 0L) {
    return(sc_empty_donor_diagnostics())
  }

  row_values <- cbind(y, x)
  row_center <- apply(row_values, 1L, stats::median)
  row_scale <- apply(row_values, 1L, stats::mad, constant = 1.4826)
  row_sd <- apply(row_values, 1L, stats::sd)
  row_scale[!is.finite(row_scale) | row_scale <= .Machine$double.eps] <- row_sd[!is.finite(row_scale) | row_scale <= .Machine$double.eps]
  row_scale[!is.finite(row_scale) | row_scale <= .Machine$double.eps] <- 1

  y_scaled <- (y - row_center) / row_scale
  x_scaled <- sweep(sweep(x, 1L, row_center, "-"), 1L, row_scale, "/")
  donor_residuals <- sweep(x_scaled, 1L, y_scaled, "-")
  scores <- sqrt(colMeans(donor_residuals^2))
  scores[!is.finite(scores)] <- max(scores[is.finite(scores)], 1)

  score_center <- stats::median(scores)
  score_scale <- stats::mad(scores, constant = 1.4826)
  if (!is.finite(score_scale) || score_scale <= .Machine$double.eps) {
    score_scale <- stats::sd(scores)
  }
  if (!is.finite(score_scale) || score_scale <= .Machine$double.eps) {
    score_scale <- 1
  }

  excess <- pmax(0, scores - score_center)
  u <- excess / max(donor_tukey_c * score_scale, .Machine$double.eps)
  reliability <- ifelse(u < 1, (1 - u^2)^2, min_donor_weight)
  reliability <- pmax(min_donor_weight, pmin(1, reliability))
  reliability[!is.finite(reliability)] <- min_donor_weight

  names(scores) <- donor_names
  names(reliability) <- donor_names

  data.frame(
    donor = donor_names,
    donor_score = as.numeric(scores),
    robust_donor_weight = as.numeric(reliability),
    stringsAsFactors = FALSE
  )
}

sc_empty_donor_diagnostics <- function() {
  data.frame(
    donor = character(),
    donor_score = numeric(),
    robust_donor_weight = numeric(),
    stringsAsFactors = FALSE
  )
}

sc_donor_penalties <- function(donor_diagnostics,
                               robust_donors = FALSE,
                               donor_penalty_lambda = 1,
                               max_donor_penalty = 1e4) {
  if (nrow(donor_diagnostics) == 0L) return(numeric())

  reliability <- donor_diagnostics$robust_donor_weight
  names(reliability) <- donor_diagnostics$donor

  if (!isTRUE(robust_donors) || donor_penalty_lambda <= 0) {
    penalties <- rep(0, length(reliability))
    names(penalties) <- names(reliability)
    return(penalties)
  }

  penalties <- donor_penalty_lambda * (1 / pmax(reliability, .Machine$double.eps) - 1)
  penalties <- pmin(max_donor_penalty, pmax(0, penalties))
  penalties[!is.finite(penalties)] <- max_donor_penalty
  names(penalties) <- names(reliability)
  penalties
}

sc_append_donor_penalty_rows <- function(y,
                                         x,
                                         donor_penalties) {
  donor_penalties <- as.numeric(donor_penalties)
  if (length(donor_penalties) == 0L || all(donor_penalties <= 0)) {
    return(list(y = y, x = x))
  }

  keep <- donor_penalties > 0
  if (!any(keep)) return(list(y = y, x = x))

  penalty_rows <- matrix(0, nrow = sum(keep), ncol = ncol(x))
  penalty_rows[cbind(seq_len(sum(keep)), which(keep))] <- sqrt(donor_penalties[keep])

  list(
    y = c(y, rep(0, sum(keep))),
    x = rbind(x, penalty_rows)
  )
}
