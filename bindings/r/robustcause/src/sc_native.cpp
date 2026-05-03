#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

namespace {

arma::vec project_simplex(arma::vec v) {
  const int n = v.n_elem;
  arma::vec u = arma::sort(v, "descend");
  double cssv = 0.0;
  int rho = 0;
  for (int i = 0; i < n; ++i) {
    cssv += u[i];
    const double t = (cssv - 1.0) / static_cast<double>(i + 1);
    if (u[i] - t > 0.0) rho = i + 1;
  }
  const double theta = (arma::sum(u.head(rho)) - 1.0) / static_cast<double>(rho);
  arma::vec out = v - theta;
  out.transform([](double x) { return x > 0.0 ? x : 0.0; });
  const double s = arma::sum(out);
  if (s <= 0.0 || !std::isfinite(s)) return arma::vec(v.n_elem, arma::fill::ones) / v.n_elem;
  return out / s;
}

double rmspe(const arma::vec& x) {
  return std::sqrt(arma::mean(arma::square(x)));
}

arma::vec fit_weighted_sc(const arma::vec& y,
                          const arma::mat& X,
                          const arma::vec& period_weights,
                          int maxit,
                          double tol,
                          double ridge) {
  const int p = X.n_cols;
  arma::vec w(p, arma::fill::ones);
  w /= static_cast<double>(p);
  arma::vec wt = period_weights;
  if (wt.n_elem != y.n_elem) wt = arma::vec(y.n_elem, arma::fill::ones);
  wt = wt / std::max(arma::mean(wt), 1e-12);

  arma::mat Xw = X.each_col() % arma::sqrt(wt);
  const double lip = 2.0 * arma::norm(Xw.t() * Xw, 2) + 2.0 * ridge + 1e-8;
  const double step = 1.0 / lip;

  for (int it = 0; it < maxit; ++it) {
    arma::vec resid = y - X * w;
    arma::vec grad = -2.0 * X.t() * (wt % resid) / static_cast<double>(y.n_elem) + 2.0 * ridge * w;
    arma::vec w_new = project_simplex(w - step * grad);
    if (arma::max(arma::abs(w_new - w)) < tol) {
      w = w_new;
      break;
    }
    w = w_new;
  }
  return w;
}

arma::vec tukey_weights(const arma::vec& residuals, double c, double min_weight) {
  double med = arma::median(residuals);
  double scale = arma::median(arma::abs(residuals - med)) / 0.6745;
  if (!std::isfinite(scale) || scale <= 1e-12) scale = rmspe(residuals) + 1e-12;
  arma::vec u = residuals / scale;
  arma::vec out(u.n_elem, arma::fill::zeros);
  for (arma::uword i = 0; i < u.n_elem; ++i) {
    double a = std::abs(u[i] / c);
    if (a < 1.0) out[i] = std::pow(1.0 - a * a, 2.0);
    if (out[i] < min_weight) out[i] = min_weight;
  }
  return out;
}

Rcpp::List build_sc_result(const arma::vec& w,
                           const arma::vec& treated_pre,
                           const arma::mat& donors_pre,
                           const arma::vec& treated_post,
                           const arma::mat& donors_post,
                           bool converged,
                           int iterations) {
  arma::vec synth_pre = donors_pre * w;
  arma::vec synth_post = donors_post * w;
  arma::vec pre_resid = treated_pre - synth_pre;
  arma::vec post_gaps = treated_post - synth_post;
  double h = arma::sum(arma::square(w));
  return Rcpp::List::create(
    Rcpp::Named("weights") = Rcpp::wrap(w),
    Rcpp::Named("synthetic_pre") = Rcpp::wrap(synth_pre),
    Rcpp::Named("synthetic_post") = Rcpp::wrap(synth_post),
    Rcpp::Named("pre_residuals") = Rcpp::wrap(pre_resid),
    Rcpp::Named("post_gaps") = Rcpp::wrap(post_gaps),
    Rcpp::Named("pre_rmspe") = rmspe(pre_resid),
    Rcpp::Named("weight_herfindahl") = h,
    Rcpp::Named("max_weight") = arma::max(w),
    Rcpp::Named("effective_donors") = 1.0 / std::max(h, 1e-12),
    Rcpp::Named("converged") = converged,
    Rcpp::Named("iterations") = iterations
  );
}

}  // namespace

extern "C" SEXP rc_r_fit_sc(SEXP treatedPreSEXP,
                            SEXP donorsPreSEXP,
                            SEXP treatedPostSEXP,
                            SEXP donorsPostSEXP,
                            SEXP maxitSEXP,
                            SEXP tolSEXP,
                            SEXP ridgeSEXP,
                            SEXP predictorTreatedSEXP,
                            SEXP predictorDonorsSEXP,
                            SEXP predictorWeightsSEXP) {
  arma::vec treated_pre = Rcpp::as<arma::vec>(treatedPreSEXP);
  arma::mat donors_pre = Rcpp::as<arma::mat>(donorsPreSEXP);
  arma::vec treated_post = Rcpp::as<arma::vec>(treatedPostSEXP);
  arma::mat donors_post = Rcpp::as<arma::mat>(donorsPostSEXP);
  int maxit = Rcpp::as<int>(maxitSEXP);
  double tol = Rcpp::as<double>(tolSEXP);
  double ridge = Rcpp::as<double>(ridgeSEXP);
  arma::vec wt(treated_pre.n_elem, arma::fill::ones);
  arma::vec w = fit_weighted_sc(treated_pre, donors_pre, wt, maxit, tol, ridge);
  return build_sc_result(w, treated_pre, donors_pre, treated_post, donors_post, true, maxit);
}

extern "C" SEXP rc_r_fit_mm_sc(SEXP treatedPreSEXP,
                               SEXP donorsPreSEXP,
                               SEXP treatedPostSEXP,
                               SEXP donorsPostSEXP,
                               SEXP startupMaxitSEXP,
                               SEXP maxitSEXP,
                               SEXP subproblemMaxitSEXP,
                               SEXP tolSEXP,
                               SEXP subproblemTolSEXP,
                               SEXP l1SmoothingSEXP,
                               SEXP tukeyCSEXP,
                               SEXP minScaleSEXP,
                               SEXP ridgeSEXP,
                               SEXP minTimeWeightSEXP,
                               SEXP predictorTreatedSEXP,
                               SEXP predictorDonorsSEXP,
                               SEXP predictorWeightsSEXP) {
  arma::vec treated_pre = Rcpp::as<arma::vec>(treatedPreSEXP);
  arma::mat donors_pre = Rcpp::as<arma::mat>(donorsPreSEXP);
  arma::vec treated_post = Rcpp::as<arma::vec>(treatedPostSEXP);
  arma::mat donors_post = Rcpp::as<arma::mat>(donorsPostSEXP);
  int startup_maxit = Rcpp::as<int>(startupMaxitSEXP);
  int maxit = Rcpp::as<int>(maxitSEXP);
  int subproblem_maxit = Rcpp::as<int>(subproblemMaxitSEXP);
  double tol = Rcpp::as<double>(tolSEXP);
  double sub_tol = Rcpp::as<double>(subproblemTolSEXP);
  double tukey_c = Rcpp::as<double>(tukeyCSEXP);
  double ridge = Rcpp::as<double>(ridgeSEXP);
  double min_time_weight = Rcpp::as<double>(minTimeWeightSEXP);

  arma::vec wt(treated_pre.n_elem, arma::fill::ones);
  arma::vec start_w = fit_weighted_sc(treated_pre, donors_pre, wt, startup_maxit, tol, ridge);
  arma::vec w = start_w;
  int iter = 0;
  for (; iter < maxit; ++iter) {
    arma::vec resid = treated_pre - donors_pre * w;
    arma::vec wt_new = tukey_weights(resid, tukey_c, min_time_weight);
    arma::vec w_new = fit_weighted_sc(treated_pre, donors_pre, wt_new, subproblem_maxit, sub_tol, ridge);
    if (arma::max(arma::abs(w_new - w)) < tol) {
      w = w_new;
      wt = wt_new;
      break;
    }
    w = w_new;
    wt = wt_new;
  }

  Rcpp::List out = build_sc_result(w, treated_pre, donors_pre, treated_post, donors_post, true, iter + 1);
  arma::vec resid = treated_pre - donors_pre * w;
  out["startup_weights"] = Rcpp::wrap(start_w);
  out["robust_time_weights"] = Rcpp::wrap(wt);
  out["downweighted_periods"] = Rcpp::wrap(arma::find(wt < 0.5));
  out["standard_rmspe"] = rmspe(treated_pre - donors_pre * start_w);
  out["startup_rmspe"] = rmspe(treated_pre - donors_pre * start_w);
  out["mm_rmspe"] = rmspe(resid);
  out["effective_periods"] = std::pow(arma::sum(wt), 2.0) / std::max(arma::sum(arma::square(wt)), 1e-12);
  return out;
}

extern "C" SEXP rc_r_fit_sc_placebos(SEXP outcomesSEXP,
                                     SEXP treatmentStartSEXP,
                                     SEXP maxitSEXP,
                                     SEXP tolSEXP,
                                     SEXP ridgeSEXP,
                                     SEXP startupMaxitSEXP,
                                     SEXP subproblemMaxitSEXP,
                                     SEXP subproblemTolSEXP,
                                     SEXP l1SmoothingSEXP,
                                     SEXP tukeyCSEXP,
                                     SEXP minScaleSEXP,
                                     SEXP minTimeWeightSEXP,
                                     SEXP useMmSEXP) {
  arma::mat outcomes = Rcpp::as<arma::mat>(outcomesSEXP);
  int treatment_start = Rcpp::as<int>(treatmentStartSEXP);
  int maxit = Rcpp::as<int>(maxitSEXP);
  double tol = Rcpp::as<double>(tolSEXP);
  double ridge = Rcpp::as<double>(ridgeSEXP);
  int startup_maxit = Rcpp::as<int>(startupMaxitSEXP);
  int subproblem_maxit = Rcpp::as<int>(subproblemMaxitSEXP);
  double sub_tol = Rcpp::as<double>(subproblemTolSEXP);
  double tukey_c = Rcpp::as<double>(tukeyCSEXP);
  double min_time_weight = Rcpp::as<double>(minTimeWeightSEXP);
  bool use_mm = Rcpp::as<bool>(useMmSEXP);

  int n_units = outcomes.n_cols;
  int n_post = outcomes.n_rows - treatment_start + 1;
  arma::mat post_gap_matrix(n_units, n_post, arma::fill::zeros);
  arma::vec pre_rmspe(n_units), post_rmspe(n_units), ratio(n_units), avg_abs_gap(n_units);

  for (int j = 0; j < n_units; ++j) {
    arma::uvec donor_idx(n_units - 1);
    int k = 0;
    for (int u = 0; u < n_units; ++u) if (u != j) donor_idx[k++] = u;
    arma::vec y_pre = outcomes.submat(0, j, treatment_start - 2, j);
    arma::mat x_pre = outcomes.submat(0, 0, treatment_start - 2, n_units - 1).cols(donor_idx);
    arma::vec y_post = outcomes.submat(treatment_start - 1, j, outcomes.n_rows - 1, j);
    arma::mat x_post = outcomes.submat(treatment_start - 1, 0, outcomes.n_rows - 1, n_units - 1).cols(donor_idx);
    arma::vec wt(y_pre.n_elem, arma::fill::ones);
    arma::vec w = fit_weighted_sc(y_pre, x_pre, wt, maxit, tol, ridge);
    if (use_mm) {
      for (int it = 0; it < maxit; ++it) {
        wt = tukey_weights(y_pre - x_pre * w, tukey_c, min_time_weight);
        arma::vec w_new = fit_weighted_sc(y_pre, x_pre, wt, subproblem_maxit, sub_tol, ridge);
        if (arma::max(arma::abs(w_new - w)) < tol) { w = w_new; break; }
        w = w_new;
      }
    }
    arma::vec pre_gap = y_pre - x_pre * w;
    arma::vec post_gap = y_post - x_post * w;
    post_gap_matrix.row(j) = post_gap.t();
    pre_rmspe[j] = std::max(rmspe(pre_gap), 1e-12);
    post_rmspe[j] = rmspe(post_gap);
    ratio[j] = post_rmspe[j] / pre_rmspe[j];
    avg_abs_gap[j] = arma::mean(arma::abs(post_gap));
  }

  return Rcpp::List::create(
    Rcpp::Named("unit_indices") = Rcpp::seq(0, n_units - 1),
    Rcpp::Named("post_gap_matrix") = Rcpp::wrap(post_gap_matrix),
    Rcpp::Named("pre_rmspe") = Rcpp::wrap(pre_rmspe),
    Rcpp::Named("post_rmspe") = Rcpp::wrap(post_rmspe),
    Rcpp::Named("rmspe_ratio") = Rcpp::wrap(ratio),
    Rcpp::Named("avg_abs_post_gap") = Rcpp::wrap(avg_abs_gap)
  );
}
