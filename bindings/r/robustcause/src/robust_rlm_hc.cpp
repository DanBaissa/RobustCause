#include "robust/robust_rlm_hc.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace robust {
namespace {

struct RlmFitState {
  vec beta;
  vec fitted;
  vec resid;
  vec resid_work;
  vec psi_weights;
  vec weights_hc;
  vec hat;
  double scale = std::numeric_limits<double>::quiet_NaN();
  bool converged = false;
  int iterations = 0;
};

inline double normal_pdf(double x) {
  static constexpr double inv_sqrt_2pi = 0.39894228040143267794;
  return inv_sqrt_2pi * std::exp(-0.5 * x * x);
}

inline double normal_cdf(double x) {
  return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

inline double proposal2_gamma(double k) {
  const double phi = normal_pdf(k);
  const double Phi = normal_cdf(k);
  const double val = 2.0 * (Phi - 0.5 - k * phi + k * k * (1.0 - Phi));
  return std::max(val, 1e-8);
}

inline double median_of_vector(std::vector<double> x) {
  if (x.empty()) {
    throw std::invalid_argument("median_of_vector: empty input");
  }
  const std::size_t n = x.size();
  const std::size_t mid = n / 2;
  std::nth_element(x.begin(), x.begin() + mid, x.end());
  double med = x[mid];
  if (n % 2 == 0) {
    auto left_max_it = std::max_element(x.begin(), x.begin() + mid);
    med = 0.5 * (med + *left_max_it);
  }
  return med;
}

inline double mad_scale(const vec& x) {
  std::vector<double> vals(x.begin(), x.end());
  const double med = median_of_vector(vals);

  std::vector<double> devs(static_cast<std::size_t>(x.n_elem));
  for (arma::uword i = 0; i < x.n_elem; ++i) {
    devs[static_cast<std::size_t>(i)] = std::abs(x[i] - med);
  }

  double mad = median_of_vector(devs) * 1.4826;
  if (!(mad > 0.0) || !std::isfinite(mad)) {
    mad = std::sqrt(std::max(1e-16, arma::dot(x, x) /
                                        std::max(1.0, static_cast<double>(x.n_elem - 1))));
  }
  return mad;
}

inline vec sanitize_weights(vec w, double min_weight) {
  w.transform([&](double wi) {
    if (!std::isfinite(wi)) return min_weight;
    return std::max(min_weight, wi);
  });
  return w;
}

inline vec sanitize_prior_weights(const vec& w, arma::uword n) {
  if (w.n_elem == 0) {
    return vec(n, arma::fill::ones);
  }
  if (w.n_elem != n) {
    throw std::invalid_argument("fit_rlm: prior_weights length must match n");
  }
  vec out = w;
  out.transform([](double wi) {
    if (!std::isfinite(wi) || wi <= 0.0) {
      throw std::invalid_argument("fit_rlm: prior_weights must be finite and positive");
    }
    return wi;
  });
  return out;
}

inline vec solve_wls_active(const mat& X,
                            const vec& y,
                            const vec& w,
                            double min_weight,
                            double ridge) {
  const arma::uvec active = arma::find(w > min_weight);
  if (active.n_elem < X.n_cols) {
    throw std::runtime_error("solve_wls_active: too few active observations");
  }

  const mat Xa = X.rows(active);
  const vec ya = y.elem(active);
  const vec sw = arma::sqrt(w.elem(active));
  const mat Xw = Xa.each_col() % sw;
  const vec yw = ya % sw;

  mat XtWX = Xw.t() * Xw;
  XtWX.diag() += ridge;
  const vec XtWy = Xw.t() * yw;

  vec beta;
  const bool ok = arma::solve(beta, XtWX, XtWy,
                              arma::solve_opts::fast + arma::solve_opts::likely_sympd);
  if (!ok || !beta.is_finite()) {
    throw std::runtime_error("solve_wls_active: weighted solve failed");
  }
  return beta;
}

inline vec hatvalues_weighted(const mat& X,
                              const vec& w,
                              double min_weight,
                              double ridge) {
  const arma::uvec active = arma::find(w > min_weight);
  if (active.n_elem == 0) {
    throw std::runtime_error("hatvalues_weighted: no active observations");
  }

  const mat Xa = X.rows(active);
  const vec sw = arma::sqrt(w.elem(active));
  const mat Z = Xa.each_col() % sw;

  mat XtWX = Z.t() * Z;
  XtWX.diag() += ridge;

  mat XtWX_inv;
  const bool ok = arma::inv_sympd(XtWX_inv, XtWX);
  if (!ok || !XtWX_inv.is_finite()) {
    throw std::runtime_error("hatvalues_weighted: inverse solve failed");
  }

  vec h(X.n_rows, arma::fill::zeros);
  for (arma::uword j = 0; j < active.n_elem; ++j) {
    const arma::uword i = active[j];
    const arma::rowvec zi = X.row(i) * std::sqrt(w[i]);
    h[i] = arma::as_scalar(zi * XtWX_inv * zi.t());
  }
  return h;
}

inline double robust_scale(const vec& resid, ScaleEstimate scale_est, double k2) {
  double scale = mad_scale(resid);
  if (scale_est == ScaleEstimate::kMAD) {
    return scale;
  }

  const double gamma = proposal2_gamma(k2);
  for (int iter = 0; iter < 100; ++iter) {
    const double cap = k2 * scale;
    vec clipped = arma::square(resid);
    clipped.transform([cap](double ri2) { return std::min(ri2, cap * cap); });
    const double updated = std::sqrt(std::max(1e-16, arma::mean(clipped) / gamma));
    if (!std::isfinite(updated) || !(updated > 0.0)) {
      break;
    }
    if (std::abs(updated - scale) <= 1e-8 * std::max(1.0, scale)) {
      scale = updated;
      break;
    }
    scale = updated;
  }
  return scale;
}

inline vec make_case_base_weights(const vec& prior_weights, WeightMethod wt_method, arma::uword n) {
  if (wt_method == WeightMethod::kCase) {
    return prior_weights;
  }
  return vec(n, arma::fill::ones);
}

inline mat make_work_matrix(const mat& X, const vec& prior_weights, WeightMethod wt_method) {
  if (wt_method == WeightMethod::kInvVar) {
    return X.each_col() % arma::sqrt(prior_weights);
  }
  return X;
}

inline vec make_work_response(const vec& y, const vec& prior_weights, WeightMethod wt_method) {
  if (wt_method == WeightMethod::kInvVar) {
    return y % arma::sqrt(prior_weights);
  }
  return y;
}

inline vec make_hc_resid(const vec& resid, const vec& prior_weights, WeightMethod wt_method) {
  if (wt_method == WeightMethod::kInvVar) {
    return resid % arma::sqrt(prior_weights);
  }
  return resid;
}

inline vec initial_beta_ols(const mat& X_work,
                            const vec& y_work,
                            const vec& base_weights,
                            double min_weight,
                            double ridge) {
  return solve_wls_active(X_work, y_work, sanitize_weights(base_weights, min_weight), min_weight, ridge);
}

inline vec initial_beta_lts(const mat& X,
                            const vec& y,
                            const RlmControl& ctl) {
  SEstControl s_ctl = ctl.mm_s_control;
  s_ctl.add_intercept = false;
  const SEstResult s_fit = fit_s_estimator(X, y, s_ctl);
  if (!s_fit.coef.is_finite()) {
    throw std::runtime_error("fit_rlm: LTS-style init failed");
  }
  return s_fit.coef;
}

inline double convergence_metric(const vec& beta,
                                 const vec& beta_new,
                                 const vec& resid_work,
                                 const vec& resid_work_new,
                                 const vec& psi_w,
                                 const vec& psi_w_new,
                                 TestVector test_vec) {
  if (test_vec == TestVector::kResid) {
    const double denom = std::max(1.0, arma::abs(resid_work).max());
    return arma::abs(resid_work_new - resid_work).max() / denom;
  }
  if (test_vec == TestVector::kWeight) {
    const double denom = std::max(1.0, arma::abs(psi_w).max());
    return arma::abs(psi_w_new - psi_w).max() / denom;
  }
  const double denom = std::max(1.0, arma::abs(beta).max());
  return arma::abs(beta_new - beta).max() / denom;
}

inline vec final_hc_weights(const vec& psi_w,
                            const vec& prior_weights,
                            WeightMethod wt_method) {
  if (wt_method == WeightMethod::kCase) {
    return psi_w % prior_weights;
  }
  return psi_w;
}

inline RlmFitState run_irls(const mat& X,
                            const vec& y,
                            const mat& X_work,
                            const vec& y_work,
                            const vec& prior_weights,
                            const RlmControl& ctl,
                            const vec& beta_init,
                            double fixed_scale) {
  vec beta = beta_init;
  vec fitted(X.n_rows, arma::fill::zeros);
  vec resid(X.n_rows, arma::fill::zeros);
  vec fitted_work(X_work.n_rows, arma::fill::zeros);
  vec resid_work(X_work.n_rows, arma::fill::zeros);
  vec psi_w(X.n_rows, arma::fill::ones);
  vec psi_w_new(X.n_rows, arma::fill::ones);
  vec combined_w(X.n_rows, arma::fill::ones);
  const vec base_weights = make_case_base_weights(prior_weights, ctl.wt_method, X.n_rows);

  const bool use_fixed_scale = std::isfinite(fixed_scale) && fixed_scale > 0.0;
  double scale = fixed_scale;
  bool converged = false;
  int iterations = 0;

  for (int iter = 0; iter < ctl.maxit; ++iter) {
    fitted = X * beta;
    resid = y - fitted;
    fitted_work = X_work * beta;
    resid_work = y_work - fitted_work;
    scale = use_fixed_scale ? fixed_scale : robust_scale(resid_work, ctl.scale_est, ctl.k2);
    if (!(scale > 0.0) || !std::isfinite(scale)) {
      throw std::runtime_error("fit_rlm: invalid scale estimate");
    }

    psi_w = sanitize_weights(psi_weights(resid_work / scale, ctl.psi, ctl.tuning), ctl.min_weight);
    combined_w = sanitize_weights(psi_w % base_weights, ctl.min_weight);
    const vec beta_new = solve_wls_active(X_work, y_work, combined_w, ctl.min_weight, ctl.ridge);

    const vec resid_work_new = y_work - X_work * beta_new;
    const double scale_new = use_fixed_scale ? fixed_scale : robust_scale(resid_work_new, ctl.scale_est, ctl.k2);
    psi_w_new = sanitize_weights(psi_weights(resid_work_new / scale_new, ctl.psi, ctl.tuning), ctl.min_weight);
    const double diff = convergence_metric(beta, beta_new, resid_work, resid_work_new, psi_w, psi_w_new, ctl.test_vec);

    beta = beta_new;
    iterations = iter + 1;
    if (diff < ctl.acc) {
      converged = true;
      break;
    }
  }

  fitted = X * beta;
  resid = y - fitted;
  resid_work = y_work - X_work * beta;
  scale = use_fixed_scale ? fixed_scale : robust_scale(resid_work, ctl.scale_est, ctl.k2);
  psi_w = sanitize_weights(psi_weights(resid_work / scale, ctl.psi, ctl.tuning), ctl.min_weight);

  RlmFitState out;
  out.beta = beta;
  out.fitted = fitted;
  out.resid = resid;
  out.resid_work = resid_work;
  out.psi_weights = psi_w;
  out.weights_hc = final_hc_weights(psi_w, prior_weights, ctl.wt_method);
  out.hat = hatvalues_weighted(make_work_matrix(X, prior_weights, ctl.wt_method), out.weights_hc, ctl.min_weight, ctl.ridge);
  out.scale = scale;
  out.converged = converged;
  out.iterations = iterations;
  return out;
}

inline vec initial_beta(const mat& X,
                        const vec& y,
                        const mat& X_work,
                        const vec& y_work,
                        const vec& prior_weights,
                        const RlmControl& ctl,
                        double* fixed_scale) {
  if (ctl.method == RlmMethod::kMM) {
    SEstControl s_ctl = ctl.mm_s_control;
    s_ctl.add_intercept = false;
    const SEstResult s_fit = fit_s_estimator(X, y, s_ctl);
    if (!s_fit.coef.is_finite() || !(s_fit.scale > 0.0) || !std::isfinite(s_fit.scale)) {
      throw std::runtime_error("fit_rlm: MM start from S-estimator failed");
    }
    *fixed_scale = s_fit.scale;
    return s_fit.coef;
  }

  if (ctl.init_method == InitMethod::kUser) {
    if (ctl.init_coef.n_elem != X.n_cols) {
      throw std::invalid_argument("fit_rlm: user init_coef has wrong length");
    }
    return ctl.init_coef;
  }
  if (ctl.init_method == InitMethod::kLTS) {
    return initial_beta_lts(X, y, ctl);
  }
  return initial_beta_ols(
    X_work,
    y_work,
    make_case_base_weights(prior_weights, ctl.wt_method, X.n_rows),
    ctl.min_weight,
    ctl.ridge
  );
}

inline vec omega_hc(const vec& e_working,
                    const vec& h,
                    double n,
                    double p,
                    HCType type) {
  vec one_minus_h = 1.0 - h;
  one_minus_h.transform([](double x) { return std::max(1e-12, x); });

  switch (type) {
    case HCType::kHC0:
      return arma::square(e_working);
    case HCType::kHC1:
      return arma::square(e_working) * (n / (n - p));
    case HCType::kHC2:
      return arma::square(e_working) / one_minus_h;
    case HCType::kHC3:
      return arma::square(e_working) / arma::square(one_minus_h);
    case HCType::kHC4: {
      vec delta = arma::min(arma::vec(h.n_elem).fill(4.0), n * h / p);
      vec denom = arma::exp(delta % arma::log(one_minus_h));
      denom.transform([](double x) { return std::max(1e-12, x); });
      return arma::square(e_working) / denom;
    }
    case HCType::kHC4m: {
      vec nhp = n * h / p;
      vec delta = arma::min(arma::vec(h.n_elem).fill(1.0), nhp) +
                  arma::min(arma::vec(h.n_elem).fill(1.5), nhp);
      vec denom = arma::exp(delta % arma::log(one_minus_h));
      denom.transform([](double x) { return std::max(1e-12, x); });
      return arma::square(e_working) / denom;
    }
    case HCType::kHC5: {
      const double k = 0.7;
      const double hmax = h.max();
      const double cap = std::max(4.0, n * k * hmax / p);
      vec delta = arma::min(arma::vec(h.n_elem).fill(cap), n * h / p);
      vec denom = arma::exp(0.5 * delta % arma::log(one_minus_h));
      denom.transform([](double x) { return std::max(1e-12, x); });
      return arma::square(e_working) / denom;
    }
    default:
      throw std::invalid_argument("omega_hc: unknown HC type");
  }
}

}  // namespace

vec psi_weights(const vec& u, PsiType psi, double tuning) {
  vec w(u.n_elem, arma::fill::ones);
  constexpr double eps = 1e-14;

  if (psi == PsiType::kHuber) {
    for (arma::uword i = 0; i < u.n_elem; ++i) {
      const double a = std::abs(u[i]);
      w[i] = (a <= tuning || a < eps) ? 1.0 : (tuning / a);
    }
    return w;
  }

  if (psi == PsiType::kHampel) {
    const double a = tuning / 4.0;
    const double b = tuning / 2.0;
    const double c = tuning;
    for (arma::uword i = 0; i < u.n_elem; ++i) {
      const double abs_u = std::abs(u[i]);
      if (abs_u < eps || abs_u <= a) {
        w[i] = 1.0;
      } else if (abs_u <= b) {
        w[i] = a / abs_u;
      } else if (abs_u <= c) {
        w[i] = a * (c - abs_u) / ((c - b) * abs_u);
      } else {
        w[i] = 0.0;
      }
    }
    return w;
  }

  for (arma::uword i = 0; i < u.n_elem; ++i) {
    const double a = std::abs(u[i]);
    if (a < eps) {
      w[i] = 1.0;
    } else if (a < tuning) {
      const double z = u[i] / tuning;
      const double t = 1.0 - z * z;
      w[i] = t * t;
    } else {
      w[i] = 0.0;
    }
  }
  return w;
}

RlmResult fit_rlm(const mat& X_in, const vec& y, const RlmControl& ctl) {
  if (X_in.n_rows != y.n_elem) {
    throw std::invalid_argument("fit_rlm: X rows must match y length");
  }
  if (X_in.n_rows == 0 || X_in.n_cols == 0) {
    throw std::invalid_argument("fit_rlm: X must be non-empty");
  }
  if (!X_in.is_finite() || !y.is_finite()) {
    throw std::invalid_argument("fit_rlm: X and y must be finite");
  }

  const mat X = ctl.add_intercept ? add_intercept(X_in) : X_in;
  const arma::uword n = X.n_rows;
  const arma::uword p = X.n_cols;
  if (n <= p) {
    throw std::invalid_argument("fit_rlm: need n > p");
  }

  const vec prior_weights = sanitize_prior_weights(ctl.prior_weights, n);
  const mat X_work = make_work_matrix(X, prior_weights, ctl.wt_method);
  const vec y_work = make_work_response(y, prior_weights, ctl.wt_method);

  double fixed_scale = std::numeric_limits<double>::quiet_NaN();
  const vec beta_init = initial_beta(X, y, X_work, y_work, prior_weights, ctl, &fixed_scale);
  const RlmFitState fit_state = run_irls(X, y, X_work, y_work, prior_weights, ctl, beta_init, fixed_scale);

  RlmResult out;
  out.coef = fit_state.beta;
  out.fitted = fit_state.fitted;
  out.resid = fit_state.resid;
  out.weights = fit_state.psi_weights;
  out.hat = fit_state.hat;
  out.scale = fit_state.scale;
  out.converged = fit_state.converged;
  out.iterations = fit_state.iterations;
  out.method = ctl.method;
  out.wt_method = ctl.wt_method;
  out.prior_weights = prior_weights;
  out.X = X;
  out.y = y;
  out.X_hc = X_work;
  out.resid_hc = make_hc_resid(fit_state.resid, prior_weights, ctl.wt_method);
  out.weights_hc = fit_state.weights_hc;
  return out;
}

mat vcov_hc(const RlmResult& fit, HCType type, double ridge) {
  const mat& X = fit.X_hc.n_elem > 0 ? fit.X_hc : fit.X;
  const vec& e = fit.resid_hc.n_elem > 0 ? fit.resid_hc : fit.resid;
  const vec& w = fit.weights_hc.n_elem > 0 ? fit.weights_hc : fit.weights;
  const vec& h = fit.hat;

  const double n = static_cast<double>(X.n_rows);
  const double p = static_cast<double>(X.n_cols);
  const vec e_working = w % e;
  const vec omega = omega_hc(e_working, h, n, p, type);

  mat XtWX = X.t() * arma::diagmat(w) * X;
  XtWX.diag() += ridge;

  mat bread;
  const bool ok = arma::inv_sympd(bread, XtWX);
  if (!ok || !bread.is_finite()) {
    throw std::runtime_error("vcov_hc: inverse solve failed");
  }

  const mat meat = X.t() * arma::diagmat(omega) * X;
  return bread * meat * bread;
}

InferenceResult confint_normal(const RlmResult& fit,
                               HCType type,
                               double zcrit,
                               double ridge) {
  InferenceResult out;
  out.vcov = vcov_hc(fit, type, ridge);
  out.se = arma::sqrt(arma::clamp(out.vcov.diag(), 0.0, std::numeric_limits<double>::infinity()));
  out.ci.set_size(fit.coef.n_elem, 2);
  out.ci.col(0) = fit.coef - zcrit * out.se;
  out.ci.col(1) = fit.coef + zcrit * out.se;
  return out;
}

std::string hc_name(HCType type) {
  switch (type) {
    case HCType::kHC0: return "HC0";
    case HCType::kHC1: return "HC1";
    case HCType::kHC2: return "HC2";
    case HCType::kHC3: return "HC3";
    case HCType::kHC4: return "HC4";
    case HCType::kHC4m: return "HC4m";
    case HCType::kHC5: return "HC5";
    default: return "UNKNOWN";
  }
}

std::string rlm_method_name(RlmMethod method) {
  switch (method) {
    case RlmMethod::kM: return "M";
    case RlmMethod::kMM: return "MM";
    default: return "UNKNOWN";
  }
}

std::string psi_name(PsiType psi) {
  switch (psi) {
    case PsiType::kHuber: return "huber";
    case PsiType::kTukeyBisquare: return "tukey_bisquare";
    case PsiType::kHampel: return "hampel";
    default: return "UNKNOWN";
  }
}

std::string weight_method_name(WeightMethod method) {
  switch (method) {
    case WeightMethod::kInvVar: return "inv.var";
    case WeightMethod::kCase: return "case";
    default: return "UNKNOWN";
  }
}

std::string init_method_name(InitMethod method) {
  switch (method) {
    case InitMethod::kOLS: return "ls";
    case InitMethod::kLTS: return "lts";
    case InitMethod::kUser: return "user";
    default: return "UNKNOWN";
  }
}

std::string scale_estimate_name(ScaleEstimate method) {
  switch (method) {
    case ScaleEstimate::kMAD: return "MAD";
    case ScaleEstimate::kHuber: return "Huber";
    case ScaleEstimate::kProposal2: return "proposal 2";
    default: return "UNKNOWN";
  }
}

std::string test_vector_name(TestVector method) {
  switch (method) {
    case TestVector::kCoef: return "coef";
    case TestVector::kResid: return "resid";
    case TestVector::kWeight: return "w";
    default: return "UNKNOWN";
  }
}

}  // namespace robust
