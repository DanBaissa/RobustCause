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
  vec weights;
  vec hat;
  double scale = std::numeric_limits<double>::quiet_NaN();
  bool converged = false;
  int iterations = 0;
};

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

inline vec sanitize_weights(vec w, double min_weight) {
  w.transform([&](double wi) {
    if (!std::isfinite(wi)) return min_weight;
    return std::max(min_weight, wi);
  });
  return w;
}

inline RlmFitState run_irls(const mat& X,
                            const vec& y,
                            const RlmControl& ctl,
                            const vec& beta_init,
                            double fixed_scale) {
  vec beta = beta_init;
  vec fitted(X.n_rows, arma::fill::zeros);
  vec resid(X.n_rows, arma::fill::zeros);
  vec w(X.n_rows, arma::fill::ones);
  vec u(X.n_rows, arma::fill::zeros);

  const bool use_fixed_scale = std::isfinite(fixed_scale) && fixed_scale > 0.0;
  double scale = fixed_scale;
  bool converged = false;
  int iterations = 0;

  for (int iter = 0; iter < ctl.maxit; ++iter) {
    fitted = X * beta;
    resid = y - fitted;
    scale = use_fixed_scale ? fixed_scale : mad_scale(resid);
    if (!(scale > 0.0) || !std::isfinite(scale)) {
      throw std::runtime_error("fit_rlm: invalid scale estimate");
    }

    u = resid / scale;
    w = sanitize_weights(psi_weights(u, ctl.psi, ctl.tuning), ctl.min_weight);

    const vec beta_new = solve_wls_active(X, y, w, ctl.min_weight, ctl.ridge);
    const double denom = std::max(1.0, arma::abs(beta).max());
    const double diff = arma::abs(beta_new - beta).max() / denom;
    beta = beta_new;
    iterations = iter + 1;

    if (diff < ctl.tol) {
      converged = true;
      break;
    }
  }

  fitted = X * beta;
  resid = y - fitted;
  scale = use_fixed_scale ? fixed_scale : mad_scale(resid);
  if (!(scale > 0.0) || !std::isfinite(scale)) {
    throw std::runtime_error("fit_rlm: invalid final scale estimate");
  }

  u = resid / scale;
  w = sanitize_weights(psi_weights(u, ctl.psi, ctl.tuning), ctl.min_weight);

  RlmFitState out;
  out.beta = beta;
  out.fitted = fitted;
  out.resid = resid;
  out.weights = w;
  out.hat = hatvalues_weighted(X, w, ctl.min_weight, ctl.ridge);
  out.scale = scale;
  out.converged = converged;
  out.iterations = iterations;
  return out;
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

  vec beta_init;
  double fixed_scale = std::numeric_limits<double>::quiet_NaN();

  if (ctl.method == RlmMethod::kMM) {
    SEstControl s_ctl = ctl.mm_s_control;
    s_ctl.add_intercept = false;
    const SEstResult s_fit = fit_s_estimator(X, y, s_ctl);
    if (!s_fit.coef.is_finite() || !(s_fit.scale > 0.0) || !std::isfinite(s_fit.scale)) {
      throw std::runtime_error("fit_rlm: MM start from S-estimator failed");
    }
    beta_init = s_fit.coef;
    fixed_scale = s_fit.scale;
  } else {
    const bool ok_ols = arma::solve(beta_init, X, y);
    if (!ok_ols || !beta_init.is_finite()) {
      throw std::runtime_error("fit_rlm: OLS init failed");
    }
  }

  const RlmFitState fit_state = run_irls(X, y, ctl, beta_init, fixed_scale);

  RlmResult out;
  out.coef = fit_state.beta;
  out.fitted = fit_state.fitted;
  out.resid = fit_state.resid;
  out.weights = fit_state.weights;
  out.hat = fit_state.hat;
  out.scale = fit_state.scale;
  out.converged = fit_state.converged;
  out.iterations = fit_state.iterations;
  out.method = ctl.method;
  out.X = X;
  out.y = y;
  return out;
}

mat vcov_hc(const RlmResult& fit, HCType type, double ridge) {
  const mat& X = fit.X;
  const vec& e = fit.resid;
  const vec& w = fit.weights;
  const vec& h = fit.hat;

  const double n = static_cast<double>(X.n_rows);
  const double p = static_cast<double>(X.n_cols);

  vec e_working = w % e;
  vec omega = omega_hc(e_working, h, n, p, type);

  mat XtWX = X.t() * arma::diagmat(w) * X;
  XtWX.diag() += ridge;

  mat bread;
  const bool ok = arma::inv_sympd(bread, XtWX);
  if (!ok || !bread.is_finite()) {
    throw std::runtime_error("vcov_hc: inverse solve failed");
  }

  mat meat = X.t() * arma::diagmat(omega) * X;
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

}  // namespace robust
