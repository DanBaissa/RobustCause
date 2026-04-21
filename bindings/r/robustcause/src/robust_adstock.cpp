#include "robust/robust_adstock.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <stdexcept>

namespace robust {
namespace {

double tukey_rho(double u, double c) {
  const double abs_u = std::abs(u);
  if (abs_u >= c) {
    return (c * c) / 6.0;
  }
  const double t = 1.0 - (u / c) * (u / c);
  return (c * c / 6.0) * (1.0 - t * t * t);
}

double tukey_weight(double u, double c) {
  const double abs_u = std::abs(u);
  if (abs_u >= c) {
    return 0.0;
  }
  const double t = 1.0 - (u / c) * (u / c);
  return t * t;
}

void validate_signal(const std::vector<double>& signal) {
  if (signal.empty()) {
    throw std::invalid_argument("signal must not be empty");
  }
}

void validate_design(const std::vector<std::vector<double>>& X,
                     const std::vector<double>* y = nullptr) {
  if (X.empty()) {
    throw std::invalid_argument("X must not be empty");
  }
  const size_t p = X.front().size();
  if (p == 0) {
    throw std::invalid_argument("X must have at least one column");
  }
  for (const auto& row : X) {
    if (row.size() != p) {
      throw std::invalid_argument("X rows must have equal length");
    }
  }
  if (y != nullptr && y->size() != X.size()) {
    throw std::invalid_argument("X and y must have the same number of rows");
  }
}

bool solve_linear_system(std::vector<std::vector<double>> A,
                         std::vector<double> b,
                         std::vector<double>& x) {
  const int p = static_cast<int>(A.size());
  x = std::vector<double>(p, 0.0);

  for (int col = 0; col < p; ++col) {
    int pivot = col;
    double best = std::abs(A[col][col]);
    for (int row = col + 1; row < p; ++row) {
      const double candidate = std::abs(A[row][col]);
      if (candidate > best) {
        best = candidate;
        pivot = row;
      }
    }
    if (best < 1e-12 || !std::isfinite(best)) {
      return false;
    }
    if (pivot != col) {
      std::swap(A[pivot], A[col]);
      std::swap(b[pivot], b[col]);
    }
    const double diag = A[col][col];
    for (int k = col; k < p; ++k) {
      A[col][k] /= diag;
    }
    b[col] /= diag;
    for (int row = 0; row < p; ++row) {
      if (row == col) {
        continue;
      }
      const double factor = A[row][col];
      for (int k = col; k < p; ++k) {
        A[row][k] -= factor * A[col][k];
      }
      b[row] -= factor * b[col];
    }
  }
  x = std::move(b);
  return true;
}

AdstockRegressionResult weighted_ls(const std::vector<std::vector<double>>& X,
                                    const std::vector<double>& y,
                                    const std::vector<double>& w) {
  const int n = static_cast<int>(X.size());
  const int p = static_cast<int>(X.front().size());

  std::vector<std::vector<double>> A(p, std::vector<double>(p, 0.0));
  std::vector<double> b(p, 0.0);

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < p; ++j) {
      b[j] += w[i] * X[i][j] * y[i];
      for (int k = 0; k < p; ++k) {
        A[j][k] += w[i] * X[i][j] * X[i][k];
      }
    }
  }

  std::vector<double> coef;
  const bool ok = solve_linear_system(A, b, coef);
  return {coef, ok, std::numeric_limits<double>::quiet_NaN()};
}

std::vector<double> residuals_from_beta(const std::vector<std::vector<double>>& X,
                                        const std::vector<double>& y,
                                        const std::vector<double>& beta) {
  const int n = static_cast<int>(X.size());
  const int p = static_cast<int>(X.front().size());
  std::vector<double> resid(n, 0.0);
  for (int i = 0; i < n; ++i) {
    double xb = 0.0;
    for (int j = 0; j < p; ++j) {
      xb += X[i][j] * beta[j];
    }
    resid[i] = y[i] - xb;
  }
  return resid;
}

double solve_s_scale_from_residuals(const std::vector<double>& resid,
                                    double c,
                                    double b_target) {
  double base = mad(resid) / 0.67448975;
  if (!std::isfinite(base) || base < 1e-8) {
    double mean_abs = 0.0;
    for (double value : resid) {
      mean_abs += std::abs(value);
    }
    base = std::max(mean_abs / std::max<size_t>(1, resid.size()), 1.0);
  }

  const auto objective = [&](double scale) {
    double acc = 0.0;
    for (double value : resid) {
      acc += tukey_rho(value / scale, c);
    }
    return acc / static_cast<double>(resid.size()) - b_target;
  };

  double lo = std::max(base * 1e-6, 1e-8);
  double hi = std::max(base, 1e-4);
  double f_hi = objective(hi);

  int expand = 0;
  while (f_hi > 0.0 && expand < 25) {
    hi *= 2.0;
    f_hi = objective(hi);
    ++expand;
  }
  if (!std::isfinite(f_hi) || f_hi > 0.0) {
    return hi;
  }

  for (int iter = 0; iter < 35; ++iter) {
    const double mid = 0.5 * (lo + hi);
    const double f_mid = objective(mid);
    if (!std::isfinite(f_mid)) {
      lo = mid;
      continue;
    }
    if (f_mid > 0.0) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  return 0.5 * (lo + hi);
}

AdstockRegressionResult local_s_refine(const std::vector<std::vector<double>>& X,
                                       const std::vector<double>& y,
                                       const std::vector<double>& beta_init,
                                       int max_iter) {
  const double c = 1.547;
  const double rho_max = (c * c) / 6.0;
  const double b_target = 0.5 * rho_max;

  std::vector<double> beta = beta_init;
  double best_scale = std::numeric_limits<double>::infinity();

  for (int iter = 0; iter < max_iter; ++iter) {
    const std::vector<double> resid = residuals_from_beta(X, y, beta);
    const double scale = solve_s_scale_from_residuals(resid, c, b_target);
    if (!std::isfinite(scale) || scale <= 0.0) {
      break;
    }
    best_scale = scale;

    std::vector<double> w(resid.size(), 0.0);
    for (size_t i = 0; i < resid.size(); ++i) {
      w[i] = tukey_weight(resid[i] / scale, c);
    }

    AdstockRegressionResult fit = weighted_ls(X, y, w);
    if (!fit.ok) {
      break;
    }

    double max_diff = 0.0;
    for (size_t j = 0; j < beta.size(); ++j) {
      max_diff = std::max(max_diff, std::abs(fit.coef[j] - beta[j]));
    }
    beta = std::move(fit.coef);
    if (max_diff < 1e-7) {
      break;
    }
  }

  return {beta, std::isfinite(best_scale), best_scale};
}

std::vector<double> subset_fit(const std::vector<std::vector<double>>& X,
                               const std::vector<double>& y,
                               const std::vector<int>& idx) {
  const int p = static_cast<int>(X.front().size());
  std::vector<std::vector<double>> A(p, std::vector<double>(p, 0.0));
  std::vector<double> b(p, 0.0);
  for (int r = 0; r < p; ++r) {
    for (int c = 0; c < p; ++c) {
      A[r][c] = X[idx[r]][c];
    }
    b[r] = y[idx[r]];
  }
  std::vector<double> beta;
  if (!solve_linear_system(A, b, beta)) {
    return {};
  }
  return beta;
}

double apply_increment_rule(double x, const AdstockConfig& cfg) {
  switch (cfg.increment_method) {
    case AdstockIncrementMethod::kPlain:
      return x;
    case AdstockIncrementMethod::kHuber:
      return huber_clip(x, cfg.clip_c);
    case AdstockIncrementMethod::kTanh:
      return tanh_bound(x, cfg.clip_c);
    case AdstockIncrementMethod::kSoftsign:
      return softsign_bound(x, cfg.clip_c);
    case AdstockIncrementMethod::kAdaptiveClip:
      return std::min(x, cfg.adaptive_upper);
  }
  return x;
}

std::vector<double> preclean_signal_generic(const std::vector<double>& signal,
                                            const std::vector<std::vector<double>>& X,
                                            const AdstockPrecleanConfig& cfg,
                                            bool use_mm) {
  validate_signal(signal);
  validate_design(X, &signal);

  AdstockRegressionResult fit =
    use_mm ? adstock_mm_regression_from_s(X,
                                          signal,
                                          cfg.s_starts,
                                          cfg.s_max_iter,
                                          cfg.mm_tukey_c,
                                          cfg.mm_max_iter,
                                          cfg.seed)
           : adstock_s_regression(X, signal, cfg.s_starts, cfg.s_max_iter, cfg.seed);

  if (!fit.ok || !std::isfinite(fit.scale) || fit.scale <= 0.0) {
    return signal;
  }

  const std::vector<double> yhat = adstock_fitted_values(X, fit.coef);
  const double clip_limit = cfg.residual_clip_multiplier * fit.scale;
  std::vector<double> clean(signal.size(), 0.0);
  for (size_t i = 0; i < signal.size(); ++i) {
    clean[i] = yhat[i] + huber_clip(signal[i] - yhat[i], clip_limit);
    if (cfg.nonnegative) {
      clean[i] = std::max(0.0, clean[i]);
    }
  }
  return clean;
}

}  // namespace

double mean(const std::vector<double>& x) {
  validate_signal(x);
  return std::accumulate(x.begin(), x.end(), 0.0) / static_cast<double>(x.size());
}

double median(std::vector<double> x) {
  validate_signal(x);
  const size_t n = x.size();
  std::nth_element(x.begin(), x.begin() + n / 2, x.end());
  double med = x[n / 2];
  if (n % 2 == 0) {
    std::nth_element(x.begin(), x.begin() + n / 2 - 1, x.end());
    med = 0.5 * (med + x[n / 2 - 1]);
  }
  return med;
}

double mad(const std::vector<double>& x) {
  const double med = median(x);
  std::vector<double> abs_dev(x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    abs_dev[i] = std::abs(x[i] - med);
  }
  return median(std::move(abs_dev));
}

double correlation(const std::vector<double>& x, const std::vector<double>& y) {
  validate_signal(x);
  if (x.size() != y.size()) {
    throw std::invalid_argument("x and y must have the same length");
  }
  const double mx = mean(x);
  const double my = mean(y);
  double num = 0.0;
  double dx = 0.0;
  double dy = 0.0;
  for (size_t i = 0; i < x.size(); ++i) {
    const double a = x[i] - mx;
    const double b = y[i] - my;
    num += a * b;
    dx += a * a;
    dy += b * b;
  }
  if (dx <= 0.0 || dy <= 0.0) {
    return 0.0;
  }
  return num / std::sqrt(dx * dy);
}

double log1p_safe(double x) {
  return std::log(1.0 + std::max(0.0, x));
}

double huber_clip(double x, double c) {
  if (std::abs(x) <= c) {
    return x;
  }
  return c * (x > 0.0 ? 1.0 : -1.0);
}

double tanh_bound(double x, double c) {
  if (c <= 0.0) {
    return x;
  }
  return c * std::tanh(x / c);
}

double softsign_bound(double x, double c) {
  if (c <= 0.0) {
    return x;
  }
  return (c * x) / (c + std::abs(x));
}

double adaptive_upper_from_signal(const std::vector<double>& signal, double k) {
  validate_signal(signal);
  const double med = median(signal);
  const double m = mad(signal);
  return med + k * 1.4826 * std::max(m, 1e-8);
}

std::string increment_method_name(AdstockIncrementMethod method) {
  switch (method) {
    case AdstockIncrementMethod::kPlain:
      return "plain";
    case AdstockIncrementMethod::kHuber:
      return "huber";
    case AdstockIncrementMethod::kTanh:
      return "tanh";
    case AdstockIncrementMethod::kSoftsign:
      return "softsign";
    case AdstockIncrementMethod::kAdaptiveClip:
      return "adaptive_clip";
  }
  return "plain";
}

AdstockRegressionResult adstock_ols(const std::vector<std::vector<double>>& X,
                                    const std::vector<double>& y) {
  validate_design(X, &y);
  return weighted_ls(X, y, std::vector<double>(y.size(), 1.0));
}

AdstockRegressionResult adstock_huber_regression(const std::vector<std::vector<double>>& X,
                                                 const std::vector<double>& y,
                                                 double c,
                                                 int max_iter) {
  validate_design(X, &y);
  AdstockRegressionResult fit = adstock_ols(X, y);
  if (!fit.ok) {
    return fit;
  }

  std::vector<double> beta = fit.coef;
  const int n = static_cast<int>(X.size());
  const int p = static_cast<int>(X.front().size());
  double scale = 1.0;

  for (int iter = 0; iter < max_iter; ++iter) {
    std::vector<double> resid(n, 0.0);
    for (int i = 0; i < n; ++i) {
      double xb = 0.0;
      for (int j = 0; j < p; ++j) {
        xb += X[i][j] * beta[j];
      }
      resid[i] = y[i] - xb;
    }

    scale = mad(resid) / 0.6745;
    if (scale < 1e-8 || !std::isfinite(scale)) {
      scale = 1.0;
    }

    std::vector<double> w(n, 1.0);
    for (int i = 0; i < n; ++i) {
      const double u = resid[i] / scale;
      w[i] = (std::abs(u) <= c) ? 1.0 : (c / std::abs(u));
    }

    AdstockRegressionResult new_fit = weighted_ls(X, y, w);
    if (!new_fit.ok) {
      break;
    }

    double max_diff = 0.0;
    for (int j = 0; j < p; ++j) {
      max_diff = std::max(max_diff, std::abs(new_fit.coef[j] - beta[j]));
    }
    beta = std::move(new_fit.coef);
    if (max_diff < 1e-8) {
      return {beta, true, scale};
    }
  }
  return {beta, true, scale};
}

AdstockRegressionResult adstock_s_regression(const std::vector<std::vector<double>>& X,
                                             const std::vector<double>& y,
                                             int n_starts,
                                             int max_iter,
                                             unsigned int seed) {
  validate_design(X, &y);
  const int n = static_cast<int>(X.size());
  const int p = static_cast<int>(X.front().size());
  if (n < p) {
    return {{}, false, std::numeric_limits<double>::quiet_NaN()};
  }

  std::mt19937 rng(seed);
  std::vector<int> all_idx(n);
  std::iota(all_idx.begin(), all_idx.end(), 0);

  AdstockRegressionResult best{{}, false, std::numeric_limits<double>::infinity()};
  for (int start = 0; start < n_starts; ++start) {
    std::shuffle(all_idx.begin(), all_idx.end(), rng);
    std::vector<int> subset(all_idx.begin(), all_idx.begin() + p);
    std::vector<double> beta0 = subset_fit(X, y, subset);
    if (beta0.empty()) {
      continue;
    }

    AdstockRegressionResult refined = local_s_refine(X, y, beta0, max_iter);
    if (refined.ok && refined.scale < best.scale) {
      best = std::move(refined);
    }
  }
  if (!best.ok) {
    return {{}, false, std::numeric_limits<double>::quiet_NaN()};
  }
  return best;
}

AdstockRegressionResult adstock_mm_regression_from_s(const std::vector<std::vector<double>>& X,
                                                     const std::vector<double>& y,
                                                     int s_starts,
                                                     int s_max_iter,
                                                     double c_m,
                                                     int mm_max_iter,
                                                     unsigned int seed) {
  validate_design(X, &y);
  AdstockRegressionResult s_fit = adstock_s_regression(X, y, s_starts, s_max_iter, seed);
  if (!s_fit.ok || !std::isfinite(s_fit.scale) || s_fit.scale <= 0.0) {
    return {s_fit.coef, false, s_fit.scale};
  }

  std::vector<double> beta = s_fit.coef;
  const double scale = s_fit.scale;

  for (int iter = 0; iter < mm_max_iter; ++iter) {
    std::vector<double> resid = residuals_from_beta(X, y, beta);
    std::vector<double> w(y.size(), 0.0);
    for (size_t i = 0; i < y.size(); ++i) {
      w[i] = tukey_weight(resid[i] / scale, c_m);
    }

    AdstockRegressionResult fit = weighted_ls(X, y, w);
    if (!fit.ok) {
      break;
    }

    double max_diff = 0.0;
    for (size_t j = 0; j < beta.size(); ++j) {
      max_diff = std::max(max_diff, std::abs(fit.coef[j] - beta[j]));
    }
    beta = std::move(fit.coef);
    if (max_diff < 1e-7) {
      break;
    }
  }

  return {beta, true, scale};
}

std::vector<double> adstock_fitted_values(const std::vector<std::vector<double>>& X,
                                          const std::vector<double>& beta) {
  validate_design(X);
  if (X.front().size() != beta.size()) {
    throw std::invalid_argument("beta length must match the number of X columns");
  }

  std::vector<double> out(X.size(), 0.0);
  for (size_t i = 0; i < X.size(); ++i) {
    for (size_t j = 0; j < beta.size(); ++j) {
      out[i] += X[i][j] * beta[j];
    }
  }
  return out;
}

std::vector<double> preclean_signal_s(const std::vector<double>& signal,
                                      const std::vector<std::vector<double>>& X,
                                      const AdstockPrecleanConfig& cfg) {
  return preclean_signal_generic(signal, X, cfg, false);
}

std::vector<double> preclean_signal_mm(const std::vector<double>& signal,
                                       const std::vector<std::vector<double>>& X,
                                       const AdstockPrecleanConfig& cfg) {
  return preclean_signal_generic(signal, X, cfg, true);
}

std::vector<double> build_adstock(const std::vector<double>& signal,
                                  const std::vector<int>& unit_ids,
                                  const AdstockConfig& cfg) {
  validate_signal(signal);
  if (signal.size() != unit_ids.size()) {
    throw std::invalid_argument("signal and unit_ids must have the same length");
  }

  AdstockConfig resolved = cfg;
  if (resolved.increment_method == AdstockIncrementMethod::kAdaptiveClip &&
      resolved.adaptive_upper <= 0.0) {
    resolved.adaptive_upper = adaptive_upper_from_signal(signal);
  }

  std::vector<double> out(signal.size(), 0.0);
  std::map<int, double> prev_by_unit;

  for (size_t i = 0; i < signal.size(); ++i) {
    const int unit = unit_ids[i];
    const double prev = prev_by_unit.count(unit) ? prev_by_unit[unit] : 0.0;
    const double increment = apply_increment_rule(signal[i], resolved);
    out[i] = resolved.rho * prev + increment;
    prev_by_unit[unit] = out[i];
  }
  return out;
}

std::vector<double> build_robust_adstock(const std::vector<double>& signal,
                                         const std::vector<int>& unit_ids,
                                         const std::vector<std::vector<double>>& X_preclean,
                                         const RobustAdstockConfig& cfg) {
  std::vector<double> cleaned = signal;
  if (cfg.preclean.enabled) {
    cleaned = cfg.preclean.use_mm ? preclean_signal_mm(signal, X_preclean, cfg.preclean)
                                  : preclean_signal_s(signal, X_preclean, cfg.preclean);
  }
  return build_adstock(cleaned, unit_ids, cfg.adstock);
}

}  // namespace robust
