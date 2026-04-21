#include "robust/robust_s_estimator.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_set>
#include <utility>
#include <vector>

namespace robust {
namespace {

struct StartCandidate {
  vec beta;
  double scale = std::numeric_limits<double>::infinity();
  bool converged = false;
  int iterations = 0;
};

struct ParsedFormula {
  std::string response;
  bool intercept = true;
  std::vector<std::vector<std::string>> terms;
  std::vector<std::string> term_labels;
  std::vector<std::string> variables_used;
};

inline bool all_finite(const vec& x) {
  return x.is_finite();
}

inline bool all_finite(const mat& x) {
  return x.is_finite();
}

inline std::string trim(const std::string& s) {
  const auto begin = s.find_first_not_of(" \t\n\r");
  if (begin == std::string::npos) return "";
  const auto end = s.find_last_not_of(" \t\n\r");
  return s.substr(begin, end - begin + 1);
}

inline std::vector<std::string> split_simple(const std::string& s, char delim) {
  std::vector<std::string> out;
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    item = trim(item);
    if (!item.empty()) out.push_back(item);
  }
  return out;
}

inline std::string join_colon(const std::vector<std::string>& vars) {
  std::string out;
  for (std::size_t i = 0; i < vars.size(); ++i) {
    if (i > 0) out += ":";
    out += vars[i];
  }
  return out;
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

inline bool full_rank(const mat& A) {
  return static_cast<int>(arma::rank(A)) == static_cast<int>(A.n_cols);
}

inline vec solve_subset_ls(const mat& X_sub, const vec& y_sub) {
  vec beta;
  const bool ok = arma::solve(beta, X_sub, y_sub);
  if (!ok || !beta.is_finite()) {
    throw std::runtime_error("solve_subset_ls: solve failed");
  }
  return beta;
}

inline std::vector<int> sample_without_replacement(int n, int k, std::mt19937_64& rng) {
  std::vector<int> idx(static_cast<std::size_t>(n));
  std::iota(idx.begin(), idx.end(), 0);
  std::shuffle(idx.begin(), idx.end(), rng);
  idx.resize(static_cast<std::size_t>(k));
  return idx;
}

inline vec solve_wls_active(const mat& X,
                            const vec& y,
                            const vec& w,
                            double min_weight,
                            double ridge) {
  const uvec active = arma::find(w > min_weight);
  if (active.n_elem < X.n_cols) {
    throw std::runtime_error("solve_wls_active: too few active observations");
  }

  const vec sw = arma::sqrt(w.elem(active));
  const mat Xa = X.rows(active);
  const vec ya = y.elem(active);

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

inline double mean_rho(const vec& r, double s, double c) {
  if (!(s > 0.0) || !std::isfinite(s)) {
    return std::numeric_limits<double>::infinity();
  }

  double acc = 0.0;
  for (arma::uword i = 0; i < r.n_elem; ++i) {
    acc += tukey_rho_normalized(r[i] / s, c);
  }
  return acc / static_cast<double>(r.n_elem);
}

inline StartCandidate refine_s_start(const mat& X,
                                     const vec& y,
                                     vec beta,
                                     const SEstControl& ctl,
                                     int max_iters) {
  StartCandidate out;
  out.beta = beta;

  vec resid(X.n_rows, arma::fill::zeros);
  vec weights(X.n_rows, arma::fill::ones);
  double prev_scale = std::numeric_limits<double>::infinity();

  for (int iter = 0; iter < max_iters; ++iter) {
    resid = y - X * beta;
    if (!resid.is_finite()) {
      throw std::runtime_error("refine_s_start: non-finite residuals");
    }

    const double scale = solve_s_scale(resid, ctl.c, ctl.b, ctl.max_scale_iter, ctl.scale_tol);
    if (!(scale >= 0.0) || !std::isfinite(scale)) {
      throw std::runtime_error("refine_s_start: non-finite scale");
    }

    if (scale == 0.0) {
      out.beta = beta;
      out.scale = 0.0;
      out.converged = true;
      out.iterations = iter + 1;
      return out;
    }

    for (arma::uword i = 0; i < resid.n_elem; ++i) {
      weights[i] = tukey_weight(resid[i] / scale, ctl.c);
    }

    if (!weights.is_finite() || arma::accu(weights > ctl.min_weight) < X.n_cols) {
      throw std::runtime_error("refine_s_start: insufficient positive weights");
    }

    const vec beta_new = solve_wls_active(X, y, weights, ctl.min_weight, ctl.ridge);

    const double beta_denom = std::max(1.0, arma::abs(beta).max());
    const double beta_diff = arma::abs(beta_new - beta).max() / beta_denom;
    const double scale_diff = std::abs(scale - prev_scale) / std::max(1.0, std::abs(prev_scale));

    beta = beta_new;
    prev_scale = scale;
    out.iterations = iter + 1;

    if (beta_diff < ctl.tol && (iter == 0 || scale_diff < ctl.tol)) {
      out.beta = beta;
      out.scale = scale;
      out.converged = true;
      return out;
    }
  }

  resid = y - X * beta;
  out.beta = beta;
  out.scale = solve_s_scale(resid, ctl.c, ctl.b, ctl.max_scale_iter, ctl.scale_tol);
  out.converged = false;
  return out;
}

inline std::vector<std::pair<char, std::string>> split_additive_rhs(const std::string& rhs) {
  std::vector<std::pair<char, std::string>> pieces;
  std::string current;
  char sign = '+';

  for (char ch : rhs) {
    if (ch == '+' || ch == '-') {
      const std::string token = trim(current);
      if (!token.empty()) pieces.push_back({sign, token});
      sign = ch;
      current.clear();
    } else {
      current.push_back(ch);
    }
  }

  const std::string token = trim(current);
  if (!token.empty()) pieces.push_back({sign, token});
  return pieces;
}

inline void add_term(std::vector<std::vector<std::string>>& terms,
                     std::vector<std::string>& labels,
                     std::unordered_set<std::string>& seen,
                     std::vector<std::string> vars) {
  for (auto& v : vars) v = trim(v);
  vars.erase(std::remove_if(vars.begin(), vars.end(), [](const std::string& s) { return s.empty(); }), vars.end());
  if (vars.empty()) return;
  std::sort(vars.begin(), vars.end());
  const std::string label = join_colon(vars);
  if (seen.insert(label).second) {
    terms.push_back(vars);
    labels.push_back(label);
  }
}

inline void expand_star_terms_rec(const std::vector<std::string>& vars,
                                  std::size_t idx,
                                  std::vector<std::string>& current,
                                  std::vector<std::vector<std::string>>& out) {
  if (idx == vars.size()) {
    if (!current.empty()) out.push_back(current);
    return;
  }
  expand_star_terms_rec(vars, idx + 1, current, out);
  current.push_back(vars[idx]);
  expand_star_terms_rec(vars, idx + 1, current, out);
  current.pop_back();
}

inline ParsedFormula parse_formula(const std::string& formula, bool default_intercept) {
  const auto tilde_pos = formula.find('~');
  if (tilde_pos == std::string::npos) {
    throw std::invalid_argument("parse_formula: formula must contain '~'");
  }

  ParsedFormula out;
  out.response = trim(formula.substr(0, tilde_pos));
  if (out.response.empty()) {
    throw std::invalid_argument("parse_formula: missing response variable");
  }

  std::string rhs = trim(formula.substr(tilde_pos + 1));
  if (rhs.empty()) {
    throw std::invalid_argument("parse_formula: missing RHS");
  }

  out.intercept = default_intercept;
  std::unordered_set<std::string> seen_terms;
  std::unordered_set<std::string> seen_vars;

  const auto pieces = split_additive_rhs(rhs);
  for (const auto& piece : pieces) {
    const char sign = piece.first;
    const std::string token = trim(piece.second);
    if (token == "1") {
      out.intercept = (sign == '+');
      continue;
    }
    if (token == "0") {
      out.intercept = false;
      continue;
    }
    if (token.find('(') != std::string::npos || token.find(')') != std::string::npos) {
      throw std::invalid_argument("parse_formula: parentheses are not supported in this lightweight parser");
    }

    if (token.find('*') != std::string::npos) {
      const auto vars = split_simple(token, '*');
      if (vars.empty()) continue;
      std::vector<std::vector<std::string>> expanded;
      std::vector<std::string> current;
      expand_star_terms_rec(vars, 0, current, expanded);
      for (auto& term : expanded) {
        add_term(out.terms, out.term_labels, seen_terms, term);
        for (const auto& var : term) seen_vars.insert(var);
      }
      continue;
    }

    if (token.find(':') != std::string::npos) {
      auto vars = split_simple(token, ':');
      add_term(out.terms, out.term_labels, seen_terms, vars);
      for (const auto& var : vars) seen_vars.insert(var);
      continue;
    }

    add_term(out.terms, out.term_labels, seen_terms, {token});
    seen_vars.insert(token);
  }

  out.variables_used.assign(seen_vars.begin(), seen_vars.end());
  std::sort(out.variables_used.begin(), out.variables_used.end());
  return out;
}

inline vec lookup_column(const NamedColumns& data, const std::string& name) {
  const auto it = data.find(name);
  if (it == data.end()) {
    throw std::invalid_argument("build_model_frame: missing variable '" + name + "'");
  }
  return it->second;
}

inline vec interaction_column(const NamedColumns& data, const std::vector<std::string>& vars) {
  vec out = lookup_column(data, vars.front());
  for (std::size_t j = 1; j < vars.size(); ++j) {
    out %= lookup_column(data, vars[j]);
  }
  return out;
}

}  // namespace

mat add_intercept(const mat& X) {
  mat out(X.n_rows, X.n_cols + 1, arma::fill::ones);
  out.cols(1, X.n_cols) = X;
  return out;
}

double tukey_rho_normalized(double u, double c) {
  const double a = std::abs(u);
  if (a >= c) return 1.0;
  const double z = u / c;
  const double t = 1.0 - z * z;
  return 1.0 - t * t * t;
}

double tukey_weight(double u, double c) {
  const double a = std::abs(u);
  if (a >= c) return 0.0;
  const double z = u / c;
  const double t = 1.0 - z * z;
  return t * t;
}

double solve_s_scale(const vec& resid, double c, double b, int max_iter, double tol) {
  const double max_abs = arma::abs(resid).max();
  if (!(max_abs > 0.0) || !std::isfinite(max_abs)) return 0.0;

  double high = std::max(1.0, max_abs);
  while (mean_rho(resid, high, c) > b && std::isfinite(high) && high < 1e16) {
    high *= 2.0;
  }

  if (!std::isfinite(high) || high >= 1e16) {
    return mad_scale(resid);
  }

  double low = 0.0;
  for (int iter = 0; iter < max_iter; ++iter) {
    double mid = 0.5 * (low + high);
    if (!(mid > 0.0) || !std::isfinite(mid)) {
      mid = std::max(1e-16, 0.5 * high);
    }

    const double val = mean_rho(resid, mid, c);
    if (std::abs(val - b) < tol) {
      return mid;
    }
    if (val > b) {
      low = mid;
    } else {
      high = mid;
    }
  }

  return high;
}

ModelFrame build_model_frame(const std::string& formula,
                             const NamedColumns& data,
                             NAMethod na_method,
                             bool default_intercept) {
  ParsedFormula parsed = parse_formula(formula, default_intercept);
  if (data.empty()) {
    throw std::invalid_argument("build_model_frame: data is empty");
  }

  const vec y_full = lookup_column(data, parsed.response);
  const arma::uword n = y_full.n_elem;
  if (n == 0) {
    throw std::invalid_argument("build_model_frame: empty response");
  }

  for (const auto& kv : data) {
    if (kv.second.n_elem != n) {
      throw std::invalid_argument("build_model_frame: all columns must have same length");
    }
  }

  std::vector<vec> columns;
  columns.reserve(parsed.terms.size());
  for (const auto& term : parsed.terms) {
    columns.push_back(interaction_column(data, term));
  }

  arma::uvec keep(n, arma::fill::ones);
  for (arma::uword i = 0; i < n; ++i) {
    bool ok = std::isfinite(y_full[i]);
    for (const auto& col : columns) {
      ok = ok && std::isfinite(col[i]);
    }
    keep[i] = ok ? 1u : 0u;
  }

  const uvec kept_rows = arma::find(keep == 1u);
  if (kept_rows.n_elem == 0) {
    throw std::invalid_argument("build_model_frame: no complete cases remain");
  }
  if (na_method == NAMethod::kFail && kept_rows.n_elem != n) {
    throw std::invalid_argument("build_model_frame: NA/non-finite values present and na_method == kFail");
  }

  ModelFrame out;
  out.response_name = parsed.response;
  out.kept_rows = kept_rows;
  out.has_intercept = parsed.intercept;
  out.y = y_full.elem(kept_rows);

  const arma::uword p = static_cast<arma::uword>(columns.size() + (parsed.intercept ? 1 : 0));
  out.X.set_size(kept_rows.n_elem, p);

  arma::uword col_idx = 0;
  if (parsed.intercept) {
    out.X.col(col_idx).ones();
    out.coef_names.push_back("(Intercept)");
    ++col_idx;
  }

  for (std::size_t j = 0; j < columns.size(); ++j, ++col_idx) {
    out.X.col(col_idx) = columns[j].elem(kept_rows);
    out.coef_names.push_back(parsed.term_labels[j]);
  }

  return out;
}

SEstResult fit_s_estimator(const mat& X_in, const vec& y, const SEstControl& ctl) {
  if (X_in.n_rows != y.n_elem) {
    throw std::invalid_argument("fit_s_estimator: X rows must match y length");
  }
  if (X_in.n_rows == 0 || X_in.n_cols == 0) {
    throw std::invalid_argument("fit_s_estimator: X must be non-empty");
  }
  if (!all_finite(X_in) || !all_finite(y)) {
    throw std::invalid_argument("fit_s_estimator: X and y must be finite");
  }
  if (ctl.n_starts <= 0 || ctl.n_best_starts <= 0 || ctl.max_refine <= 0 || ctl.max_scale_iter <= 0) {
    throw std::invalid_argument("fit_s_estimator: iteration / start counts must be positive");
  }
  if (!(ctl.b > 0.0 && ctl.b < 1.0)) {
    throw std::invalid_argument("fit_s_estimator: b must be in (0, 1)");
  }
  if (!(ctl.c > 0.0)) {
    throw std::invalid_argument("fit_s_estimator: c must be > 0");
  }
  if (!(ctl.tol > 0.0) || !(ctl.scale_tol > 0.0)) {
    throw std::invalid_argument("fit_s_estimator: tolerances must be > 0");
  }
  if (ctl.use_fast_s && (ctl.fast_s_screen_subsets <= 0 || ctl.fast_s_screen_iters <= 0 || ctl.fast_s_keep <= 0)) {
    throw std::invalid_argument("fit_s_estimator: fast-s control values must be positive");
  }

  const mat X = ctl.add_intercept ? add_intercept(X_in) : X_in;
  const int n = static_cast<int>(X.n_rows);
  const int p = static_cast<int>(X.n_cols);

  if (n <= p) {
    throw std::invalid_argument("fit_s_estimator: need n > p");
  }
  if (!full_rank(X)) {
    throw std::invalid_argument("fit_s_estimator: X must have full column rank");
  }

  std::mt19937_64 rng(ctl.seed);
  std::vector<StartCandidate> candidates;
  candidates.reserve(static_cast<std::size_t>(ctl.n_starts + 1));

  int starts_tried = 0;
  int starts_used = 0;

  if (ctl.include_ols_start) {
    try {
      vec beta_ols;
      const bool ok = arma::solve(beta_ols, X, y);
      if (ok && beta_ols.is_finite()) {
        StartCandidate cand = refine_s_start(X, y, beta_ols, ctl,
                                             ctl.use_fast_s ? ctl.fast_s_screen_iters : ctl.max_refine);
        if (std::isfinite(cand.scale)) {
          candidates.push_back(cand);
          ++starts_used;
        }
      }
    } catch (...) {
    }
  }

  const int target_starts = ctl.use_fast_s ? ctl.fast_s_screen_subsets : ctl.n_starts;
  const int refine_screen_iters = ctl.use_fast_s ? ctl.fast_s_screen_iters : ctl.max_refine;
  while (starts_used < target_starts && starts_tried < target_starts * 25) {
    ++starts_tried;

    const std::vector<int> idx = sample_without_replacement(n, p, rng);
    mat X_sub(p, p);
    vec y_sub(p);

    for (int j = 0; j < p; ++j) {
      const arma::uword row = static_cast<arma::uword>(idx[static_cast<std::size_t>(j)]);
      X_sub.row(static_cast<arma::uword>(j)) = X.row(row);
      y_sub[static_cast<arma::uword>(j)] = y[row];
    }

    if (!full_rank(X_sub)) continue;

    try {
      const vec beta0 = solve_subset_ls(X_sub, y_sub);
      StartCandidate cand = refine_s_start(X, y, beta0, ctl, refine_screen_iters);
      if (std::isfinite(cand.scale)) {
        candidates.push_back(cand);
        ++starts_used;
      }
    } catch (...) {
      continue;
    }
  }

  if (candidates.empty()) {
    throw std::runtime_error("fit_s_estimator: no valid starts found");
  }

  std::sort(candidates.begin(), candidates.end(),
            [](const StartCandidate& a, const StartCandidate& b) {
              return a.scale < b.scale;
            });

  const int keep = ctl.use_fast_s
      ? std::min<int>(ctl.fast_s_keep, static_cast<int>(candidates.size()))
      : std::min<int>(ctl.n_best_starts, static_cast<int>(candidates.size()));

  StartCandidate best = candidates.front();
  for (int i = 0; i < keep; ++i) {
    try {
      StartCandidate cand = refine_s_start(X, y, candidates[static_cast<std::size_t>(i)].beta,
                                           ctl, ctl.max_refine);
      if (cand.scale < best.scale) {
        best = cand;
      }
    } catch (...) {
      continue;
    }
  }

  const vec fitted = X * best.beta;
  const vec resid = y - fitted;
  const double scale = solve_s_scale(resid, ctl.c, ctl.b, ctl.max_scale_iter, ctl.scale_tol);

  vec weights(resid.n_elem, arma::fill::ones);
  if (scale > 0.0) {
    for (arma::uword i = 0; i < resid.n_elem; ++i) {
      weights[i] = tukey_weight(resid[i] / scale, ctl.c);
    }
  }

  SEstResult out;
  out.coef = best.beta;
  out.fitted = fitted;
  out.resid = resid;
  out.weights = weights;
  out.scale = scale;
  out.objective = scale;
  out.converged = best.converged;
  out.iterations = best.iterations;
  out.starts_tried = starts_tried;
  out.starts_used = starts_used;
  out.status = best.converged ? SStatus::kOk : SStatus::kMaxRefineReached;
  out.message = best.converged ? "ok" : "maximum refinement iterations reached before strict convergence";
  out.X = X;
  out.y = y;
  return out;
}

SEstResult fit_s_estimator(const std::string& formula,
                           const NamedColumns& data,
                           const SEstControl& ctl,
                           NAMethod na_method) {
  ModelFrame mf = build_model_frame(formula, data, na_method, ctl.add_intercept);

  SEstControl inner_ctl = ctl;
  inner_ctl.add_intercept = false;

  SEstResult out = fit_s_estimator(mf.X, mf.y, inner_ctl);
  out.coef_names = mf.coef_names;
  out.kept_rows = mf.kept_rows;
  out.rows_dropped_na = static_cast<int>(data.begin()->second.n_elem - mf.kept_rows.n_elem);
  out.formula = formula;
  return out;
}

}  // namespace robust
