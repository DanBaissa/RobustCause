#include "robust/robust_sc.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace robust {
namespace {

struct QpControl {
  int maxit = 500;
  double tol = 1e-8;
  double ridge = 1e-10;
};

void validate_predictor_inputs(const arma::vec& predictor_treated,
                               const arma::mat& predictor_donors,
                               const arma::vec& predictor_weights,
                               arma::uword donor_count) {
  const bool has_treated = predictor_treated.n_elem > 0;
  const bool has_donors = predictor_donors.n_elem > 0;
  const bool has_weights = predictor_weights.n_elem > 0;
  if (!(has_treated || has_donors || has_weights)) {
    return;
  }
  if (!has_treated || !has_donors) {
    throw std::invalid_argument("predictor_treated and predictor_donors must both be supplied when predictors are used");
  }
  if (predictor_donors.n_rows != predictor_treated.n_elem || predictor_donors.n_cols != donor_count) {
    throw std::invalid_argument("predictor_donors must align with predictor_treated and donor count");
  }
  if (has_weights && predictor_weights.n_elem != predictor_treated.n_elem) {
    throw std::invalid_argument("predictor_weights length must match predictor_treated");
  }
  if (!predictor_treated.is_finite() || !predictor_donors.is_finite() || !predictor_weights.is_finite()) {
    throw std::invalid_argument("predictor inputs must be finite");
  }
}

void validate_sc_inputs(const arma::vec& treated_pre,
                        const arma::mat& donors_pre,
                        const arma::vec& treated_post,
                        const arma::mat& donors_post,
                        const arma::vec& predictor_treated,
                        const arma::mat& predictor_donors,
                        const arma::vec& predictor_weights) {
  if (treated_pre.n_elem == 0) {
    throw std::invalid_argument("treated_pre must be non-empty");
  }
  if (!treated_pre.is_finite() || !donors_pre.is_finite() ||
      !treated_post.is_finite() || !donors_post.is_finite()) {
    throw std::invalid_argument("synthetic-control inputs must be finite");
  }
  if (donors_pre.n_rows != treated_pre.n_elem || donors_pre.n_cols == 0) {
    throw std::invalid_argument("donors_pre must have one row per pre-treatment period and at least one donor");
  }
  if (treated_post.n_elem > 0 && donors_post.n_rows != treated_post.n_elem) {
    throw std::invalid_argument("donors_post must have one row per post-treatment period");
  }
  if (donors_post.n_elem > 0 && donors_post.n_cols != donors_pre.n_cols) {
    throw std::invalid_argument("donors_post must have the same number of donor columns as donors_pre");
  }
  validate_predictor_inputs(
    predictor_treated, predictor_donors, predictor_weights, donors_pre.n_cols);
}

template <class Control>
void validate_controls(const Control& control) {
  if (control.maxit <= 0) {
    throw std::invalid_argument("iteration counts must be positive");
  }
  if (control.tol <= 0.0 || control.ridge < 0.0) {
    throw std::invalid_argument("tolerances must be positive and ridge must be nonnegative");
  }
}

void validate_controls(const MmscControl& control) {
  if (control.startup_maxit <= 0 || control.maxit <= 0 || control.subproblem_maxit <= 0) {
    throw std::invalid_argument("iteration counts must be positive");
  }
  if (control.tukey_c <= 0.0 || control.tol <= 0.0 || control.subproblem_tol <= 0.0 ||
      control.min_scale <= 0.0 || control.l1_smoothing <= 0.0 || control.min_time_weight <= 0.0 ||
      control.ridge < 0.0) {
    throw std::invalid_argument("invalid MM-SC tuning parameters");
  }
}

double mean_square(const arma::vec& x) {
  if (x.n_elem == 0) {
    return 0.0;
  }
  return arma::dot(x, x) / static_cast<double>(x.n_elem);
}

double root_mean_square(const arma::vec& x) {
  return std::sqrt(mean_square(x));
}

double smooth_l1_objective(const arma::vec& resid, double eps) {
  return arma::accu(arma::sqrt(arma::square(resid) + eps * eps));
}

double tukey_rho_scalar(double u, double c) {
  const double au = std::abs(u);
  if (au >= c) {
    return (c * c) / 6.0;
  }
  const double x = u / c;
  const double one_minus = 1.0 - x * x;
  return (c * c / 6.0) * (1.0 - one_minus * one_minus * one_minus);
}

arma::vec tukey_weights(const arma::vec& u, double c, double min_weight) {
  arma::vec out(u.n_elem, arma::fill::ones);
  for (arma::uword i = 0; i < u.n_elem; ++i) {
    const double au = std::abs(u[i]);
    if (au >= c) {
      out[i] = min_weight;
      continue;
    }
    const double x = u[i] / c;
    const double one_minus = 1.0 - x * x;
    out[i] = std::max(one_minus * one_minus, min_weight);
  }
  return out;
}

double tukey_objective(const arma::vec& resid, double scale, double c) {
  double out = 0.0;
  for (arma::uword i = 0; i < resid.n_elem; ++i) {
    out += tukey_rho_scalar(resid[i] / scale, c);
  }
  return out;
}

QpControl to_qp_control(const ScControl& control) {
  QpControl out;
  out.maxit = control.maxit;
  out.tol = control.tol;
  out.ridge = control.ridge;
  return out;
}

QpControl to_qp_control(const MmscControl& control) {
  QpControl out;
  out.maxit = control.subproblem_maxit;
  out.tol = control.subproblem_tol;
  out.ridge = control.ridge;
  return out;
}

struct AugmentedDesign {
  arma::vec y;
  arma::mat X;
  arma::vec obs_weights;
  arma::uword pre_rows = 0;
};

AugmentedDesign make_design(const arma::vec& treated_pre,
                            const arma::mat& donors_pre,
                            const arma::vec& predictor_treated,
                            const arma::mat& predictor_donors,
                            const arma::vec& predictor_weights,
                            const arma::vec* time_weights = nullptr) {
  AugmentedDesign out;
  out.pre_rows = treated_pre.n_elem;

  const arma::uword feature_rows = predictor_treated.n_elem;
  out.y.set_size(treated_pre.n_elem + feature_rows);
  out.X.set_size(donors_pre.n_rows + feature_rows, donors_pre.n_cols);
  out.y.head(treated_pre.n_elem) = treated_pre;
  out.X.rows(0, donors_pre.n_rows - 1) = donors_pre;

  out.obs_weights = arma::vec(out.y.n_elem, arma::fill::ones);
  if (time_weights != nullptr) {
    out.obs_weights.head(treated_pre.n_elem) = *time_weights;
  }

  if (feature_rows > 0) {
    arma::vec weights = predictor_weights.n_elem > 0
      ? predictor_weights
      : arma::vec(feature_rows, arma::fill::ones);
    for (arma::uword i = 0; i < feature_rows; ++i) {
      const arma::uword row = treated_pre.n_elem + i;
      out.y[row] = predictor_treated[i];
      out.X.row(row) = predictor_donors.row(i);
      out.obs_weights[row] = weights[i];
    }
  }

  return out;
}

arma::vec solve_weighted_simplex_qp(const arma::mat& X,
                                    const arma::vec& y,
                                    const arma::vec& obs_weights,
                                    const arma::vec& start,
                                    const QpControl& control) {
  const arma::vec sqrt_w = arma::sqrt(obs_weights);
  arma::mat Xw = X.each_col() % sqrt_w;
  arma::vec yw = y % sqrt_w;
  arma::mat H = 2.0 * (Xw.t() * Xw);
  H.diag() += 2.0 * control.ridge;
  arma::vec g = -2.0 * (Xw.t() * yw);

  double lipschitz = arma::norm(H, 2);
  if (!std::isfinite(lipschitz) || lipschitz <= 0.0) {
    lipschitz = 1.0;
  }
  const double step = 1.0 / lipschitz;

  arma::vec w = project_simplex(start);
  for (int iter = 0; iter < control.maxit; ++iter) {
    const arma::vec grad = H * w + g;
    const arma::vec w_new = project_simplex(w - step * grad);
    if (arma::norm(w_new - w, 2) <= control.tol * std::max(1.0, arma::norm(w, 2))) {
      return w_new;
    }
    w = w_new;
  }
  return w;
}

arma::uvec downweighted_periods(const arma::vec& time_weights) {
  std::vector<arma::uword> idx;
  for (arma::uword i = 0; i < time_weights.n_elem; ++i) {
    if (time_weights[i] < 0.5) {
      idx.push_back(i);
    }
  }
  arma::uvec out(idx.size());
  for (arma::uword i = 0; i < out.n_elem; ++i) {
    out[i] = idx[i];
  }
  return out;
}

double donor_effective_size(const arma::vec& weights) {
  const double hhi = arma::dot(weights, weights);
  if (hhi <= 0.0) {
    return 0.0;
  }
  return 1.0 / hhi;
}

ScResult finish_sc_result(const arma::vec& treated_pre,
                          const arma::mat& donors_pre,
                          const arma::vec& treated_post,
                          const arma::mat& donors_post,
                          const arma::vec& weights,
                          double objective_value,
                          bool converged,
                          int iterations) {
  ScResult out;
  out.weights = weights;
  out.synthetic_pre = donors_pre * weights;
  out.synthetic_post = donors_post.n_elem > 0 ? donors_post * weights : arma::vec();
  out.pre_residuals = treated_pre - out.synthetic_pre;
  out.post_gaps = treated_post.n_elem > 0 ? treated_post - out.synthetic_post : arma::vec();
  out.objective_value = objective_value;
  out.pre_rmspe = root_mean_square(out.pre_residuals);
  out.weight_herfindahl = arma::dot(weights, weights);
  out.max_weight = arma::max(weights);
  out.effective_donors = donor_effective_size(weights);
  out.converged = converged;
  out.iterations = iterations;
  return out;
}

arma::vec startup_l1_fit(const arma::vec& treated_pre,
                         const arma::mat& donors_pre,
                         const arma::vec& predictor_treated,
                         const arma::mat& predictor_donors,
                         const arma::vec& predictor_weights,
                         const MmscControl& control) {
  arma::vec w(donors_pre.n_cols, arma::fill::ones);
  w /= static_cast<double>(donors_pre.n_cols);

  for (int iter = 0; iter < control.startup_maxit; ++iter) {
    const arma::vec resid = treated_pre - donors_pre * w;
    arma::vec time_weights = 1.0 / arma::sqrt(arma::square(resid) + control.l1_smoothing * control.l1_smoothing);
    time_weights = arma::clamp(time_weights, control.min_time_weight, std::numeric_limits<double>::infinity());
    const AugmentedDesign design = make_design(
      treated_pre, donors_pre, predictor_treated, predictor_donors, predictor_weights, &time_weights);
    const arma::vec w_new = solve_weighted_simplex_qp(
      design.X, design.y, design.obs_weights, w, to_qp_control(control));
    if (arma::norm(w_new - w, 2) <= control.tol * std::max(1.0, arma::norm(w, 2))) {
      return w_new;
    }
    w = w_new;
  }
  return w;
}

ScPlaceboResult finish_placebo_result(const arma::uvec& unit_indices,
                                      const arma::mat& weight_matrix,
                                      const arma::mat& post_gap_matrix,
                                      const arma::vec& pre_rmspe) {
  ScPlaceboResult out;
  out.unit_indices = unit_indices;
  out.weight_matrix = weight_matrix;
  out.post_gap_matrix = post_gap_matrix;
  out.pre_rmspe = pre_rmspe;
  out.post_rmspe.set_size(unit_indices.n_elem);
  out.rmspe_ratio.set_size(unit_indices.n_elem);
  for (arma::uword i = 0; i < unit_indices.n_elem; ++i) {
    out.post_rmspe[i] = post_gap_matrix.n_cols > 0 ? root_mean_square(post_gap_matrix.row(i).t()) : 0.0;
    out.rmspe_ratio[i] = out.pre_rmspe[i] > 0.0 ? out.post_rmspe[i] / out.pre_rmspe[i] : 0.0;
  }
  return out;
}

}  // namespace

arma::vec project_simplex(const arma::vec& v) {
  if (v.n_elem == 0) {
    return arma::vec();
  }

  arma::vec u = arma::sort(v, "descend");
  double cssv = 0.0;
  arma::uword rho = 0;
  bool found = false;
  for (arma::uword i = 0; i < u.n_elem; ++i) {
    cssv += u[i];
    const double t = (cssv - 1.0) / static_cast<double>(i + 1);
    if (u[i] - t > 0.0) {
      rho = i;
      found = true;
    }
  }
  if (!found) {
    return arma::vec(v.n_elem, arma::fill::value(1.0 / static_cast<double>(v.n_elem)));
  }
  const double theta = (arma::accu(u.head(rho + 1)) - 1.0) / static_cast<double>(rho + 1);
  arma::vec w = v - theta;
  w.transform([](double x) { return x > 0.0 ? x : 0.0; });
  const double sum_w = arma::accu(w);
  if (sum_w <= 0.0 || !std::isfinite(sum_w)) {
    return arma::vec(v.n_elem, arma::fill::value(1.0 / static_cast<double>(v.n_elem)));
  }
  return w / sum_w;
}

double mad_scale(const arma::vec& resid, double min_scale) {
  if (resid.n_elem == 0) {
    return min_scale;
  }
  const double med = arma::median(resid);
  const double mad = 1.4826 * arma::median(arma::abs(resid - med));
  if (!std::isfinite(mad) || mad < min_scale) {
    return std::max(root_mean_square(resid), min_scale);
  }
  return mad;
}

ScData prepare_sc_data(const arma::mat& outcomes,
                       arma::uword treated_index,
                       arma::uword treatment_start,
                       const arma::uvec& donor_indices,
                       const arma::vec& predictor_treated,
                       const arma::mat& predictor_donors,
                       const arma::vec& predictor_weights) {
  if (outcomes.n_rows == 0 || outcomes.n_cols < 2) {
    throw std::invalid_argument("outcomes must have at least one period and at least two units");
  }
  if (!outcomes.is_finite()) {
    throw std::invalid_argument("outcomes must be finite");
  }
  if (treated_index >= outcomes.n_cols) {
    throw std::invalid_argument("treated_index out of range");
  }
  if (treatment_start == 0 || treatment_start >= outcomes.n_rows) {
    throw std::invalid_argument("treatment_start must leave at least one pre- and one post-treatment period");
  }

  arma::uvec donors = donor_indices;
  if (donors.n_elem == 0) {
    donors.set_size(outcomes.n_cols - 1);
    arma::uword k = 0;
    for (arma::uword j = 0; j < outcomes.n_cols; ++j) {
      if (j != treated_index) {
        donors[k++] = j;
      }
    }
  }
  if (arma::any(donors == treated_index)) {
    throw std::invalid_argument("donor_indices must exclude the treated unit");
  }
  for (arma::uword i = 0; i < donors.n_elem; ++i) {
    if (donors[i] >= outcomes.n_cols) {
      throw std::invalid_argument("donor_indices out of range");
    }
  }

  validate_predictor_inputs(predictor_treated, predictor_donors, predictor_weights, donors.n_elem);

  ScData out;
  out.treated_index = treated_index;
  out.donor_indices = donors;
  out.pre_period_indices = arma::regspace<arma::uvec>(0, treatment_start - 1);
  out.post_period_indices = arma::regspace<arma::uvec>(treatment_start, outcomes.n_rows - 1);
  out.treated_pre = outcomes.col(treated_index).head(treatment_start);
  out.treated_post = outcomes.col(treated_index).tail(outcomes.n_rows - treatment_start);
  const arma::mat pre_block = outcomes.rows(0, treatment_start - 1);
  const arma::mat post_block = outcomes.rows(treatment_start, outcomes.n_rows - 1);
  out.donors_pre = pre_block.cols(donors);
  out.donors_post = post_block.cols(donors);
  out.predictor_treated = predictor_treated;
  out.predictor_donors = predictor_donors;
  out.predictor_weights = predictor_weights;
  return out;
}

ScData prepare_sc_data(const ScPanel& panel) {
  return prepare_sc_data(
    panel.outcomes,
    panel.treated_index,
    panel.treatment_start,
    panel.donor_indices,
    panel.predictor_treated,
    panel.predictor_donors,
    panel.predictor_weights);
}

ScResult fit_sc(const arma::vec& treated_pre,
                const arma::mat& donors_pre,
                const arma::vec& treated_post,
                const arma::mat& donors_post,
                const ScControl& control,
                const arma::vec& predictor_treated,
                const arma::mat& predictor_donors,
                const arma::vec& predictor_weights) {
  validate_sc_inputs(
    treated_pre, donors_pre, treated_post, donors_post, predictor_treated, predictor_donors, predictor_weights);
  validate_controls(control);

  arma::vec start(donors_pre.n_cols, arma::fill::ones);
  start /= static_cast<double>(donors_pre.n_cols);
  const AugmentedDesign design = make_design(
    treated_pre, donors_pre, predictor_treated, predictor_donors, predictor_weights);
  const arma::vec weights = solve_weighted_simplex_qp(
    design.X, design.y, design.obs_weights, start, to_qp_control(control));
  const double objective = arma::dot(design.obs_weights, arma::square(design.y - design.X * weights));
  return finish_sc_result(treated_pre, donors_pre, treated_post, donors_post, weights, objective, true, control.maxit);
}

ScResult fit_sc(const ScData& data, const ScControl& control) {
  return fit_sc(
    data.treated_pre,
    data.donors_pre,
    data.treated_post,
    data.donors_post,
    control,
    data.predictor_treated,
    data.predictor_donors,
    data.predictor_weights);
}

MmscResult fit_mm_sc(const arma::vec& treated_pre,
                     const arma::mat& donors_pre,
                     const arma::vec& treated_post,
                     const arma::mat& donors_post,
                     const MmscControl& control,
                     const arma::vec& predictor_treated,
                     const arma::mat& predictor_donors,
                     const arma::vec& predictor_weights) {
  validate_sc_inputs(
    treated_pre, donors_pre, treated_post, donors_post, predictor_treated, predictor_donors, predictor_weights);
  validate_controls(control);

  MmscResult out;
  out.standard_fit = fit_sc(
    treated_pre, donors_pre, treated_post, donors_post, ScControl{control.subproblem_maxit, control.subproblem_tol, control.ridge},
    predictor_treated, predictor_donors, predictor_weights);
  out.standard_rmspe = out.standard_fit.pre_rmspe;

  out.startup_weights = startup_l1_fit(
    treated_pre, donors_pre, predictor_treated, predictor_donors, predictor_weights, control);
  out.startup_objective = smooth_l1_objective(treated_pre - donors_pre * out.startup_weights, control.l1_smoothing);
  out.scale = mad_scale(treated_pre - donors_pre * out.startup_weights, control.min_scale);

  if (out.scale <= control.min_scale) {
    out.weights = out.startup_weights;
    out.synthetic_pre = donors_pre * out.weights;
    out.synthetic_post = donors_post.n_elem > 0 ? donors_post * out.weights : arma::vec();
    out.pre_residuals = treated_pre - out.synthetic_pre;
    out.post_gaps = treated_post.n_elem > 0 ? treated_post - out.synthetic_post : arma::vec();
    out.robust_time_weights = arma::vec(treated_pre.n_elem, arma::fill::ones);
    out.objective_value = tukey_objective(out.pre_residuals, control.min_scale, control.tukey_c);
    out.startup_rmspe = root_mean_square(treated_pre - donors_pre * out.startup_weights);
    out.mm_rmspe = out.startup_rmspe;
    out.effective_periods = static_cast<double>(treated_pre.n_elem);
    out.weight_herfindahl = arma::dot(out.weights, out.weights);
    out.max_weight = arma::max(out.weights);
    out.effective_donors = donor_effective_size(out.weights);
    out.downweighted_periods = arma::uvec();
    out.converged = true;
    return out;
  }

  arma::vec w = out.startup_weights;
  arma::vec time_weights(treated_pre.n_elem, arma::fill::ones);
  double prev_obj = tukey_objective(treated_pre - donors_pre * w, out.scale, control.tukey_c);

  for (int iter = 0; iter < control.maxit; ++iter) {
    const arma::vec resid = treated_pre - donors_pre * w;
    time_weights = tukey_weights(resid / out.scale, control.tukey_c, control.min_time_weight);
    const AugmentedDesign design = make_design(
      treated_pre, donors_pre, predictor_treated, predictor_donors, predictor_weights, &time_weights);
    const arma::vec w_new = solve_weighted_simplex_qp(
      design.X, design.y, design.obs_weights, w, to_qp_control(control));
    const double obj = tukey_objective(treated_pre - donors_pre * w_new, out.scale, control.tukey_c);
    out.iterations = iter + 1;
    if (arma::norm(w_new - w, 2) <= control.tol * std::max(1.0, arma::norm(w, 2)) &&
        std::abs(obj - prev_obj) <= control.tol * std::max(1.0, std::abs(prev_obj))) {
      w = w_new;
      prev_obj = obj;
      out.converged = true;
      break;
    }
    w = w_new;
    prev_obj = obj;
  }

  out.weights = w;
  out.synthetic_pre = donors_pre * out.weights;
  out.synthetic_post = donors_post.n_elem > 0 ? donors_post * out.weights : arma::vec();
  out.pre_residuals = treated_pre - out.synthetic_pre;
  out.post_gaps = treated_post.n_elem > 0 ? treated_post - out.synthetic_post : arma::vec();
  out.robust_time_weights = tukey_weights(out.pre_residuals / out.scale, control.tukey_c, control.min_time_weight);
  out.downweighted_periods = downweighted_periods(out.robust_time_weights);
  out.objective_value = tukey_objective(out.pre_residuals, out.scale, control.tukey_c);
  out.startup_rmspe = root_mean_square(treated_pre - donors_pre * out.startup_weights);
  out.mm_rmspe = root_mean_square(out.pre_residuals);
  out.effective_periods = std::pow(arma::accu(out.robust_time_weights), 2.0) /
    std::max(arma::accu(arma::square(out.robust_time_weights)), control.min_time_weight);
  out.weight_herfindahl = arma::dot(out.weights, out.weights);
  out.max_weight = arma::max(out.weights);
  out.effective_donors = donor_effective_size(out.weights);
  return out;
}

MmscResult fit_mm_sc(const ScData& data, const MmscControl& control) {
  return fit_mm_sc(
    data.treated_pre,
    data.donors_pre,
    data.treated_post,
    data.donors_post,
    control,
    data.predictor_treated,
    data.predictor_donors,
    data.predictor_weights);
}

ScPlaceboResult fit_sc_placebos(const arma::mat& outcomes,
                                int treatment_start,
                                const ScControl& control) {
  const arma::uword treatment_start_uword = static_cast<arma::uword>(treatment_start);
  if (outcomes.n_cols < 2) {
    throw std::invalid_argument("placebo analysis requires at least two units");
  }
  arma::uvec unit_indices = arma::regspace<arma::uvec>(0, outcomes.n_cols - 1);
  arma::mat weight_matrix(unit_indices.n_elem, outcomes.n_cols - 1, arma::fill::zeros);
  arma::mat post_gap_matrix(unit_indices.n_elem, outcomes.n_rows - treatment_start_uword, arma::fill::zeros);
  arma::vec pre_rmspe(unit_indices.n_elem, arma::fill::zeros);

  for (arma::uword i = 0; i < unit_indices.n_elem; ++i) {
    const ScData data = prepare_sc_data(outcomes, i, treatment_start_uword);
    const ScResult fit = fit_sc(data, control);
    weight_matrix.row(i) = fit.weights.t();
    if (fit.post_gaps.n_elem > 0) {
      post_gap_matrix.row(i) = fit.post_gaps.t();
    }
    pre_rmspe[i] = fit.pre_rmspe;
  }
  return finish_placebo_result(unit_indices, weight_matrix, post_gap_matrix, pre_rmspe);
}

ScPlaceboResult fit_mm_sc_placebos(const arma::mat& outcomes,
                                   int treatment_start,
                                   const MmscControl& control) {
  const arma::uword treatment_start_uword = static_cast<arma::uword>(treatment_start);
  if (outcomes.n_cols < 2) {
    throw std::invalid_argument("placebo analysis requires at least two units");
  }
  arma::uvec unit_indices = arma::regspace<arma::uvec>(0, outcomes.n_cols - 1);
  arma::mat weight_matrix(unit_indices.n_elem, outcomes.n_cols - 1, arma::fill::zeros);
  arma::mat post_gap_matrix(unit_indices.n_elem, outcomes.n_rows - treatment_start_uword, arma::fill::zeros);
  arma::vec pre_rmspe(unit_indices.n_elem, arma::fill::zeros);

  for (arma::uword i = 0; i < unit_indices.n_elem; ++i) {
    const ScData data = prepare_sc_data(outcomes, i, treatment_start_uword);
    const MmscResult fit = fit_mm_sc(data, control);
    weight_matrix.row(i) = fit.weights.t();
    if (fit.post_gaps.n_elem > 0) {
      post_gap_matrix.row(i) = fit.post_gaps.t();
    }
    pre_rmspe[i] = fit.mm_rmspe;
  }
  return finish_placebo_result(unit_indices, weight_matrix, post_gap_matrix, pre_rmspe);
}

}  // namespace robust
