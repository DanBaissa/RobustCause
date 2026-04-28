#include "robustcause/capi.h"
#include "robustcause/robustcause.hpp"

#include <armadillo>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

namespace {

bool approx(double lhs, double rhs, double tol) {
  return std::abs(lhs - rhs) <= tol;
}

void require(bool condition, const char* message) {
  if (!condition) {
    std::cerr << message << "\n";
    std::exit(1);
  }
}

}  // namespace

int main() {
  arma::mat X(8, 2);
  X.col(0) = arma::vec({-3, -2, -1, 0, 1, 2, 3, 4});
  X.col(1) = arma::vec({1, 1, 2, 2, 3, 3, 4, 4});

  arma::vec y = 0.75 + 1.5 * X.col(0) - 0.4 * X.col(1);
  y[2] += 3.5;
  y[6] -= 2.0;

  robust::RlmControl rlm_ctl;
  robust::RlmResult rlm_fit = robust::fit_rlm(X, y, rlm_ctl);
  require(rlm_fit.coef.n_elem == 3, "RLM coefficient size mismatch");
  require(rlm_fit.converged, "RLM failed to converge in smoke test");
  require(approx(rlm_fit.coef[1], 1.5, 0.5), "RLM slope estimate is not in the expected range");

  robust::SEstControl s_ctl;
  s_ctl.n_starts = 50;
  s_ctl.fast_s_screen_subsets = 50;
  s_ctl.fast_s_keep = 10;
  robust::SEstResult s_fit = robust::fit_s_estimator(X, y, s_ctl);
  require(s_fit.coef.n_elem == 3, "S-estimator coefficient size mismatch");
  require(s_fit.scale >= 0.0, "S-estimator returned an invalid scale");

  const std::vector<double> signal{10.0, 12.0, 100.0, 5.0, 7.0, 80.0};
  const std::vector<int> unit_ids{1, 1, 1, 2, 2, 2};
  const std::vector<std::vector<double>> X_preclean{
    {1.0, 0.0},
    {1.0, 1.0},
    {1.0, 2.0},
    {1.0, 0.0},
    {1.0, 1.0},
    {1.0, 2.0},
  };

  robust::RobustAdstockConfig adstock_ctl;
  adstock_ctl.preclean.enabled = true;
  adstock_ctl.preclean.use_mm = true;
  adstock_ctl.preclean.nonnegative = true;
  adstock_ctl.adstock.rho = 0.8;
  adstock_ctl.adstock.increment_method = robust::AdstockIncrementMethod::kHuber;
  adstock_ctl.adstock.clip_c = 20.0;

  const std::vector<double> robust_stock =
    robust::build_robust_adstock(signal, unit_ids, X_preclean, adstock_ctl);
  require(robust_stock.size() == signal.size(), "Robust adstock size mismatch");
  require(robust_stock[2] < 100.0 + 0.8 * robust_stock[1], "Robust adstock did not damp the outlier");

  robust::AdstockConfig plain_ctl;
  plain_ctl.rho = 0.8;
  const std::vector<double> plain_stock = robust::build_adstock(signal, unit_ids, plain_ctl);
  require(approx(plain_stock[1], 0.8 * plain_stock[0] + signal[1], 1e-9),
          "Plain adstock recurrence failed");
  require(approx(plain_stock[3], signal[3], 1e-9),
          "Plain adstock failed to reset for a new unit");

  const std::vector<double> y_mmm{5.5, 6.2, 7.0, 4.8, 5.0, 5.8};
  const std::vector<std::vector<double>> channels_mmm{
    {10.0, 12.0, 20.0, 8.0, 7.0, 15.0},
    {5.0, 6.0, 9.0, 4.0, 3.0, 7.0}
  };
  const std::vector<std::vector<double>> preclean_mmm{
    {1.0, 0.0},
    {1.0, 1.0},
    {1.0, 2.0},
    {1.0, 0.0},
    {1.0, 1.0},
    {1.0, 2.0}
  };
  const std::vector<std::vector<double>> controls_mmm{
    {0.2},
    {0.4},
    {0.6},
    {0.1},
    {0.3},
    {0.5}
  };
  robust::MmmConfig mmm_ctl;
  mmm_ctl.fit_method = robust::MmmFitMethod::kHuber;
  mmm_ctl.channels.resize(2);
  for (int j = 0; j < 2; ++j) {
    mmm_ctl.channels[j].preclean.enabled = true;
    mmm_ctl.channels[j].preclean.use_mm = true;
    mmm_ctl.channels[j].adstock.rho = 0.7;
    mmm_ctl.channels[j].adstock.increment_method = robust::AdstockIncrementMethod::kHuber;
    mmm_ctl.channels[j].adstock.clip_c = 10.0;
    mmm_ctl.channels[j].hill_lambda = 0.1;
  }
  const robust::MmmFitResult mmm_fit =
    robust::fit_mmm(y_mmm, unit_ids, channels_mmm, preclean_mmm, controls_mmm, mmm_ctl);
  require(mmm_fit.ok, "MMM fit failed");
  require(mmm_fit.channel_coef.size() == channels_mmm.size(), "MMM channel coefficient size mismatch");
  require(mmm_fit.design_matrix.size() == y_mmm.size(), "MMM design size mismatch");

  arma::mat outcomes_sc(7, 4);
  outcomes_sc.col(0) = arma::vec{1.0, 2.1, 9.0, 4.0, 5.1, 6.0, 7.0};
  outcomes_sc.col(1) = arma::vec{1.1, 2.0, 3.0, 4.1, 5.0, 6.1, 7.2};
  outcomes_sc.col(2) = arma::vec{0.8, 1.9, 3.2, 3.9, 5.2, 5.9, 6.9};
  outcomes_sc.col(3) = arma::vec{1.3, 2.2, 2.8, 4.2, 4.9, 6.0, 7.1};

  robust::ScPanel sc_panel;
  sc_panel.outcomes = outcomes_sc;
  sc_panel.treated_index = 0;
  sc_panel.treatment_start = 5;
  sc_panel.predictor_treated = arma::vec{2.0, 4.0};
  sc_panel.predictor_donors.set_size(2, 3);
  sc_panel.predictor_donors.row(0) = arma::rowvec{2.0, 1.9, 2.1};
  sc_panel.predictor_donors.row(1) = arma::rowvec{4.1, 3.9, 4.0};
  sc_panel.predictor_weights = arma::vec{5.0, 5.0};

  const robust::ScData sc_data = robust::prepare_sc_data(sc_panel);
  const robust::ScResult sc_standard_fit = robust::fit_sc(sc_data);
  require(sc_standard_fit.weights.n_elem == 3, "SC weight length mismatch");
  require(approx(arma::accu(sc_standard_fit.weights), 1.0, 1e-8), "SC weights do not sum to one");

  robust::MmscControl sc_ctl;
  const robust::MmscResult sc_fit = robust::fit_mm_sc(sc_data, sc_ctl);
  require(sc_fit.weights.n_elem == 3, "MM-SC weight length mismatch");
  require(approx(arma::accu(sc_fit.weights), 1.0, 1e-8), "MM-SC weights do not sum to one");
  require(arma::all(sc_fit.weights >= -1e-10), "MM-SC weights must be nonnegative");
  require(sc_fit.robust_time_weights[2] < 0.5, "MM-SC did not downweight the contaminated period");
  require(sc_fit.standard_fit.weights.n_elem == 3, "Embedded standard SC result missing");
  require(sc_fit.mm_rmspe <= sc_fit.startup_rmspe + 1e-8, "MM-SC refinement should not worsen startup RMSPE materially");

  const robust::ScPlaceboResult placebo_fit = robust::fit_mm_sc_placebos(outcomes_sc, 5, sc_ctl);
  require(placebo_fit.unit_indices.n_elem == 4, "Placebo unit count mismatch");
  require(placebo_fit.post_gap_matrix.n_cols == 2, "Placebo post gap width mismatch");
  require(placebo_fit.rmspe_ratio.n_elem == 4, "Placebo RMSPE ratios missing");

  rc_rlm_options_t capi_ctl = rc_default_rlm_options();
  std::vector<double> coef(3);
  std::vector<double> fitted(8);
  std::vector<double> resid(8);
  std::vector<double> weights(8);
  std::vector<double> hat(8);
  double scale = 0.0;
  int converged = 0;
  int iterations = 0;
  char error[256];

  const rc_status_t status = rc_fit_rlm(
    X.memptr(),
    X.n_rows,
    X.n_cols,
    y.memptr(),
    &capi_ctl,
    coef.data(),
    fitted.data(),
    resid.data(),
    weights.data(),
    hat.data(),
    &scale,
    &converged,
    &iterations,
    error,
    sizeof(error)
  );

  require(status == RC_STATUS_OK, error);
  require(converged == 1, "C API RLM fit did not converge");
  require(approx(coef[1], rlm_fit.coef[1], 1e-6), "C API and C++ coefficient mismatch");

  return 0;
}
