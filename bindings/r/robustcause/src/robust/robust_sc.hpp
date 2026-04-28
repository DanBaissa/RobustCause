#ifndef ROBUST_SC_HPP
#define ROBUST_SC_HPP

#include <RcppArmadillo.h>

namespace robust {

struct ScControl {
  int maxit = 500;
  double tol = 1e-8;
  double ridge = 1e-10;
};

struct MmscControl {
  int startup_maxit = 75;
  int maxit = 100;
  int subproblem_maxit = 500;
  double tol = 1e-8;
  double subproblem_tol = 1e-8;
  double l1_smoothing = 1e-6;
  double tukey_c = 4.685;
  double min_scale = 1e-8;
  double ridge = 1e-10;
  double min_time_weight = 1e-8;
};

struct ScPanel {
  arma::mat outcomes;
  arma::uword treated_index = 0;
  arma::uword treatment_start = 0;
  arma::uvec donor_indices;
  arma::vec predictor_treated;
  arma::mat predictor_donors;
  arma::vec predictor_weights;
};

struct ScData {
  arma::vec treated_pre;
  arma::mat donors_pre;
  arma::vec treated_post;
  arma::mat donors_post;
  arma::vec predictor_treated;
  arma::mat predictor_donors;
  arma::vec predictor_weights;
  arma::uword treated_index = 0;
  arma::uvec donor_indices;
  arma::uvec pre_period_indices;
  arma::uvec post_period_indices;
};

struct ScResult {
  arma::vec weights;
  arma::vec synthetic_pre;
  arma::vec synthetic_post;
  arma::vec pre_residuals;
  arma::vec post_gaps;
  double objective_value = 0.0;
  double pre_rmspe = 0.0;
  double weight_herfindahl = 0.0;
  double max_weight = 0.0;
  double effective_donors = 0.0;
  bool converged = false;
  int iterations = 0;
};

struct MmscResult {
  arma::vec weights;
  arma::vec startup_weights;
  arma::vec synthetic_pre;
  arma::vec synthetic_post;
  arma::vec pre_residuals;
  arma::vec post_gaps;
  arma::vec robust_time_weights;
  arma::uvec downweighted_periods;
  double scale = 0.0;
  double objective_value = 0.0;
  double startup_objective = 0.0;
  double standard_rmspe = 0.0;
  double startup_rmspe = 0.0;
  double mm_rmspe = 0.0;
  double effective_periods = 0.0;
  double weight_herfindahl = 0.0;
  double max_weight = 0.0;
  double effective_donors = 0.0;
  bool converged = false;
  int iterations = 0;
  ScResult standard_fit;
};

struct ScPlaceboResult {
  arma::uvec unit_indices;
  arma::mat weight_matrix;
  arma::mat post_gap_matrix;
  arma::vec pre_rmspe;
  arma::vec post_rmspe;
  arma::vec rmspe_ratio;
};

arma::vec project_simplex(const arma::vec& v);
double mad_scale(const arma::vec& resid, double min_scale = 1e-8);
ScData prepare_sc_data(const arma::mat& outcomes,
                       arma::uword treated_index,
                       arma::uword treatment_start,
                       const arma::uvec& donor_indices = arma::uvec(),
                       const arma::vec& predictor_treated = arma::vec(),
                       const arma::mat& predictor_donors = arma::mat(),
                       const arma::vec& predictor_weights = arma::vec());
ScData prepare_sc_data(const ScPanel& panel);

ScResult fit_sc(const arma::vec& treated_pre,
                const arma::mat& donors_pre,
                const arma::vec& treated_post,
                const arma::mat& donors_post,
                const ScControl& control = {},
                const arma::vec& predictor_treated = arma::vec(),
                const arma::mat& predictor_donors = arma::mat(),
                const arma::vec& predictor_weights = arma::vec());
ScResult fit_sc(const ScData& data, const ScControl& control = {});

MmscResult fit_mm_sc(const arma::vec& treated_pre,
                     const arma::mat& donors_pre,
                     const arma::vec& treated_post,
                     const arma::mat& donors_post,
                     const MmscControl& control = {},
                     const arma::vec& predictor_treated = arma::vec(),
                     const arma::mat& predictor_donors = arma::mat(),
                     const arma::vec& predictor_weights = arma::vec());
MmscResult fit_mm_sc(const ScData& data, const MmscControl& control = {});

ScPlaceboResult fit_sc_placebos(const arma::mat& outcomes,
                                int treatment_start,
                                const ScControl& control = {});
ScPlaceboResult fit_mm_sc_placebos(const arma::mat& outcomes,
                                   int treatment_start,
                                   const MmscControl& control = {});

}  // namespace robust

#endif
