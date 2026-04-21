#ifndef ROBUST_S_ESTIMATOR_HPP
#define ROBUST_S_ESTIMATOR_HPP

#include <armadillo>

#include <cstdint>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

namespace robust {

using arma::mat;
using arma::uvec;
using arma::vec;
using NamedColumns = std::unordered_map<std::string, vec>;

enum class SStatus {
  kOk = 0,
  kMaxRefineReached,
  kNoValidStarts,
  kInvalidArgument,
  kNumericalFailure
};

enum class NAMethod {
  kFail = 0,
  kOmit
};

struct SEstControl {
  bool add_intercept = true;
  int n_starts = 500;
  int n_best_starts = 25;
  int max_refine = 100;
  int max_scale_iter = 100;
  double tol = 1e-8;
  double scale_tol = 1e-10;
  double c = 1.54764;
  double b = 0.5;
  double ridge = 1e-10;
  double min_weight = 1e-12;
  std::uint64_t seed = 123456789ULL;

  // Fast-S-style screening.
  // This is a practical concentration-step implementation for speed.
  // It is not intended to claim exact parity with robustbase internals.
  bool use_fast_s = true;
  bool include_ols_start = true;
  int fast_s_screen_subsets = 500;
  int fast_s_screen_iters = 2;
  int fast_s_keep = 25;
};

struct ModelFrame {
  mat X;
  vec y;
  std::vector<std::string> coef_names;
  uvec kept_rows;
  bool has_intercept = true;
  std::string response_name;
};

struct SEstResult {
  vec coef;
  vec fitted;
  vec resid;
  vec weights;
  double scale = std::numeric_limits<double>::quiet_NaN();
  double objective = std::numeric_limits<double>::quiet_NaN();
  bool converged = false;
  int iterations = 0;
  int starts_tried = 0;
  int starts_used = 0;
  int rows_dropped_na = 0;
  SStatus status = SStatus::kInvalidArgument;
  std::string message;
  mat X;
  vec y;
  std::vector<std::string> coef_names;
  uvec kept_rows;
  std::string formula;
};

SEstResult fit_s_estimator(const mat& X, const vec& y, const SEstControl& ctl = {});
SEstResult fit_s_estimator(const std::string& formula,
                           const NamedColumns& data,
                           const SEstControl& ctl = {},
                           NAMethod na_method = NAMethod::kOmit);

ModelFrame build_model_frame(const std::string& formula,
                             const NamedColumns& data,
                             NAMethod na_method = NAMethod::kOmit,
                             bool default_intercept = true);

mat add_intercept(const mat& X);
double tukey_rho_normalized(double u, double c);
double tukey_weight(double u, double c);
double solve_s_scale(const vec& resid, double c, double b, int max_iter, double tol);

}  // namespace robust

#endif
