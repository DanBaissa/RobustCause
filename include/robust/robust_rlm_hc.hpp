#ifndef ROBUST_RLM_HC_HPP
#define ROBUST_RLM_HC_HPP

#include "robust/robust_s_estimator.hpp"

#include <armadillo>

#include <limits>
#include <string>

namespace robust {

using arma::mat;
using arma::vec;

enum class PsiType {
  kHuber = 0,
  kTukeyBisquare,
  kHampel
};

enum class RlmMethod {
  kM = 0,
  kMM
};

enum class HCType {
  kHC0 = 0,
  kHC1,
  kHC2,
  kHC3,
  kHC4,
  kHC4m,
  kHC5
};

enum class WeightMethod {
  kInvVar = 0,
  kCase
};

enum class InitMethod {
  kOLS = 0,
  kLTS,
  kUser
};

enum class ScaleEstimate {
  kMAD = 0,
  kHuber,
  kProposal2
};

enum class TestVector {
  kCoef = 0,
  kResid,
  kWeight
};

struct RlmControl {
  RlmMethod method = RlmMethod::kM;
  PsiType psi = PsiType::kHuber;
  double tuning = 1.345;
  int maxit = 100;
  double tol = 1e-8;
  double acc = 1e-8;
  bool add_intercept = true;
  double ridge = 1e-10;
  double min_weight = 1e-12;
  vec prior_weights;
  WeightMethod wt_method = WeightMethod::kCase;
  InitMethod init_method = InitMethod::kOLS;
  vec init_coef;
  ScaleEstimate scale_est = ScaleEstimate::kMAD;
  double k2 = 1.345;
  TestVector test_vec = TestVector::kCoef;
  SEstControl mm_s_control = {};
};

struct RlmResult {
  vec coef;
  vec fitted;
  vec resid;
  vec weights;
  vec hat;
  double scale = std::numeric_limits<double>::quiet_NaN();
  bool converged = false;
  int iterations = 0;
  RlmMethod method = RlmMethod::kM;
  WeightMethod wt_method = WeightMethod::kCase;
  vec prior_weights;
  mat X;
  vec y;
  mat X_hc;
  vec resid_hc;
  vec weights_hc;
};

struct InferenceResult {
  mat vcov;
  vec se;
  mat ci;  // columns: lower, upper
};

RlmResult fit_rlm(const mat& X, const vec& y, const RlmControl& ctl = {});
mat vcov_hc(const RlmResult& fit, HCType type = HCType::kHC3, double ridge = 1e-10);
InferenceResult confint_normal(const RlmResult& fit,
                               HCType type = HCType::kHC3,
                               double zcrit = 1.959963984540054,
                               double ridge = 1e-10);

mat add_intercept(const mat& X);
vec psi_weights(const vec& u, PsiType psi, double tuning);
std::string rlm_method_name(RlmMethod method);
std::string hc_name(HCType type);
std::string psi_name(PsiType psi);
std::string weight_method_name(WeightMethod method);
std::string init_method_name(InitMethod method);
std::string scale_estimate_name(ScaleEstimate method);
std::string test_vector_name(TestVector method);

}  // namespace robust

#endif
