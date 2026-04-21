#ifndef ROBUST_RLM_HC_HPP
#define ROBUST_RLM_HC_HPP

#include <armadillo>

#include <limits>
#include <string>

namespace robust {

using arma::mat;
using arma::vec;

enum class PsiType {
  kHuber = 0,
  kTukeyBisquare
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

struct RlmControl {
  PsiType psi = PsiType::kHuber;
  double tuning = 1.345;
  int maxit = 100;
  double tol = 1e-8;
  bool add_intercept = true;
  double ridge = 1e-10;
  double min_weight = 1e-12;
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
  mat X;
  vec y;
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
std::string hc_name(HCType type);

}  // namespace robust

#endif
