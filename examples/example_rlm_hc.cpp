#include "robust/robust_rlm_hc.hpp"

#include <armadillo>
#include <iostream>
#include <random>

int main() {
  using robust::HCType;
  using robust::InferenceResult;
  using robust::RlmControl;
  using robust::RlmResult;
  using robust::fit_rlm;
  using robust::confint_normal;
  using robust::hc_name;

  constexpr int n = 400;
  arma::mat X(n, 3);
  arma::vec y(n);

  std::mt19937_64 rng(123);
  std::normal_distribution<double> norm(0.0, 1.0);

  for (int i = 0; i < n; ++i) {
    const double x1 = norm(rng);
    const double x2 = norm(rng);
    const double x3 = norm(rng);
    X(i, 0) = x1;
    X(i, 1) = x2;
    X(i, 2) = x3;
    y[i] = 1.0 + 2.0 * x1 - 1.5 * x2 + 0.8 * x3 + 0.35 * norm(rng);
  }

  for (int i : {10, 75, 160, 320}) {
    y[i] += (i % 2 == 0 ? 10.0 : -9.0);
  }

  RlmControl ctl;
  ctl.psi = robust::PsiType::kHuber;
  ctl.tuning = 1.345;

  RlmResult fit = fit_rlm(X, y, ctl);

  std::cout << "converged: " << (fit.converged ? "true" : "false") << "\n";
  std::cout << "iterations: " << fit.iterations << "\n";
  std::cout << "scale: " << fit.scale << "\n";
  std::cout << "coef:\n" << fit.coef << "\n";

  for (HCType hc : {HCType::kHC0, HCType::kHC1, HCType::kHC2, HCType::kHC3, HCType::kHC4, HCType::kHC4m, HCType::kHC5}) {
    InferenceResult inf = confint_normal(fit, hc);
    std::cout << "=== " << hc_name(hc) << " ===\n";
    std::cout << "se:\n" << inf.se << "\n";
    std::cout << "ci:\n" << inf.ci << "\n";
  }

  return 0;
}
