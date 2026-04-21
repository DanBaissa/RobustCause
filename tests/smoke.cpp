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
