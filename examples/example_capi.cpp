#include "robustcause/capi.h"

#include <iostream>
#include <vector>

int main() {
  const std::size_t n = 6;
  const std::size_t p = 2;

  // Column-major layout to match Armadillo and the C ABI.
  const std::vector<double> x = {
    -2.0, -1.0, 0.0, 1.0, 2.0, 3.0,
     1.0,  0.5, 1.0, 2.0, 1.5, 2.5
  };
  const std::vector<double> y = {0.1, 0.9, 1.8, 3.2, 3.9, 5.1};

  rc_rlm_options_t options = rc_default_rlm_options();
  std::vector<double> coef(p + 1);
  std::vector<double> fitted(n);
  std::vector<double> resid(n);
  std::vector<double> weights(n);
  std::vector<double> hat(n);
  char error[256];
  double scale = 0.0;
  int converged = 0;
  int iterations = 0;

  const rc_status_t status = rc_fit_rlm(
    x.data(),
    n,
    p,
    y.data(),
    &options,
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

  if (status != RC_STATUS_OK) {
    std::cerr << "rc_fit_rlm failed: " << error << "\n";
    return 1;
  }

  std::cout << "converged: " << converged << "\n";
  std::cout << "iterations: " << iterations << "\n";
  std::cout << "scale: " << scale << "\n";
  std::cout << "coef:";
  for (double value : coef) {
    std::cout << " " << value;
  }
  std::cout << "\n";
  return 0;
}
