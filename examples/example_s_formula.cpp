#include "robust/robust_s_estimator.hpp"

#include <armadillo>
#include <iostream>
#include <limits>
#include <random>

int main() {
  using robust::NAMethod;
  using robust::NamedColumns;
  using robust::SEstControl;
  using robust::SEstResult;
  using robust::fit_s_estimator;

  constexpr int n = 300;
  arma::vec x1(n), x2(n), x3(n), y(n);

  std::mt19937_64 rng(42);
  std::normal_distribution<double> norm(0.0, 1.0);
  std::uniform_real_distribution<double> uni(0.0, 1.0);

  for (int i = 0; i < n; ++i) {
    x1[i] = norm(rng);
    x2[i] = norm(rng);
    x3[i] = uni(rng);
    y[i] = 1.0 + 2.0 * x1[i] - 1.5 * x2[i] + 0.75 * x3[i] + 1.25 * x1[i] * x2[i] + 0.25 * norm(rng);
  }

  for (int i : {5, 15, 60, 125, 200}) {
    y[i] += (i % 2 == 0 ? 12.0 : -10.0);
  }

  x2[7] = std::numeric_limits<double>::quiet_NaN();
  x3[9] = std::numeric_limits<double>::quiet_NaN();

  NamedColumns data;
  data["y"] = y;
  data["x1"] = x1;
  data["x2"] = x2;
  data["x3"] = x3;

  SEstControl ctl;
  ctl.use_fast_s = true;
  ctl.fast_s_screen_subsets = 600;
  ctl.fast_s_screen_iters = 2;
  ctl.fast_s_keep = 30;
  ctl.max_refine = 100;

  SEstResult fit = fit_s_estimator("y ~ x1 + x2 + x3 + x1:x2", data, ctl, NAMethod::kOmit);

  std::cout << "formula: " << fit.formula << "\n";
  std::cout << "rows dropped for NA: " << fit.rows_dropped_na << "\n";
  std::cout << "status: " << fit.message << "\n";
  std::cout << "scale: " << fit.scale << "\n";
  std::cout << "coef:\n";
  for (arma::uword j = 0; j < fit.coef.n_elem; ++j) {
    const std::string name = j < fit.coef_names.size() ? fit.coef_names[j] : ("b" + std::to_string(j));
    std::cout << "  " << name << ": " << fit.coef[j] << "\n";
  }

  return 0;
}
