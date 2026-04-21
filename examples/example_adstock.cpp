#include "robust/robust_adstock.hpp"

#include <iomanip>
#include <iostream>
#include <vector>

int main() {
  const std::vector<double> signal{10.0, 12.0, 55.0, 9.0, 8.0, 45.0};
  const std::vector<int> unit_ids{1, 1, 1, 2, 2, 2};
  const std::vector<std::vector<double>> X_preclean{
    {1.0, 0.0},
    {1.0, 1.0},
    {1.0, 2.0},
    {1.0, 0.0},
    {1.0, 1.0},
    {1.0, 2.0},
  };

  robust::RobustAdstockConfig cfg;
  cfg.preclean.enabled = true;
  cfg.preclean.use_mm = true;
  cfg.preclean.nonnegative = true;
  cfg.adstock.rho = 0.8;
  cfg.adstock.increment_method = robust::AdstockIncrementMethod::kHuber;
  cfg.adstock.clip_c = 15.0;

  const std::vector<double> stock = robust::build_robust_adstock(signal, unit_ids, X_preclean, cfg);

  std::cout << "Robust adstock:\n";
  for (double value : stock) {
    std::cout << std::fixed << std::setprecision(3) << "  " << value << "\n";
  }

  return 0;
}
