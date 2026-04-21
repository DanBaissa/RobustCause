#include "robust/robust_mmm.hpp"

#include <cmath>
#include <iostream>
#include <random>
#include <vector>

int main() {
  const int n_units = 80;
  const int n_periods = 8;
  const double rho_true = 0.8;
  const std::vector<double> beta_true{1.0, 0.8, 0.6};
  const std::vector<double> lambda_true{0.15, 0.10, 0.12};

  std::mt19937 rng(2468);
  std::normal_distribution<double> norm01(0.0, 1.0);
  std::student_t_distribution<double> t3(3.0);
  std::uniform_real_distribution<double> unif01(0.0, 1.0);

  std::vector<int> unit_ids;
  std::vector<double> y;
  std::vector<std::vector<double>> channels(3);
  std::vector<std::vector<double>> preclean_X;
  std::vector<std::vector<double>> controls;

  std::vector<double> q(n_units);
  for (int i = 0; i < n_units; ++i) {
    q[i] = norm01(rng);
  }

  for (int i = 0; i < n_units; ++i) {
    std::vector<double> stock_prev(3, 0.0);
    for (int t = 0; t < n_periods; ++t) {
      const double sin_t = std::sin(static_cast<double>(t) / 2.5);
      const double x1 = q[i] + 0.3 * norm01(rng);
      const double x2 = 0.4 * sin_t + 0.3 * norm01(rng);

      std::vector<double> clean_signal(3, 0.0);
      clean_signal[0] = std::exp(0.5 + 0.45 * x1 + 0.15 * sin_t + 0.20 * norm01(rng));
      clean_signal[1] = std::exp(0.4 + 0.30 * x2 - 0.10 * sin_t + 0.20 * norm01(rng));
      clean_signal[2] = std::exp(0.45 + 0.20 * x1 + 0.20 * x2 + 0.20 * norm01(rng));

      std::vector<double> observed_signal = clean_signal;
      for (int j = 0; j < 3; ++j) {
        if (unif01(rng) < 0.025) {
          observed_signal[j] *= (8.0 + 25.0 * unif01(rng));
        }
        if (unif01(rng) < 0.020) {
          observed_signal[j] = 0.0;
        }
        stock_prev[j] = rho_true * stock_prev[j] + clean_signal[j];
      }

      double outcome = 0.0;
      for (int j = 0; j < 3; ++j) {
        outcome += beta_true[j] * robust::hill_transform(stock_prev[j], lambda_true[j]);
      }
      outcome += 0.5 * x1 - 0.25 * x2 + t3(rng);

      unit_ids.push_back(i);
      y.push_back(outcome);
      for (int j = 0; j < 3; ++j) {
        channels[j].push_back(observed_signal[j]);
      }
      preclean_X.push_back({1.0, x1, x2, sin_t});
      controls.push_back({x1, x2});
    }
  }

  robust::MmmConfig cfg;
  cfg.fit_method = robust::MmmFitMethod::kHuber;
  cfg.channels.resize(3);
  for (int j = 0; j < 3; ++j) {
    cfg.channels[j].name = "channel_" + std::to_string(j + 1);
    cfg.channels[j].preclean.enabled = true;
    cfg.channels[j].preclean.use_mm = true;
    cfg.channels[j].preclean.nonnegative = true;
    cfg.channels[j].preclean.seed = 1000 + j;
    cfg.channels[j].adstock.rho = rho_true;
    cfg.channels[j].adstock.increment_method = robust::AdstockIncrementMethod::kAdaptiveClip;
    cfg.channels[j].hill_lambda = lambda_true[j];
  }

  const robust::MmmFitResult fit = robust::fit_mmm(y, unit_ids, channels, preclean_X, controls, cfg);

  std::cout << "MMM fit ok: " << (fit.ok ? "true" : "false") << "\n";
  for (size_t j = 0; j < fit.channel_coef.size(); ++j) {
    std::cout << "channel_" << (j + 1) << ": " << fit.channel_coef[j] << "\n";
  }
  return 0;
}
