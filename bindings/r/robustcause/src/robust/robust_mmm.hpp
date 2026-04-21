#ifndef ROBUST_MMM_HPP
#define ROBUST_MMM_HPP

#include "robust_adstock.hpp"

#include <string>
#include <vector>

namespace robust {

enum class MmmFitMethod {
  kOLS = 0,
  kHuber,
  kS,
  kMM
};

struct MmmChannelConfig {
  std::string name;
  AdstockPrecleanConfig preclean;
  AdstockConfig adstock;
  double hill_lambda = 0.10;
  bool apply_hill = true;
};

struct MmmConfig {
  std::vector<MmmChannelConfig> channels;
  MmmFitMethod fit_method = MmmFitMethod::kHuber;
  double huber_c = 1.345;
  int huber_max_iter = 75;
  int s_starts = 20;
  int s_max_iter = 50;
  int mm_max_iter = 50;
  double mm_tukey_c = 4.685;
  unsigned int seed = 1234;
};

struct MmmFitResult {
  bool ok = false;
  std::vector<double> full_coef;
  std::vector<double> channel_coef;
  std::vector<std::vector<double>> cleaned_signals;
  std::vector<std::vector<double>> stocks;
  std::vector<std::vector<double>> transformed_channels;
  std::vector<std::vector<double>> design_matrix;
  AdstockRegressionResult regression;
};

double hill_transform(double x, double lambda);
std::string mmm_fit_method_name(MmmFitMethod method);
std::vector<std::vector<double>> build_mmm_design(
  const std::vector<std::vector<double>>& transformed_channels,
  const std::vector<std::vector<double>>& controls
);
MmmFitResult fit_mmm(const std::vector<double>& y,
                     const std::vector<int>& unit_ids,
                     const std::vector<std::vector<double>>& channel_signals,
                     const std::vector<std::vector<double>>& preclean_X,
                     const std::vector<std::vector<double>>& controls,
                     const MmmConfig& cfg = {});

}  // namespace robust

#endif
