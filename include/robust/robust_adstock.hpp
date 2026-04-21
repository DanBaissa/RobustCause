#ifndef ROBUST_ADSTOCK_HPP
#define ROBUST_ADSTOCK_HPP

#include <string>
#include <vector>

namespace robust {

enum class AdstockIncrementMethod {
  kPlain = 0,
  kHuber,
  kTanh,
  kSoftsign,
  kAdaptiveClip
};

struct AdstockRegressionResult {
  std::vector<double> coef;
  bool ok = false;
  double scale = 0.0;
};

struct AdstockPrecleanConfig {
  bool enabled = false;
  bool nonnegative = false;
  bool use_mm = true;
  double residual_clip_multiplier = 2.5;
  int s_starts = 20;
  int s_max_iter = 50;
  int mm_max_iter = 50;
  double mm_tukey_c = 4.685;
  unsigned int seed = 1234;
};

struct AdstockConfig {
  double rho = 0.85;
  AdstockIncrementMethod increment_method = AdstockIncrementMethod::kPlain;
  double clip_c = 1.5;
  double adaptive_upper = 0.0;
};

struct RobustAdstockConfig {
  AdstockPrecleanConfig preclean;
  AdstockConfig adstock;
};

double mean(const std::vector<double>& x);
double median(std::vector<double> x);
double mad(const std::vector<double>& x);
double correlation(const std::vector<double>& x, const std::vector<double>& y);
double log1p_safe(double x);
double huber_clip(double x, double c);
double tanh_bound(double x, double c);
double softsign_bound(double x, double c);
double adaptive_upper_from_signal(const std::vector<double>& signal, double k = 3.0);
std::string increment_method_name(AdstockIncrementMethod method);

AdstockRegressionResult adstock_ols(const std::vector<std::vector<double>>& X,
                                    const std::vector<double>& y);
AdstockRegressionResult adstock_huber_regression(const std::vector<std::vector<double>>& X,
                                                 const std::vector<double>& y,
                                                 double c = 1.345,
                                                 int max_iter = 75);
AdstockRegressionResult adstock_s_regression(const std::vector<std::vector<double>>& X,
                                             const std::vector<double>& y,
                                             int n_starts = 20,
                                             int max_iter = 50,
                                             unsigned int seed = 1234);
AdstockRegressionResult adstock_mm_regression_from_s(const std::vector<std::vector<double>>& X,
                                                     const std::vector<double>& y,
                                                     int s_starts = 20,
                                                     int s_max_iter = 50,
                                                     double c_m = 4.685,
                                                     int mm_max_iter = 50,
                                                     unsigned int seed = 1234);

std::vector<double> adstock_fitted_values(const std::vector<std::vector<double>>& X,
                                          const std::vector<double>& beta);
std::vector<double> preclean_signal_s(const std::vector<double>& signal,
                                      const std::vector<std::vector<double>>& X,
                                      const AdstockPrecleanConfig& cfg);
std::vector<double> preclean_signal_mm(const std::vector<double>& signal,
                                       const std::vector<std::vector<double>>& X,
                                       const AdstockPrecleanConfig& cfg);
std::vector<double> build_adstock(const std::vector<double>& signal,
                                  const std::vector<int>& unit_ids,
                                  const AdstockConfig& cfg = {});
std::vector<double> build_robust_adstock(const std::vector<double>& signal,
                                         const std::vector<int>& unit_ids,
                                         const std::vector<std::vector<double>>& X_preclean,
                                         const RobustAdstockConfig& cfg = {});

}  // namespace robust

#endif
