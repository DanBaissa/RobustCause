#include <Rcpp.h>

#include <string>
#include <vector>

#include "robust/robust_adstock.hpp"

using namespace Rcpp;

namespace {

robust::AdstockIncrementMethod parse_increment_method(const std::string& method) {
  if (method == "plain") return robust::AdstockIncrementMethod::kPlain;
  if (method == "huber") return robust::AdstockIncrementMethod::kHuber;
  if (method == "tanh") return robust::AdstockIncrementMethod::kTanh;
  if (method == "softsign") return robust::AdstockIncrementMethod::kSoftsign;
  if (method == "adaptive_clip") return robust::AdstockIncrementMethod::kAdaptiveClip;
  throw std::invalid_argument("Unknown adstock increment method: " + method);
}

std::vector<std::vector<double>> matrix_to_nested(const NumericMatrix& x) {
  std::vector<std::vector<double>> out(x.nrow(), std::vector<double>(x.ncol()));
  for (int i = 0; i < x.nrow(); ++i) {
    for (int j = 0; j < x.ncol(); ++j) {
      out[i][j] = x(i, j);
    }
  }
  return out;
}

}  // namespace

extern "C" SEXP rc_r_build_adstock(SEXP signalSEXP,
                                   SEXP unitIdsSEXP,
                                   SEXP xPrecleanSEXP,
                                   SEXP controlSEXP) {
  BEGIN_RCPP

  NumericVector signal_in(signalSEXP);
  IntegerVector unit_ids_in(unitIdsSEXP);
  NumericMatrix x_preclean_in(xPrecleanSEXP);
  List control(controlSEXP);

  std::vector<double> signal(signal_in.begin(), signal_in.end());
  std::vector<int> unit_ids(unit_ids_in.begin(), unit_ids_in.end());

  robust::RobustAdstockConfig cfg;
  cfg.preclean.enabled = as<bool>(control["preclean_enabled"]);
  cfg.preclean.nonnegative = as<bool>(control["nonnegative"]);
  cfg.preclean.use_mm = as<bool>(control["use_mm"]);
  cfg.preclean.residual_clip_multiplier = as<double>(control["residual_clip_multiplier"]);
  cfg.preclean.s_starts = as<int>(control["s_starts"]);
  cfg.preclean.s_max_iter = as<int>(control["s_max_iter"]);
  cfg.preclean.mm_max_iter = as<int>(control["mm_max_iter"]);
  cfg.preclean.mm_tukey_c = as<double>(control["mm_tukey_c"]);
  cfg.preclean.seed = static_cast<unsigned int>(as<int>(control["seed"]));
  cfg.adstock.rho = as<double>(control["rho"]);
  cfg.adstock.increment_method = parse_increment_method(as<std::string>(control["increment_method"]));
  cfg.adstock.clip_c = as<double>(control["clip_c"]);
  cfg.adstock.adaptive_upper = as<double>(control["adaptive_upper"]);

  std::vector<double> cleaned = signal;
  if (cfg.preclean.enabled) {
    const std::vector<std::vector<double>> x_preclean = matrix_to_nested(x_preclean_in);
    cleaned = cfg.preclean.use_mm
      ? robust::preclean_signal_mm(signal, x_preclean, cfg.preclean)
      : robust::preclean_signal_s(signal, x_preclean, cfg.preclean);
  }

  const std::vector<double> stock = robust::build_adstock(cleaned, unit_ids, cfg.adstock);

  return List::create(
    _["stock"] = wrap(stock),
    _["cleaned_signal"] = wrap(cleaned),
    _["increment_method"] = robust::increment_method_name(cfg.adstock.increment_method),
    _["preclean_enabled"] = cfg.preclean.enabled,
    _["preclean_method"] = cfg.preclean.enabled ? (cfg.preclean.use_mm ? "mm" : "s") : "none",
    _["rho"] = cfg.adstock.rho,
    _["clip_c"] = cfg.adstock.clip_c,
    _["adaptive_upper"] = cfg.adstock.adaptive_upper
  );

  END_RCPP
}
