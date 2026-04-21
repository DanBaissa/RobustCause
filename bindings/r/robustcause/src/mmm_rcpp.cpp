#include <Rcpp.h>

#include <string>
#include <vector>

#include "robust/robust_mmm.hpp"

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

robust::MmmFitMethod parse_fit_method(const std::string& method) {
  if (method == "ols") return robust::MmmFitMethod::kOLS;
  if (method == "huber") return robust::MmmFitMethod::kHuber;
  if (method == "s") return robust::MmmFitMethod::kS;
  if (method == "mm") return robust::MmmFitMethod::kMM;
  throw std::invalid_argument("Unknown MMM fit method: " + method);
}

std::vector<std::vector<double>> matrix_to_rows(const NumericMatrix& x) {
  std::vector<std::vector<double>> out(x.nrow(), std::vector<double>(x.ncol()));
  for (int i = 0; i < x.nrow(); ++i) {
    for (int j = 0; j < x.ncol(); ++j) {
      out[i][j] = x(i, j);
    }
  }
  return out;
}

std::vector<std::vector<double>> matrix_to_columns(const NumericMatrix& x) {
  std::vector<std::vector<double>> out(x.ncol(), std::vector<double>(x.nrow()));
  for (int j = 0; j < x.ncol(); ++j) {
    for (int i = 0; i < x.nrow(); ++i) {
      out[j][i] = x(i, j);
    }
  }
  return out;
}

NumericMatrix columns_to_matrix(const std::vector<std::vector<double>>& cols) {
  if (cols.empty()) return NumericMatrix(0, 0);
  NumericMatrix out(cols.front().size(), cols.size());
  for (size_t j = 0; j < cols.size(); ++j) {
    for (size_t i = 0; i < cols[j].size(); ++i) {
      out(i, j) = cols[j][i];
    }
  }
  return out;
}

NumericMatrix rows_to_matrix(const std::vector<std::vector<double>>& rows) {
  if (rows.empty()) return NumericMatrix(0, 0);
  NumericMatrix out(rows.size(), rows.front().size());
  for (size_t i = 0; i < rows.size(); ++i) {
    for (size_t j = 0; j < rows[i].size(); ++j) {
      out(i, j) = rows[i][j];
    }
  }
  return out;
}

}  // namespace

extern "C" SEXP rc_r_fit_mmm(SEXP ySEXP,
                             SEXP unitIdsSEXP,
                             SEXP channelsSEXP,
                             SEXP precleanXSEXP,
                             SEXP controlsSEXP,
                             SEXP controlSEXP) {
  BEGIN_RCPP

  NumericVector y_in(ySEXP);
  IntegerVector unit_ids_in(unitIdsSEXP);
  NumericMatrix channels_in(channelsSEXP);
  NumericMatrix preclean_x_in(precleanXSEXP);
  NumericMatrix controls_in(controlsSEXP);
  List control(controlSEXP);

  std::vector<double> y(y_in.begin(), y_in.end());
  std::vector<int> unit_ids(unit_ids_in.begin(), unit_ids_in.end());

  robust::MmmConfig cfg;
  cfg.fit_method = parse_fit_method(as<std::string>(control["fit_method"]));
  cfg.huber_c = as<double>(control["huber_c"]);
  cfg.huber_max_iter = as<int>(control["huber_max_iter"]);
  cfg.s_starts = as<int>(control["s_starts"]);
  cfg.s_max_iter = as<int>(control["s_max_iter"]);
  cfg.mm_max_iter = as<int>(control["mm_max_iter"]);
  cfg.mm_tukey_c = as<double>(control["mm_tukey_c"]);
  cfg.seed = static_cast<unsigned int>(as<int>(control["seed"]));

  CharacterVector channel_names = control["channel_names"];
  LogicalVector preclean_enabled = control["preclean_enabled"];
  CharacterVector preclean_method = control["preclean_method"];
  LogicalVector nonnegative = control["nonnegative"];
  NumericVector residual_clip_multiplier = control["residual_clip_multiplier"];
  IntegerVector preclean_seed = control["preclean_seed"];
  NumericVector rho = control["rho"];
  CharacterVector increment_method = control["increment_method"];
  NumericVector clip_c = control["clip_c"];
  NumericVector adaptive_upper = control["adaptive_upper"];
  NumericVector hill_lambda = control["hill_lambda"];
  LogicalVector apply_hill = control["apply_hill"];

  cfg.channels.resize(channels_in.ncol());
  for (int j = 0; j < channels_in.ncol(); ++j) {
    cfg.channels[j].name = as<std::string>(channel_names[j]);
    cfg.channels[j].preclean.enabled = preclean_enabled[j];
    cfg.channels[j].preclean.use_mm = as<std::string>(preclean_method[j]) == "mm";
    cfg.channels[j].preclean.nonnegative = nonnegative[j];
    cfg.channels[j].preclean.residual_clip_multiplier = residual_clip_multiplier[j];
    cfg.channels[j].preclean.s_starts = cfg.s_starts;
    cfg.channels[j].preclean.s_max_iter = cfg.s_max_iter;
    cfg.channels[j].preclean.mm_max_iter = cfg.mm_max_iter;
    cfg.channels[j].preclean.mm_tukey_c = cfg.mm_tukey_c;
    cfg.channels[j].preclean.seed = static_cast<unsigned int>(preclean_seed[j]);
    cfg.channels[j].adstock.rho = rho[j];
    cfg.channels[j].adstock.increment_method = parse_increment_method(as<std::string>(increment_method[j]));
    cfg.channels[j].adstock.clip_c = clip_c[j];
    cfg.channels[j].adstock.adaptive_upper = adaptive_upper[j];
    cfg.channels[j].hill_lambda = hill_lambda[j];
    cfg.channels[j].apply_hill = apply_hill[j];
  }

  const robust::MmmFitResult fit = robust::fit_mmm(
    y,
    unit_ids,
    matrix_to_columns(channels_in),
    matrix_to_rows(preclean_x_in),
    matrix_to_rows(controls_in),
    cfg
  );

  return List::create(
    _["ok"] = fit.ok,
    _["full_coef"] = wrap(fit.full_coef),
    _["channel_coef"] = wrap(fit.channel_coef),
    _["cleaned_signals"] = columns_to_matrix(fit.cleaned_signals),
    _["stocks"] = columns_to_matrix(fit.stocks),
    _["transformed_channels"] = columns_to_matrix(fit.transformed_channels),
    _["design_matrix"] = rows_to_matrix(fit.design_matrix),
    _["regression_coef"] = wrap(fit.regression.coef),
    _["regression_ok"] = fit.regression.ok,
    _["regression_scale"] = fit.regression.scale,
    _["fit_method"] = mmm_fit_method_name(cfg.fit_method)
  );

  END_RCPP
}
