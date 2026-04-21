#include "robust/robust_mmm.hpp"

#include <stdexcept>

namespace robust {
namespace {

void validate_design_rows(const std::vector<std::vector<double>>& X, size_t n, const char* name) {
  if (X.empty()) {
    throw std::invalid_argument(std::string(name) + " must not be empty");
  }
  const size_t p = X.front().size();
  if (p == 0) {
    throw std::invalid_argument(std::string(name) + " must have at least one column");
  }
  if (X.size() != n) {
    throw std::invalid_argument(std::string(name) + " must have one row per observation");
  }
  for (const auto& row : X) {
    if (row.size() != p) {
      throw std::invalid_argument(std::string(name) + " rows must have equal length");
    }
  }
}

void validate_channel_matrix(const std::vector<std::vector<double>>& channels, size_t n) {
  if (channels.empty()) {
    throw std::invalid_argument("channel_signals must not be empty");
  }
  for (const auto& ch : channels) {
    if (ch.size() != n) {
      throw std::invalid_argument("all channel signals must match the outcome length");
    }
  }
}

}  // namespace

double hill_transform(double x, double lambda) {
  if (lambda <= 0.0) {
    return x;
  }
  return x / (1.0 + lambda * x);
}

std::string mmm_fit_method_name(MmmFitMethod method) {
  switch (method) {
    case MmmFitMethod::kOLS:
      return "ols";
    case MmmFitMethod::kHuber:
      return "huber";
    case MmmFitMethod::kS:
      return "s";
    case MmmFitMethod::kMM:
      return "mm";
  }
  return "unknown";
}

std::vector<std::vector<double>> build_mmm_design(
  const std::vector<std::vector<double>>& transformed_channels,
  const std::vector<std::vector<double>>& controls
) {
  if (transformed_channels.empty()) {
    throw std::invalid_argument("transformed_channels must not be empty");
  }
  const size_t n = transformed_channels.front().size();
  validate_channel_matrix(transformed_channels, n);
  validate_design_rows(controls, n, "controls");

  const size_t k = transformed_channels.size();
  const size_t c = controls.front().size();
  std::vector<std::vector<double>> X(n, std::vector<double>(1 + k + c, 1.0));
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < k; ++j) {
      X[i][1 + j] = transformed_channels[j][i];
    }
    for (size_t q = 0; q < c; ++q) {
      X[i][1 + k + q] = controls[i][q];
    }
  }
  return X;
}

MmmFitResult fit_mmm(const std::vector<double>& y,
                     const std::vector<int>& unit_ids,
                     const std::vector<std::vector<double>>& channel_signals,
                     const std::vector<std::vector<double>>& preclean_X,
                     const std::vector<std::vector<double>>& controls,
                     const MmmConfig& cfg) {
  if (y.empty()) {
    throw std::invalid_argument("y must not be empty");
  }
  const size_t n = y.size();
  if (unit_ids.size() != n) {
    throw std::invalid_argument("unit_ids must match the outcome length");
  }
  validate_channel_matrix(channel_signals, n);
  validate_design_rows(preclean_X, n, "preclean_X");
  validate_design_rows(controls, n, "controls");
  if (cfg.channels.size() != channel_signals.size()) {
    throw std::invalid_argument("cfg.channels must match the number of channels");
  }

  const size_t k = channel_signals.size();
  MmmFitResult out;
  out.cleaned_signals.resize(k);
  out.stocks.resize(k);
  out.transformed_channels.resize(k);

  for (size_t j = 0; j < k; ++j) {
    const MmmChannelConfig& ch_cfg = cfg.channels[j];
    std::vector<double> signal = channel_signals[j];
    if (ch_cfg.preclean.enabled) {
      signal = ch_cfg.preclean.use_mm
        ? preclean_signal_mm(signal, preclean_X, ch_cfg.preclean)
        : preclean_signal_s(signal, preclean_X, ch_cfg.preclean);
    }
    out.cleaned_signals[j] = signal;
    out.stocks[j] = build_adstock(signal, unit_ids, ch_cfg.adstock);
    out.transformed_channels[j].resize(n, 0.0);
    for (size_t i = 0; i < n; ++i) {
      out.transformed_channels[j][i] =
        ch_cfg.apply_hill ? hill_transform(out.stocks[j][i], ch_cfg.hill_lambda) : out.stocks[j][i];
    }
  }

  out.design_matrix = build_mmm_design(out.transformed_channels, controls);

  switch (cfg.fit_method) {
    case MmmFitMethod::kOLS:
      out.regression = adstock_ols(out.design_matrix, y);
      break;
    case MmmFitMethod::kHuber:
      out.regression = adstock_huber_regression(out.design_matrix, y, cfg.huber_c, cfg.huber_max_iter);
      break;
    case MmmFitMethod::kS:
      out.regression = adstock_s_regression(out.design_matrix, y, cfg.s_starts, cfg.s_max_iter, cfg.seed);
      break;
    case MmmFitMethod::kMM:
      out.regression = adstock_mm_regression_from_s(
        out.design_matrix, y, cfg.s_starts, cfg.s_max_iter, cfg.mm_tukey_c, cfg.mm_max_iter, cfg.seed);
      break;
  }

  out.ok = out.regression.ok;
  out.full_coef = out.regression.coef;
  if (out.ok && out.full_coef.size() >= 1 + k) {
    out.channel_coef.assign(out.full_coef.begin() + 1, out.full_coef.begin() + 1 + k);
  }
  return out;
}

}  // namespace robust
