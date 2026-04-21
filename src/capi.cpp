#include "robustcause/capi.h"

#include "robust/robust_rlm_hc.hpp"
#include "robust/robust_s_estimator.hpp"

#include <armadillo>

#include <algorithm>
#include <cstring>
#include <stdexcept>
#include <string>

namespace {

void set_error_message(char* dst, size_t dst_size, const std::string& message) {
  if (dst == nullptr || dst_size == 0) {
    return;
  }

  const size_t count = std::min(dst_size - 1, message.size());
  std::memcpy(dst, message.data(), count);
  dst[count] = '\0';
}

template <typename T>
void require_not_null(const T* ptr, const char* name) {
  if (ptr == nullptr) {
    throw std::invalid_argument(std::string(name) + " must not be null");
  }
}

robust::SEstControl to_cpp(const rc_s_options_t& options);

robust::PsiType to_cpp(rc_psi_type_t psi) {
  switch (psi) {
    case RC_PSI_HUBER:
      return robust::PsiType::kHuber;
    case RC_PSI_TUKEY_BISQUARE:
      return robust::PsiType::kTukeyBisquare;
    default:
      throw std::invalid_argument("unknown psi type");
  }
}

robust::RlmMethod to_cpp(rc_rlm_method_t method) {
  switch (method) {
    case RC_RLM_M:
      return robust::RlmMethod::kM;
    case RC_RLM_MM:
      return robust::RlmMethod::kMM;
    default:
      throw std::invalid_argument("unknown RLM method");
  }
}

robust::HCType to_cpp(rc_hc_type_t hc) {
  switch (hc) {
    case RC_HC0:
      return robust::HCType::kHC0;
    case RC_HC1:
      return robust::HCType::kHC1;
    case RC_HC2:
      return robust::HCType::kHC2;
    case RC_HC3:
      return robust::HCType::kHC3;
    case RC_HC4:
      return robust::HCType::kHC4;
    case RC_HC4M:
      return robust::HCType::kHC4m;
    case RC_HC5:
      return robust::HCType::kHC5;
    default:
      throw std::invalid_argument("unknown HC type");
  }
}

robust::RlmControl to_cpp(const rc_rlm_options_t& options) {
  robust::RlmControl ctl;
  ctl.method = to_cpp(options.method);
  ctl.psi = to_cpp(options.psi);
  ctl.tuning = options.tuning;
  ctl.maxit = options.maxit;
  ctl.tol = options.tol;
  ctl.add_intercept = options.add_intercept != 0;
  ctl.ridge = options.ridge;
  ctl.min_weight = options.min_weight;
  ctl.mm_s_control = to_cpp(options.mm_s_options);
  return ctl;
}

robust::SEstControl to_cpp(const rc_s_options_t& options) {
  robust::SEstControl ctl;
  ctl.add_intercept = options.add_intercept != 0;
  ctl.n_starts = options.n_starts;
  ctl.n_best_starts = options.n_best_starts;
  ctl.max_refine = options.max_refine;
  ctl.max_scale_iter = options.max_scale_iter;
  ctl.tol = options.tol;
  ctl.scale_tol = options.scale_tol;
  ctl.c = options.c;
  ctl.b = options.b;
  ctl.ridge = options.ridge;
  ctl.min_weight = options.min_weight;
  ctl.seed = options.seed;
  ctl.use_fast_s = options.use_fast_s != 0;
  ctl.include_ols_start = options.include_ols_start != 0;
  ctl.fast_s_screen_subsets = options.fast_s_screen_subsets;
  ctl.fast_s_screen_iters = options.fast_s_screen_iters;
  ctl.fast_s_keep = options.fast_s_keep;
  return ctl;
}

void copy_vec(const arma::vec& src, double* dst, const char* name) {
  require_not_null(dst, name);
  std::memcpy(dst, src.memptr(), sizeof(double) * static_cast<size_t>(src.n_elem));
}

void copy_mat(const arma::mat& src, double* dst, const char* name) {
  require_not_null(dst, name);
  std::memcpy(dst, src.memptr(), sizeof(double) * static_cast<size_t>(src.n_elem));
}

arma::mat view_matrix(const double* x, size_t n_rows, size_t n_cols) {
  require_not_null(x, "x");
  if (n_rows == 0 || n_cols == 0) {
    throw std::invalid_argument("x must have at least one row and one column");
  }
  return arma::mat(const_cast<double*>(x), n_rows, n_cols, false, true);
}

arma::vec view_vector(const double* x, size_t n) {
  require_not_null(x, "y");
  if (n == 0) {
    throw std::invalid_argument("y must be non-empty");
  }
  return arma::vec(const_cast<double*>(x), n, false, true);
}

template <typename Fn>
rc_status_t wrap_call(Fn&& fn, char* error_message, size_t error_message_size) {
  try {
    fn();
    set_error_message(error_message, error_message_size, "");
    return RC_STATUS_OK;
  } catch (const std::invalid_argument& ex) {
    set_error_message(error_message, error_message_size, ex.what());
    return RC_STATUS_INVALID_ARGUMENT;
  } catch (const std::exception& ex) {
    set_error_message(error_message, error_message_size, ex.what());
    return RC_STATUS_RUNTIME_ERROR;
  } catch (...) {
    set_error_message(error_message, error_message_size, "unknown error");
    return RC_STATUS_RUNTIME_ERROR;
  }
}

}  // namespace

extern "C" {

rc_rlm_options_t rc_default_rlm_options(void) {
  rc_rlm_options_t options;
  options.add_intercept = 1;
  options.maxit = 100;
  options.tol = 1e-8;
  options.tuning = 1.345;
  options.ridge = 1e-10;
  options.min_weight = 1e-12;
  options.psi = RC_PSI_HUBER;
  options.method = RC_RLM_M;
  options.mm_s_options = rc_default_s_options();
  return options;
}

rc_s_options_t rc_default_s_options(void) {
  rc_s_options_t options;
  options.add_intercept = 1;
  options.n_starts = 500;
  options.n_best_starts = 25;
  options.max_refine = 100;
  options.max_scale_iter = 100;
  options.use_fast_s = 1;
  options.include_ols_start = 1;
  options.fast_s_screen_subsets = 500;
  options.fast_s_screen_iters = 2;
  options.fast_s_keep = 25;
  options.tol = 1e-8;
  options.scale_tol = 1e-10;
  options.c = 1.54764;
  options.b = 0.5;
  options.ridge = 1e-10;
  options.min_weight = 1e-12;
  options.seed = 123456789ULL;
  return options;
}

rc_status_t rc_fit_rlm(
  const double* x,
  size_t n_rows,
  size_t n_cols,
  const double* y,
  const rc_rlm_options_t* options,
  double* coef_out,
  double* fitted_out,
  double* resid_out,
  double* weights_out,
  double* hat_out,
  double* scale_out,
  int* converged_out,
  int* iterations_out,
  char* error_message,
  size_t error_message_size
) {
  return wrap_call([&]() {
    require_not_null(options, "options");
    require_not_null(scale_out, "scale_out");
    require_not_null(converged_out, "converged_out");
    require_not_null(iterations_out, "iterations_out");

    const arma::mat X = view_matrix(x, n_rows, n_cols);
    const arma::vec Y = view_vector(y, n_rows);
    const robust::RlmResult fit = robust::fit_rlm(X, Y, to_cpp(*options));

    copy_vec(fit.coef, coef_out, "coef_out");
    copy_vec(fit.fitted, fitted_out, "fitted_out");
    copy_vec(fit.resid, resid_out, "resid_out");
    copy_vec(fit.weights, weights_out, "weights_out");
    copy_vec(fit.hat, hat_out, "hat_out");
    *scale_out = fit.scale;
    *converged_out = fit.converged ? 1 : 0;
    *iterations_out = fit.iterations;
  }, error_message, error_message_size);
}

rc_status_t rc_fit_s_estimator(
  const double* x,
  size_t n_rows,
  size_t n_cols,
  const double* y,
  const rc_s_options_t* options,
  double* coef_out,
  double* fitted_out,
  double* resid_out,
  double* weights_out,
  double* scale_out,
  int* converged_out,
  int* iterations_out,
  char* error_message,
  size_t error_message_size
) {
  return wrap_call([&]() {
    require_not_null(options, "options");
    require_not_null(scale_out, "scale_out");
    require_not_null(converged_out, "converged_out");
    require_not_null(iterations_out, "iterations_out");

    const arma::mat X = view_matrix(x, n_rows, n_cols);
    const arma::vec Y = view_vector(y, n_rows);
    const robust::SEstResult fit = robust::fit_s_estimator(X, Y, to_cpp(*options));

    copy_vec(fit.coef, coef_out, "coef_out");
    copy_vec(fit.fitted, fitted_out, "fitted_out");
    copy_vec(fit.resid, resid_out, "resid_out");
    copy_vec(fit.weights, weights_out, "weights_out");
    *scale_out = fit.scale;
    *converged_out = fit.converged ? 1 : 0;
    *iterations_out = fit.iterations;
  }, error_message, error_message_size);
}

rc_status_t rc_rlm_confint_normal(
  const double* x,
  size_t n_rows,
  size_t n_cols,
  const double* y,
  const rc_rlm_options_t* fit_options,
  rc_hc_type_t hc_type,
  double zcrit,
  double* vcov_out,
  double* se_out,
  double* ci_out,
  char* error_message,
  size_t error_message_size
) {
  return wrap_call([&]() {
    require_not_null(fit_options, "fit_options");
    const arma::mat X = view_matrix(x, n_rows, n_cols);
    const arma::vec Y = view_vector(y, n_rows);
    const robust::RlmResult fit = robust::fit_rlm(X, Y, to_cpp(*fit_options));
    const robust::InferenceResult inf = robust::confint_normal(fit, to_cpp(hc_type), zcrit, fit_options->ridge);

    copy_mat(inf.vcov, vcov_out, "vcov_out");
    copy_vec(inf.se, se_out, "se_out");
    copy_mat(inf.ci, ci_out, "ci_out");
  }, error_message, error_message_size);
}

}  // extern "C"
