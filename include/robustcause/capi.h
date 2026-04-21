#ifndef ROBUSTCAUSE_CAPI_H
#define ROBUSTCAUSE_CAPI_H

#include <stddef.h>

#if defined(_WIN32) && defined(ROBUSTCAUSE_BUILD_SHARED)
#  ifdef ROBUSTCAUSE_EXPORTS
#    define ROBUSTCAUSE_API __declspec(dllexport)
#  else
#    define ROBUSTCAUSE_API __declspec(dllimport)
#  endif
#else
#  define ROBUSTCAUSE_API
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef enum rc_status_t {
  RC_STATUS_OK = 0,
  RC_STATUS_INVALID_ARGUMENT = 1,
  RC_STATUS_RUNTIME_ERROR = 2
} rc_status_t;

typedef enum rc_psi_type_t {
  RC_PSI_HUBER = 0,
  RC_PSI_TUKEY_BISQUARE = 1
} rc_psi_type_t;

typedef enum rc_hc_type_t {
  RC_HC0 = 0,
  RC_HC1 = 1,
  RC_HC2 = 2,
  RC_HC3 = 3,
  RC_HC4 = 4,
  RC_HC4M = 5,
  RC_HC5 = 6
} rc_hc_type_t;

typedef struct rc_rlm_options_t {
  int add_intercept;
  int maxit;
  double tol;
  double tuning;
  double ridge;
  double min_weight;
  rc_psi_type_t psi;
} rc_rlm_options_t;

typedef struct rc_s_options_t {
  int add_intercept;
  int n_starts;
  int n_best_starts;
  int max_refine;
  int max_scale_iter;
  int use_fast_s;
  int include_ols_start;
  int fast_s_screen_subsets;
  int fast_s_screen_iters;
  int fast_s_keep;
  double tol;
  double scale_tol;
  double c;
  double b;
  double ridge;
  double min_weight;
  unsigned long long seed;
} rc_s_options_t;

ROBUSTCAUSE_API rc_rlm_options_t rc_default_rlm_options(void);
ROBUSTCAUSE_API rc_s_options_t rc_default_s_options(void);

ROBUSTCAUSE_API rc_status_t rc_fit_rlm(
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
);

ROBUSTCAUSE_API rc_status_t rc_fit_s_estimator(
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
);

ROBUSTCAUSE_API rc_status_t rc_rlm_confint_normal(
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
);

#ifdef __cplusplus
}
#endif

#endif
