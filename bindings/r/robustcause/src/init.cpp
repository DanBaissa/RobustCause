#include <R_ext/Rdynload.h>
#include <Rinternals.h>

extern "C" {
SEXP rc_r_confint_rlm(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP rc_r_build_adstock(SEXP, SEXP, SEXP, SEXP);
SEXP rc_r_fit_mm_dml(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP rc_r_fit_mmm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP rc_r_fit_rlm(SEXP, SEXP, SEXP);
SEXP rc_r_fit_s_estimator(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP rc_r_vcov_rlm(SEXP, SEXP, SEXP, SEXP);
}

static const R_CallMethodDef call_methods[] = {
  {"rc_r_confint_rlm", reinterpret_cast<DL_FUNC>(&rc_r_confint_rlm), 5},
  {"rc_r_build_adstock", reinterpret_cast<DL_FUNC>(&rc_r_build_adstock), 4},
  {"rc_r_fit_mm_dml", reinterpret_cast<DL_FUNC>(&rc_r_fit_mm_dml), 11},
  {"rc_r_fit_mmm", reinterpret_cast<DL_FUNC>(&rc_r_fit_mmm), 6},
  {"rc_r_fit_rlm", reinterpret_cast<DL_FUNC>(&rc_r_fit_rlm), 3},
  {"rc_r_fit_s_estimator", reinterpret_cast<DL_FUNC>(&rc_r_fit_s_estimator), 19},
  {"rc_r_vcov_rlm", reinterpret_cast<DL_FUNC>(&rc_r_vcov_rlm), 4},
  {nullptr, nullptr, 0}
};

extern "C" void R_init_robustcause(DllInfo* dll) {
  R_registerRoutines(dll, nullptr, call_methods, nullptr, nullptr);
  R_useDynamicSymbols(dll, FALSE);
}
