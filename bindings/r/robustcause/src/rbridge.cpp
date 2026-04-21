#include <RcppArmadillo.h>

#include <cstdint>
#include <stdexcept>
#include <string>

#include "robust/robust_rlm_hc.hpp"
#include "robust/robust_s_estimator.hpp"

// [[Rcpp::depends(RcppArmadillo)]]

namespace {

robust::PsiType parse_psi(const std::string& psi) {
  if (psi == "huber") {
    return robust::PsiType::kHuber;
  }
  if (psi == "tukey_bisquare") {
    return robust::PsiType::kTukeyBisquare;
  }
  throw std::invalid_argument("Unknown psi value: " + psi);
}

robust::HCType parse_hc(const std::string& hc_type) {
  if (hc_type == "HC0") return robust::HCType::kHC0;
  if (hc_type == "HC1") return robust::HCType::kHC1;
  if (hc_type == "HC2") return robust::HCType::kHC2;
  if (hc_type == "HC3") return robust::HCType::kHC3;
  if (hc_type == "HC4") return robust::HCType::kHC4;
  if (hc_type == "HC4m") return robust::HCType::kHC4m;
  if (hc_type == "HC5") return robust::HCType::kHC5;
  throw std::invalid_argument("Unknown HC type: " + hc_type);
}

robust::RlmMethod parse_method(const std::string& method) {
  if (method == "m") {
    return robust::RlmMethod::kM;
  }
  if (method == "mm") {
    return robust::RlmMethod::kMM;
  }
  throw std::invalid_argument("Unknown RLM method: " + method);
}

robust::SEstControl parse_s_control(const Rcpp::List& control) {
  robust::SEstControl s_ctl;
  s_ctl.add_intercept = false;
  s_ctl.n_starts = Rcpp::as<int>(control["n_starts"]);
  s_ctl.n_best_starts = Rcpp::as<int>(control["n_best_starts"]);
  s_ctl.max_refine = Rcpp::as<int>(control["max_refine"]);
  s_ctl.max_scale_iter = Rcpp::as<int>(control["max_scale_iter"]);
  s_ctl.tol = Rcpp::as<double>(control["tol"]);
  s_ctl.scale_tol = Rcpp::as<double>(control["scale_tol"]);
  s_ctl.c = Rcpp::as<double>(control["c"]);
  s_ctl.b = Rcpp::as<double>(control["b"]);
  s_ctl.ridge = Rcpp::as<double>(control["ridge"]);
  s_ctl.min_weight = Rcpp::as<double>(control["min_weight"]);
  s_ctl.seed = static_cast<std::uint64_t>(Rcpp::as<double>(control["seed"]));
  s_ctl.use_fast_s = Rcpp::as<bool>(control["use_fast_s"]);
  s_ctl.include_ols_start = Rcpp::as<bool>(control["include_ols_start"]);
  s_ctl.fast_s_screen_subsets = Rcpp::as<int>(control["fast_s_screen_subsets"]);
  s_ctl.fast_s_screen_iters = Rcpp::as<int>(control["fast_s_screen_iters"]);
  s_ctl.fast_s_keep = Rcpp::as<int>(control["fast_s_keep"]);
  return s_ctl;
}

robust::RlmControl parse_rlm_control(SEXP controlSEXP) {
  const Rcpp::List control(controlSEXP);
  robust::RlmControl ctl;
  ctl.method = parse_method(Rcpp::as<std::string>(control["method"]));
  ctl.psi = parse_psi(Rcpp::as<std::string>(control["psi"]));
  ctl.tuning = Rcpp::as<double>(control["tuning"]);
  ctl.maxit = Rcpp::as<int>(control["maxit"]);
  ctl.tol = Rcpp::as<double>(control["tol"]);
  ctl.add_intercept = Rcpp::as<bool>(control["add_intercept"]);
  ctl.ridge = Rcpp::as<double>(control["ridge"]);
  ctl.min_weight = Rcpp::as<double>(control["min_weight"]);
  ctl.mm_s_control = parse_s_control(control["mm_s_control"]);
  return ctl;
}

}  // namespace

extern "C" SEXP rc_r_fit_rlm(SEXP xSEXP,
                             SEXP ySEXP,
                             SEXP controlSEXP) {
  arma::mat X = Rcpp::as<arma::mat>(xSEXP);
  arma::vec y = Rcpp::as<arma::vec>(ySEXP);
  robust::RlmControl ctl = parse_rlm_control(controlSEXP);

  robust::RlmResult fit = robust::fit_rlm(X, y, ctl);

  return Rcpp::List::create(
    Rcpp::Named("coef") = Rcpp::wrap(fit.coef),
    Rcpp::Named("fitted") = Rcpp::wrap(fit.fitted),
    Rcpp::Named("resid") = Rcpp::wrap(fit.resid),
    Rcpp::Named("weights") = Rcpp::wrap(fit.weights),
    Rcpp::Named("hat") = Rcpp::wrap(fit.hat),
    Rcpp::Named("scale") = fit.scale,
    Rcpp::Named("converged") = fit.converged,
    Rcpp::Named("iterations") = fit.iterations,
    Rcpp::Named("method") = robust::rlm_method_name(fit.method),
    Rcpp::Named("psi") = Rcpp::as<std::string>(Rcpp::List(controlSEXP)["psi"]),
    Rcpp::Named("coef_names") = R_NilValue
  );
}

extern "C" SEXP rc_r_vcov_rlm(SEXP xSEXP,
                              SEXP ySEXP,
                              SEXP controlSEXP,
                              SEXP hcTypeSEXP) {
  arma::mat X = Rcpp::as<arma::mat>(xSEXP);
  arma::vec y = Rcpp::as<arma::vec>(ySEXP);
  robust::RlmControl ctl = parse_rlm_control(controlSEXP);

  robust::RlmResult fit = robust::fit_rlm(X, y, ctl);
  return Rcpp::wrap(robust::vcov_hc(fit, parse_hc(Rcpp::as<std::string>(hcTypeSEXP)), ctl.ridge));
}

extern "C" SEXP rc_r_confint_rlm(SEXP xSEXP,
                                 SEXP ySEXP,
                                 SEXP controlSEXP,
                                 SEXP hcTypeSEXP,
                                 SEXP zcritSEXP) {
  arma::mat X = Rcpp::as<arma::mat>(xSEXP);
  arma::vec y = Rcpp::as<arma::vec>(ySEXP);
  robust::RlmControl ctl = parse_rlm_control(controlSEXP);
  robust::RlmResult fit = robust::fit_rlm(X, y, ctl);
  robust::InferenceResult inf = robust::confint_normal(fit,
                                                       parse_hc(Rcpp::as<std::string>(hcTypeSEXP)),
                                                       Rcpp::as<double>(zcritSEXP),
                                                       ctl.ridge);
  return Rcpp::List::create(Rcpp::Named("vcov") = Rcpp::wrap(inf.vcov),
                            Rcpp::Named("se") = Rcpp::wrap(inf.se),
                            Rcpp::Named("ci") = Rcpp::wrap(inf.ci));
}

extern "C" SEXP rc_r_fit_s_estimator(SEXP xSEXP,
                                     SEXP ySEXP,
                                     SEXP addInterceptSEXP,
                                     SEXP nStartsSEXP,
                                     SEXP nBestStartsSEXP,
                                     SEXP maxRefineSEXP,
                                     SEXP maxScaleIterSEXP,
                                     SEXP tolSEXP,
                                     SEXP scaleTolSEXP,
                                     SEXP cSEXP,
                                     SEXP bSEXP,
                                     SEXP ridgeSEXP,
                                     SEXP minWeightSEXP,
                                     SEXP seedSEXP,
                                     SEXP useFastSSEXP,
                                     SEXP includeOlsStartSEXP,
                                     SEXP fastSScreenSubsetsSEXP,
                                     SEXP fastSScreenItersSEXP,
                                     SEXP fastSKeepSEXP) {
  arma::mat X = Rcpp::as<arma::mat>(xSEXP);
  arma::vec y = Rcpp::as<arma::vec>(ySEXP);

  robust::SEstControl ctl;
  ctl.add_intercept = Rcpp::as<bool>(addInterceptSEXP);
  ctl.n_starts = Rcpp::as<int>(nStartsSEXP);
  ctl.n_best_starts = Rcpp::as<int>(nBestStartsSEXP);
  ctl.max_refine = Rcpp::as<int>(maxRefineSEXP);
  ctl.max_scale_iter = Rcpp::as<int>(maxScaleIterSEXP);
  ctl.tol = Rcpp::as<double>(tolSEXP);
  ctl.scale_tol = Rcpp::as<double>(scaleTolSEXP);
  ctl.c = Rcpp::as<double>(cSEXP);
  ctl.b = Rcpp::as<double>(bSEXP);
  ctl.ridge = Rcpp::as<double>(ridgeSEXP);
  ctl.min_weight = Rcpp::as<double>(minWeightSEXP);
  ctl.seed = static_cast<std::uint64_t>(Rcpp::as<double>(seedSEXP));
  ctl.use_fast_s = Rcpp::as<bool>(useFastSSEXP);
  ctl.include_ols_start = Rcpp::as<bool>(includeOlsStartSEXP);
  ctl.fast_s_screen_subsets = Rcpp::as<int>(fastSScreenSubsetsSEXP);
  ctl.fast_s_screen_iters = Rcpp::as<int>(fastSScreenItersSEXP);
  ctl.fast_s_keep = Rcpp::as<int>(fastSKeepSEXP);

  robust::SEstResult fit = robust::fit_s_estimator(X, y, ctl);

  return Rcpp::List::create(
    Rcpp::Named("coef") = Rcpp::wrap(fit.coef),
    Rcpp::Named("fitted") = Rcpp::wrap(fit.fitted),
    Rcpp::Named("resid") = Rcpp::wrap(fit.resid),
    Rcpp::Named("weights") = Rcpp::wrap(fit.weights),
    Rcpp::Named("scale") = fit.scale,
    Rcpp::Named("objective") = fit.objective,
    Rcpp::Named("converged") = fit.converged,
    Rcpp::Named("iterations") = fit.iterations,
    Rcpp::Named("starts_tried") = fit.starts_tried,
    Rcpp::Named("starts_used") = fit.starts_used,
    Rcpp::Named("status") = static_cast<int>(fit.status),
    Rcpp::Named("message") = fit.message,
    Rcpp::Named("coef_names") = R_NilValue
  );
}
