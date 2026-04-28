#include <RcppArmadillo.h>
#include <Rinternals.h>

#include <cstdint>
#include <stdexcept>
#include <string>

#include "robust/robust_rlm_hc.hpp"
#include "robust/robust_sc.hpp"
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
  if (psi == "hampel") {
    return robust::PsiType::kHampel;
  }
  throw std::invalid_argument("Unknown psi value: " + psi);
}

robust::WeightMethod parse_weight_method(const std::string& wt_method) {
  if (wt_method == "inv.var") {
    return robust::WeightMethod::kInvVar;
  }
  if (wt_method == "case") {
    return robust::WeightMethod::kCase;
  }
  throw std::invalid_argument("Unknown weight method: " + wt_method);
}

robust::InitMethod parse_init_method(const std::string& init_method) {
  if (init_method == "ls") {
    return robust::InitMethod::kOLS;
  }
  if (init_method == "lts") {
    return robust::InitMethod::kLTS;
  }
  if (init_method == "user") {
    return robust::InitMethod::kUser;
  }
  throw std::invalid_argument("Unknown init method: " + init_method);
}

robust::ScaleEstimate parse_scale_estimate(const std::string& scale_est) {
  if (scale_est == "MAD") {
    return robust::ScaleEstimate::kMAD;
  }
  if (scale_est == "Huber") {
    return robust::ScaleEstimate::kHuber;
  }
  if (scale_est == "proposal 2") {
    return robust::ScaleEstimate::kProposal2;
  }
  throw std::invalid_argument("Unknown scale estimate: " + scale_est);
}

robust::TestVector parse_test_vector(const std::string& test_vec) {
  if (test_vec == "coef") {
    return robust::TestVector::kCoef;
  }
  if (test_vec == "resid") {
    return robust::TestVector::kResid;
  }
  if (test_vec == "w") {
    return robust::TestVector::kWeight;
  }
  throw std::invalid_argument("Unknown test vector: " + test_vec);
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
  ctl.acc = Rcpp::as<double>(control["acc"]);
  ctl.add_intercept = Rcpp::as<bool>(control["add_intercept"]);
  ctl.ridge = Rcpp::as<double>(control["ridge"]);
  ctl.min_weight = Rcpp::as<double>(control["min_weight"]);
  ctl.prior_weights = Rcpp::as<arma::vec>(control["prior_weights"]);
  ctl.wt_method = parse_weight_method(Rcpp::as<std::string>(control["wt_method"]));
  ctl.init_method = parse_init_method(Rcpp::as<std::string>(control["init_method"]));
  ctl.init_coef = Rcpp::as<arma::vec>(control["init_coef"]);
  ctl.scale_est = parse_scale_estimate(Rcpp::as<std::string>(control["scale_est"]));
  ctl.k2 = Rcpp::as<double>(control["k2"]);
  ctl.test_vec = parse_test_vector(Rcpp::as<std::string>(control["test_vec"]));
  ctl.mm_s_control = parse_s_control(control["mm_s_control"]);
  return ctl;
}

robust::ScControl make_sc_control(SEXP maxitSEXP,
                                  SEXP tolSEXP,
                                  SEXP ridgeSEXP) {
  robust::ScControl ctl;
  ctl.maxit = Rcpp::as<int>(maxitSEXP);
  ctl.tol = Rcpp::as<double>(tolSEXP);
  ctl.ridge = Rcpp::as<double>(ridgeSEXP);
  return ctl;
}

robust::MmscControl make_mmsc_control(SEXP startupMaxitSEXP,
                                      SEXP maxitSEXP,
                                      SEXP subproblemMaxitSEXP,
                                      SEXP tolSEXP,
                                      SEXP subproblemTolSEXP,
                                      SEXP l1SmoothingSEXP,
                                      SEXP tukeyCSEXP,
                                      SEXP minScaleSEXP,
                                      SEXP ridgeSEXP,
                                      SEXP minTimeWeightSEXP) {
  robust::MmscControl ctl;
  ctl.startup_maxit = Rcpp::as<int>(startupMaxitSEXP);
  ctl.maxit = Rcpp::as<int>(maxitSEXP);
  ctl.subproblem_maxit = Rcpp::as<int>(subproblemMaxitSEXP);
  ctl.tol = Rcpp::as<double>(tolSEXP);
  ctl.subproblem_tol = Rcpp::as<double>(subproblemTolSEXP);
  ctl.l1_smoothing = Rcpp::as<double>(l1SmoothingSEXP);
  ctl.tukey_c = Rcpp::as<double>(tukeyCSEXP);
  ctl.min_scale = Rcpp::as<double>(minScaleSEXP);
  ctl.ridge = Rcpp::as<double>(ridgeSEXP);
  ctl.min_time_weight = Rcpp::as<double>(minTimeWeightSEXP);
  return ctl;
}

Rcpp::List wrap_sc_result(const robust::ScResult& fit) {
  return Rcpp::List::create(
    Rcpp::Named("weights") = Rcpp::wrap(fit.weights),
    Rcpp::Named("synthetic_pre") = Rcpp::wrap(fit.synthetic_pre),
    Rcpp::Named("synthetic_post") = Rcpp::wrap(fit.synthetic_post),
    Rcpp::Named("pre_residuals") = Rcpp::wrap(fit.pre_residuals),
    Rcpp::Named("post_gaps") = Rcpp::wrap(fit.post_gaps),
    Rcpp::Named("objective_value") = fit.objective_value,
    Rcpp::Named("pre_rmspe") = fit.pre_rmspe,
    Rcpp::Named("weight_herfindahl") = fit.weight_herfindahl,
    Rcpp::Named("max_weight") = fit.max_weight,
    Rcpp::Named("effective_donors") = fit.effective_donors,
    Rcpp::Named("converged") = fit.converged,
    Rcpp::Named("iterations") = fit.iterations
  );
}

Rcpp::List wrap_mmsc_result(const robust::MmscResult& fit) {
  Rcpp::List out;
  out["weights"] = Rcpp::wrap(fit.weights);
  out["startup_weights"] = Rcpp::wrap(fit.startup_weights);
  out["synthetic_pre"] = Rcpp::wrap(fit.synthetic_pre);
  out["synthetic_post"] = Rcpp::wrap(fit.synthetic_post);
  out["pre_residuals"] = Rcpp::wrap(fit.pre_residuals);
  out["post_gaps"] = Rcpp::wrap(fit.post_gaps);
  out["robust_time_weights"] = Rcpp::wrap(fit.robust_time_weights);
  out["downweighted_periods"] = Rcpp::wrap(fit.downweighted_periods);
  out["scale"] = fit.scale;
  out["objective_value"] = fit.objective_value;
  out["startup_objective"] = fit.startup_objective;
  out["standard_rmspe"] = fit.standard_rmspe;
  out["startup_rmspe"] = fit.startup_rmspe;
  out["mm_rmspe"] = fit.mm_rmspe;
  out["effective_periods"] = fit.effective_periods;
  out["weight_herfindahl"] = fit.weight_herfindahl;
  out["max_weight"] = fit.max_weight;
  out["effective_donors"] = fit.effective_donors;
  out["converged"] = fit.converged;
  out["iterations"] = fit.iterations;
  out["standard_fit"] = wrap_sc_result(fit.standard_fit);
  return out;
}

Rcpp::List wrap_placebo_result(const robust::ScPlaceboResult& fit) {
  return Rcpp::List::create(
    Rcpp::Named("unit_indices") = Rcpp::wrap(fit.unit_indices),
    Rcpp::Named("weight_matrix") = Rcpp::wrap(fit.weight_matrix),
    Rcpp::Named("post_gap_matrix") = Rcpp::wrap(fit.post_gap_matrix),
    Rcpp::Named("pre_rmspe") = Rcpp::wrap(fit.pre_rmspe),
    Rcpp::Named("post_rmspe") = Rcpp::wrap(fit.post_rmspe),
    Rcpp::Named("rmspe_ratio") = Rcpp::wrap(fit.rmspe_ratio)
  );
}

}  // namespace

extern "C" SEXP rc_r_fit_rlm(SEXP xSEXP,
                             SEXP ySEXP,
                             SEXP controlSEXP) {
  arma::mat X = Rcpp::as<arma::mat>(xSEXP);
  arma::vec y = Rcpp::as<arma::vec>(ySEXP);
  robust::RlmControl ctl = parse_rlm_control(controlSEXP);

  robust::RlmResult fit = robust::fit_rlm(X, y, ctl);
  const std::string psi_name = robust::psi_name(parse_psi(Rcpp::as<std::string>(Rcpp::List(controlSEXP)["psi"])));

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
    Rcpp::Named("psi") = psi_name,
    Rcpp::Named("prior_weights") = Rcpp::wrap(fit.prior_weights),
    Rcpp::Named("wt_method") = robust::weight_method_name(fit.wt_method),
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

extern "C" SEXP rc_r_fit_sc(SEXP treatedPreSEXP,
                            SEXP donorsPreSEXP,
                            SEXP treatedPostSEXP,
                            SEXP donorsPostSEXP,
                            SEXP maxitSEXP,
                            SEXP tolSEXP,
                            SEXP ridgeSEXP,
                            SEXP predictorTreatedSEXP,
                            SEXP predictorDonorsSEXP,
                            SEXP predictorWeightsSEXP) {
  try {
    arma::vec treated_pre = Rcpp::as<arma::vec>(treatedPreSEXP);
    arma::mat donors_pre = Rcpp::as<arma::mat>(donorsPreSEXP);
    arma::vec treated_post = Rcpp::as<arma::vec>(treatedPostSEXP);
    arma::mat donors_post = Rcpp::as<arma::mat>(donorsPostSEXP);
    arma::vec predictor_treated = Rcpp::as<arma::vec>(predictorTreatedSEXP);
    arma::mat predictor_donors = Rcpp::as<arma::mat>(predictorDonorsSEXP);
    arma::vec predictor_weights = Rcpp::as<arma::vec>(predictorWeightsSEXP);
    robust::ScControl ctl = make_sc_control(maxitSEXP, tolSEXP, ridgeSEXP);
    return wrap_sc_result(robust::fit_sc(
      treated_pre, donors_pre, treated_post, donors_post, ctl,
      predictor_treated, predictor_donors, predictor_weights));
  } catch (const std::exception& e) {
    Rcpp::stop("rc_r_fit_sc failed: %s", e.what());
  } catch (...) {
    Rcpp::stop("rc_r_fit_sc failed with an unknown C++ exception");
  }
}

extern "C" SEXP rc_r_fit_mm_sc(SEXP treatedPreSEXP,
                               SEXP donorsPreSEXP,
                               SEXP treatedPostSEXP,
                               SEXP donorsPostSEXP,
                               SEXP startupMaxitSEXP,
                               SEXP maxitSEXP,
                               SEXP subproblemMaxitSEXP,
                               SEXP tolSEXP,
                               SEXP subproblemTolSEXP,
                               SEXP l1SmoothingSEXP,
                               SEXP tukeyCSEXP,
                               SEXP minScaleSEXP,
                               SEXP ridgeSEXP,
                               SEXP minTimeWeightSEXP,
                               SEXP predictorTreatedSEXP,
                               SEXP predictorDonorsSEXP,
                               SEXP predictorWeightsSEXP) {
  try {
    arma::vec treated_pre = Rcpp::as<arma::vec>(treatedPreSEXP);
    arma::mat donors_pre = Rcpp::as<arma::mat>(donorsPreSEXP);
    arma::vec treated_post = Rcpp::as<arma::vec>(treatedPostSEXP);
    arma::mat donors_post = Rcpp::as<arma::mat>(donorsPostSEXP);
    arma::vec predictor_treated = Rcpp::as<arma::vec>(predictorTreatedSEXP);
    arma::mat predictor_donors = Rcpp::as<arma::mat>(predictorDonorsSEXP);
    arma::vec predictor_weights = Rcpp::as<arma::vec>(predictorWeightsSEXP);
    robust::MmscControl ctl = make_mmsc_control(
      startupMaxitSEXP,
      maxitSEXP,
      subproblemMaxitSEXP,
      tolSEXP,
      subproblemTolSEXP,
      l1SmoothingSEXP,
      tukeyCSEXP,
      minScaleSEXP,
      ridgeSEXP,
      minTimeWeightSEXP);
    return wrap_mmsc_result(robust::fit_mm_sc(
      treated_pre, donors_pre, treated_post, donors_post, ctl,
      predictor_treated, predictor_donors, predictor_weights));
  } catch (const std::exception& e) {
    Rcpp::stop("rc_r_fit_mm_sc failed: %s", e.what());
  } catch (...) {
    Rcpp::stop("rc_r_fit_mm_sc failed with an unknown C++ exception");
  }
}

extern "C" SEXP rc_r_fit_sc_placebos(SEXP outcomesSEXP,
                                     SEXP treatmentStartSEXP,
                                     SEXP maxitSEXP,
                                     SEXP tolSEXP,
                                     SEXP ridgeSEXP,
                                     SEXP startupMaxitSEXP,
                                     SEXP subproblemMaxitSEXP,
                                     SEXP subproblemTolSEXP,
                                     SEXP l1SmoothingSEXP,
                                     SEXP tukeyCSEXP,
                                     SEXP minScaleSEXP,
                                     SEXP minTimeWeightSEXP,
                                     SEXP useMmSEXP) {
  try {
    arma::mat outcomes = Rcpp::as<arma::mat>(outcomesSEXP);
    const int treatment_start = Rcpp::as<int>(treatmentStartSEXP);
    const bool use_mm = Rcpp::as<bool>(useMmSEXP);
    if (use_mm) {
      robust::MmscControl ctl = make_mmsc_control(
        startupMaxitSEXP, maxitSEXP, subproblemMaxitSEXP, tolSEXP, subproblemTolSEXP,
        l1SmoothingSEXP, tukeyCSEXP, minScaleSEXP, ridgeSEXP, minTimeWeightSEXP);
      return wrap_placebo_result(robust::fit_mm_sc_placebos(outcomes, treatment_start, ctl));
    }
    robust::ScControl ctl = make_sc_control(maxitSEXP, tolSEXP, ridgeSEXP);
    return wrap_placebo_result(robust::fit_sc_placebos(outcomes, treatment_start, ctl));
  } catch (const std::exception& e) {
    Rcpp::stop("rc_r_fit_sc_placebos failed: %s", e.what());
  } catch (...) {
    Rcpp::stop("rc_r_fit_sc_placebos failed with an unknown C++ exception");
  }
}
