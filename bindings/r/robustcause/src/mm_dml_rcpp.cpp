#include <RcppEigen.h>

#include "mm_dml/mm_dml.hpp"

using namespace Rcpp;

extern "C" SEXP rc_r_fit_mm_dml(
    SEXP x,
    SEXP d,
    SEXP y,
    SEXP learner,
    SEXP folds,
    SEXP seed,
    SEXP mm_c,
    SEXP max_iter,
    SEXP tolerance,
    SEXP ci_level) {
  BEGIN_RCPP
  NumericMatrix x_in(x);
  NumericVector d_in(d);
  NumericVector y_in(y);
  Eigen::Map<Eigen::MatrixXd> x_map(REAL(x_in), x_in.nrow(), x_in.ncol());
  Eigen::Map<Eigen::VectorXd> d_map(REAL(d_in), d_in.size());
  Eigen::Map<Eigen::VectorXd> y_map(REAL(y_in), y_in.size());

  mm_dml::FitInput input{x_map, d_map, y_map};
  mm_dml::DmlConfig config;
  config.folds = as<int>(folds);
  config.seed = as<int>(seed);
  config.mm_c = as<double>(mm_c);
  config.max_iter = as<int>(max_iter);
  config.tolerance = as<double>(tolerance);
  config.ci_level = as<double>(ci_level);
  config.learner.kind = mm_dml::learner_kind_from_name(as<std::string>(learner));
  const auto result = mm_dml::fit_mm_dml(input, config);

  return List::create(
      _["tau_hat"] = result.tau_hat,
      _["intercept_hat"] = result.intercept_hat,
      _["std_error"] = result.std_error,
      _["ci_lower"] = result.ci_lower,
      _["ci_upper"] = result.ci_upper,
      _["iterations"] = result.iterations,
      _["converged"] = result.converged,
      _["objective"] = result.objective,
      _["kept_fraction"] = result.kept_fraction,
      _["residual_outcome_scale"] = result.residual_outcome_scale,
      _["residual_treatment_scale"] = result.residual_treatment_scale,
      _["y_residuals"] = result.y_residuals,
      _["d_residuals"] = result.d_residuals,
      _["stage_residuals"] = result.stage_residuals,
      _["learner_name"] = as<std::string>(learner),
      _["backend"] = "cpp"
  );
  END_RCPP
}
