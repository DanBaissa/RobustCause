#include <RcppEigen.h>

#include "mm_dml/mm_dml.hpp"

using namespace Rcpp;

namespace {

void assign_if_present(const List& params, const char* name, int* target) {
  if (params.containsElementNamed(name) && !Rf_isNull(params[name])) {
    *target = as<int>(params[name]);
  }
}

void assign_if_present(const List& params, const char* name, double* target) {
  if (params.containsElementNamed(name) && !Rf_isNull(params[name])) {
    *target = as<double>(params[name]);
  }
}

void assign_if_present(const List& params, const char* name, bool* target) {
  if (params.containsElementNamed(name) && !Rf_isNull(params[name])) {
    *target = as<bool>(params[name]);
  }
}

void assign_if_present(const List& params, const char* name, std::vector<double>* target) {
  if (params.containsElementNamed(name) && !Rf_isNull(params[name])) {
    *target = as<std::vector<double>>(params[name]);
  }
}

List diagnostics_to_list(
    const mm_dml::MmDmlResult& result,
    const std::string& learner_name,
    const char* role) {
  const bool is_lasso = learner_name == "lasso" || learner_name == "elastic_net";
  const bool is_forest = learner_name == "random_forest";
  const bool is_y = std::string(role) == "y";
  SEXP selected_lambda = R_NilValue;
  SEXP lambda_min = R_NilValue;
  SEXP lambda_1se = R_NilValue;
  SEXP cv_error = R_NilValue;
  SEXP nonzero = R_NilValue;
  SEXP oob_error = R_NilValue;
  SEXP num_trees = R_NilValue;
  SEXP feature_importance = R_NilValue;
  SEXP feature_use_count = R_NilValue;

  if (is_lasso) {
    selected_lambda = wrap(is_y ? result.y_selected_lambda : result.d_selected_lambda);
    lambda_min = wrap(is_y ? result.y_lambda_min : result.d_lambda_min);
    lambda_1se = wrap(is_y ? result.y_lambda_1se : result.d_lambda_1se);
    cv_error = wrap(is_y ? result.y_cv_error : result.d_cv_error);
    nonzero = wrap(is_y ? result.y_nonzero : result.d_nonzero);
  }
  if (is_forest) {
    oob_error = wrap(is_y ? result.y_oob_error : result.d_oob_error);
    num_trees = wrap(is_y ? result.y_num_trees : result.d_num_trees);
    feature_importance = wrap(is_y ? result.y_feature_importance : result.d_feature_importance);
    feature_use_count = wrap(is_y ? result.y_feature_use_count : result.d_feature_use_count);
  }

  return List::create(
      _["role"] = role,
      _["selected_lambda"] = selected_lambda,
      _["lambda_min"] = lambda_min,
      _["lambda_1se"] = lambda_1se,
      _["cv_error"] = cv_error,
      _["nonzero"] = nonzero,
      _["oob_error"] = oob_error,
      _["num_trees"] = num_trees,
      _["feature_importance"] = feature_importance,
      _["feature_use_count"] = feature_use_count
  );
}

}  // namespace

extern "C" SEXP rc_r_fit_mm_dml(
    SEXP x,
    SEXP d,
    SEXP y,
    SEXP learner,
    SEXP learner_params,
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

  const std::string learner_name = as<std::string>(learner);
  config.learner.kind = mm_dml::learner_kind_from_name(learner_name);
  List params(learner_params);
  assign_if_present(params, "cv_folds", &config.learner.cv_folds);
  assign_if_present(params, "n_estimators", &config.learner.n_estimators);
  assign_if_present(params, "max_depth", &config.learner.max_depth);
  assign_if_present(params, "min_samples_split", &config.learner.min_samples_split);
  assign_if_present(params, "min_samples_leaf", &config.learner.min_samples_leaf);
  assign_if_present(params, "max_features", &config.learner.max_features);
  assign_if_present(params, "max_features_fraction", &config.learner.max_features_fraction);
  assign_if_present(params, "sample_fraction", &config.learner.sample_fraction);
  assign_if_present(params, "bootstrap", &config.learner.bootstrap);
  assign_if_present(params, "replacement", &config.learner.replacement);
  assign_if_present(params, "split_candidates", &config.learner.split_candidates);
  assign_if_present(params, "compute_oob", &config.learner.compute_oob);
  assign_if_present(params, "compute_importance", &config.learner.compute_importance);
  assign_if_present(params, "max_iter", &config.learner.max_iter);
  assign_if_present(params, "n_lambda", &config.learner.n_lambda);
  assign_if_present(params, "learning_rate", &config.learner.learning_rate);
  assign_if_present(params, "lambda_min_ratio", &config.learner.lambda_min_ratio);
  assign_if_present(params, "tolerance", &config.learner.tolerance);
  assign_if_present(params, "standardize", &config.learner.standardize);
  assign_if_present(params, "use_lambda_1se", &config.learner.use_lambda_1se);
  assign_if_present(params, "lambda_grid", &config.learner.lambda_grid);
  assign_if_present(params, "l1_ratios", &config.learner.l1_ratios);

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
      _["learner_name"] = learner_name,
      _["learner_details"] = List::create(
          _["y"] = diagnostics_to_list(result, learner_name, "y"),
          _["d"] = diagnostics_to_list(result, learner_name, "d")
      ),
      _["backend"] = "cpp"
  );
  END_RCPP
}
