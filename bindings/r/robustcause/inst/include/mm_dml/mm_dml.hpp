#pragma once

#include <string>
#include <vector>

#include <Eigen/Dense>

namespace mm_dml {

enum class LearnerKind {
    kLasso,
    kElasticNet,
    kRandomForest,
    kHistGradientBoosting,
};

struct LearnerSpec {
    LearnerKind kind{LearnerKind::kLasso};
    int cv_folds{3};
    int n_estimators{300};
    int max_depth{8};
    int min_samples_split{10};
    int min_samples_leaf{5};
    int max_features{0};
    double max_features_fraction{0.0};
    double sample_fraction{0.632};
    bool bootstrap{true};
    bool replacement{true};
    int split_candidates{16};
    bool compute_oob{true};
    bool compute_importance{true};
    int max_iter{250};
    int n_lambda{100};
    double learning_rate{0.05};
    double lambda_min_ratio{1e-4};
    double tolerance{1e-6};
    bool standardize{true};
    bool use_lambda_1se{false};
    std::vector<double> lambda_grid{};
    std::vector<double> l1_ratios{0.1, 0.5, 0.9, 0.95, 1.0};
};

struct DmlConfig {
    int folds{2};
    int seed{123};
    LearnerSpec learner{};
    double mm_c{4.685};
    int max_iter{100};
    double tolerance{1e-7};
    double ci_level{0.95};
};

struct FitInput {
    Eigen::MatrixXd x;
    Eigen::VectorXd d;
    Eigen::VectorXd y;
};

struct MmDmlResult {
    double tau_hat{0.0};
    double intercept_hat{0.0};
    double std_error{0.0};
    double ci_lower{0.0};
    double ci_upper{0.0};
    int iterations{0};
    bool converged{false};
    double objective{0.0};
    double kept_fraction{1.0};
    double residual_outcome_scale{0.0};
    double residual_treatment_scale{0.0};
    double y_selected_lambda{0.0};
    double d_selected_lambda{0.0};
    double y_lambda_min{0.0};
    double d_lambda_min{0.0};
    double y_lambda_1se{0.0};
    double d_lambda_1se{0.0};
    double y_cv_error{0.0};
    double d_cv_error{0.0};
    double y_oob_error{0.0};
    double d_oob_error{0.0};
    int y_nonzero{0};
    int d_nonzero{0};
    int y_num_trees{0};
    int d_num_trees{0};
    std::vector<double> y_feature_importance;
    std::vector<double> d_feature_importance;
    std::vector<int> y_feature_use_count;
    std::vector<int> d_feature_use_count;
    std::vector<double> y_residuals;
    std::vector<double> d_residuals;
    std::vector<double> stage_residuals;
};

LearnerKind learner_kind_from_name(const std::string& name);
MmDmlResult fit_mm_dml(const FitInput& input, const DmlConfig& config);

}  // namespace mm_dml
