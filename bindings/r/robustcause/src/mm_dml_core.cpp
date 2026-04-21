#include "mm_dml/mm_dml.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <numeric>
#include <random>
#include <stdexcept>
#include <utility>

namespace mm_dml {
namespace {

using Eigen::MatrixXd;
using Eigen::VectorXd;

constexpr double kEps = 1e-8;

double standard_deviation(const VectorXd& values) {
    if (values.size() < 2) {
        return 0.0;
    }
    const double mu = values.mean();
    return std::sqrt(std::max((values.array() - mu).square().mean(), 0.0));
}

double median(std::vector<double> values) {
    if (values.empty()) {
        return 0.0;
    }
    const size_t mid = values.size() / 2;
    std::nth_element(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(mid), values.end());
    double m = values[mid];
    if (values.size() % 2 == 0) {
        std::nth_element(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(mid - 1), values.end());
        m = 0.5 * (m + values[mid - 1]);
    }
    return m;
}

double mad_scale(const VectorXd& residuals) {
    std::vector<double> values(residuals.data(), residuals.data() + residuals.size());
    const double med = median(values);
    std::vector<double> abs_dev(residuals.size());
    for (Eigen::Index i = 0; i < residuals.size(); ++i) {
        abs_dev[static_cast<size_t>(i)] = std::abs(residuals[i] - med);
    }
    return std::max(1.4826 * median(abs_dev), kEps);
}

double normal_quantile(double p) {
    if (p <= 0.0 || p >= 1.0) {
        throw std::invalid_argument("normal_quantile requires 0 < p < 1");
    }
    static const double a[] = {-3.969683028665376e+01, 2.209460984245205e+02,
        -2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01,
        2.506628277459239e+00};
    static const double b[] = {-5.447609879822406e+01, 1.615858368580409e+02,
        -1.556989798598866e+02, 6.680131188771972e+01, -1.328068155288572e+01};
    static const double c[] = {-7.784894002430293e-03, -3.223964580411365e-01,
        -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00,
        2.938163982698783e+00};
    static const double d[] = {7.784695709041462e-03, 3.224671290700398e-01,
        2.445134137142996e+00, 3.754408661907416e+00};
    const double plow = 0.02425;
    const double phigh = 1.0 - plow;
    if (p < plow) {
        const double q = std::sqrt(-2.0 * std::log(p));
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
               ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
    }
    if (p > phigh) {
        const double q = std::sqrt(-2.0 * std::log(1.0 - p));
        return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
               ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
    }
    const double q = p - 0.5;
    const double r = q * q;
    return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
           (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0);
}

VectorXd select_rows(const VectorXd& vector, const std::vector<int>& indices) {
    VectorXd out(indices.size());
    for (size_t i = 0; i < indices.size(); ++i) {
        out[static_cast<Eigen::Index>(i)] = vector[indices[i]];
    }
    return out;
}

MatrixXd select_rows(const MatrixXd& matrix, const std::vector<int>& indices) {
    MatrixXd out(indices.size(), matrix.cols());
    for (size_t i = 0; i < indices.size(); ++i) {
        out.row(static_cast<Eigen::Index>(i)) = matrix.row(indices[i]);
    }
    return out;
}

std::vector<std::vector<int>> make_folds(int n, int folds, int seed) {
    std::vector<int> indices(n);
    std::iota(indices.begin(), indices.end(), 0);
    std::mt19937 rng(seed);
    std::shuffle(indices.begin(), indices.end(), rng);
    std::vector<std::vector<int>> out(static_cast<size_t>(folds));
    for (int i = 0; i < n; ++i) {
        out[static_cast<size_t>(i % folds)].push_back(indices[static_cast<size_t>(i)]);
    }
    return out;
}

struct StandardizedDesign {
    MatrixXd x_scaled;
    VectorXd means;
    VectorXd scales;
    double y_mean{0.0};
};

StandardizedDesign standardize(const MatrixXd& x, const VectorXd& y) {
    StandardizedDesign out;
    out.x_scaled = x;
    out.means = x.colwise().mean().transpose();
    out.scales = VectorXd::Ones(x.cols());
    for (Eigen::Index j = 0; j < x.cols(); ++j) {
        VectorXd centered = x.col(j).array() - out.means[j];
        const double scale = std::sqrt(std::max(centered.array().square().mean(), 0.0));
        out.scales[j] = scale > kEps ? scale : 1.0;
        out.x_scaled.col(j) = centered / out.scales[j];
    }
    out.y_mean = y.mean();
    return out;
}

struct RegressionModel {
    virtual ~RegressionModel() = default;
    virtual VectorXd predict(const MatrixXd& x) const = 0;
};

struct MeanModel final : RegressionModel {
    explicit MeanModel(double value) : value(value) {}
    double value;
    VectorXd predict(const MatrixXd& x) const override {
        return VectorXd::Constant(x.rows(), value);
    }
};

struct LinearModel final : RegressionModel {
    VectorXd beta;
    double intercept{0.0};
    VectorXd predict(const MatrixXd& x) const override {
        return (x * beta).array() + intercept;
    }
};

std::unique_ptr<LinearModel> fit_elastic_net_internal(
    const MatrixXd& x,
    const VectorXd& y,
    double l1_ratio,
    int cv_folds,
    int max_iter) {
    const auto solve_for_lambda = [&](const MatrixXd& x_train, const VectorXd& y_train, double lambda) {
        const auto prep = standardize(x_train, y_train);
        const MatrixXd& xs = prep.x_scaled;
        const VectorXd yc = y_train.array() - prep.y_mean;
        const VectorXd z = xs.array().square().colwise().sum().transpose();
        VectorXd beta = VectorXd::Zero(x_train.cols());
        VectorXd residual = yc;
        for (int iter = 0; iter < max_iter; ++iter) {
            double max_delta = 0.0;
            for (Eigen::Index j = 0; j < x_train.cols(); ++j) {
                residual.array() += xs.col(j).array() * beta[j];
                const double rho = xs.col(j).dot(residual);
                const double denom = z[j] + lambda * (1.0 - l1_ratio) * x_train.rows();
                const double thresh = lambda * l1_ratio * x_train.rows();
                double updated = 0.0;
                if (rho > thresh) {
                    updated = (rho - thresh) / std::max(denom, kEps);
                } else if (rho < -thresh) {
                    updated = (rho + thresh) / std::max(denom, kEps);
                }
                residual.array() -= xs.col(j).array() * updated;
                max_delta = std::max(max_delta, std::abs(updated - beta[j]));
                beta[j] = updated;
            }
            if (max_delta < 1e-6) {
                break;
            }
        }
        auto model = std::make_unique<LinearModel>();
        model->beta = beta.array() / prep.scales.array();
        model->intercept = prep.y_mean - model->beta.dot(prep.means);
        return model;
    };

    if (x.rows() < 3 || x.cols() == 0) {
        auto model = std::make_unique<LinearModel>();
        model->beta = VectorXd::Zero(x.cols());
        model->intercept = y.size() ? y.mean() : 0.0;
        return model;
    }

    const auto prep = standardize(x, y);
    const VectorXd yc = y.array() - prep.y_mean;
    const double max_lambda = std::max((prep.x_scaled.transpose() * yc).cwiseAbs().maxCoeff() /
                                       std::max(1.0, static_cast<double>(x.rows()) * std::max(l1_ratio, 1e-3)),
                                       1e-3);
    std::vector<double> lambda_grid;
    for (int i = 0; i < 20; ++i) {
        lambda_grid.push_back(max_lambda * std::pow(1e-3, static_cast<double>(i) / 19.0));
    }

    const auto folds = make_folds(static_cast<int>(x.rows()), std::min(cv_folds, static_cast<int>(x.rows())), 911 + static_cast<int>(std::round(100 * l1_ratio)));
    double best_loss = std::numeric_limits<double>::infinity();
    double best_lambda = lambda_grid.front();
    for (double lambda : lambda_grid) {
        double loss = 0.0;
        for (const auto& test_idx : folds) {
            std::vector<char> mask(static_cast<size_t>(x.rows()), 1);
            for (int idx : test_idx) {
                mask[static_cast<size_t>(idx)] = 0;
            }
            std::vector<int> train_idx;
            train_idx.reserve(x.rows() - static_cast<int>(test_idx.size()));
            for (int row = 0; row < x.rows(); ++row) {
                if (mask[static_cast<size_t>(row)]) {
                    train_idx.push_back(row);
                }
            }
            auto model = solve_for_lambda(select_rows(x, train_idx), select_rows(y, train_idx), lambda);
            const VectorXd err = select_rows(y, test_idx) - model->predict(select_rows(x, test_idx));
            loss += err.array().square().mean();
        }
        loss /= static_cast<double>(folds.size());
        if (loss < best_loss) {
            best_loss = loss;
            best_lambda = lambda;
        }
    }
    return solve_for_lambda(x, y, best_lambda);
}

struct RegressionTree final : RegressionModel {
    struct Node {
        bool leaf{true};
        int feature{-1};
        double threshold{0.0};
        double value{0.0};
        std::unique_ptr<Node> left;
        std::unique_ptr<Node> right;
    };

    std::unique_ptr<Node> root;

    double predict_row(const Node* node, const Eigen::Ref<const Eigen::RowVectorXd>& row) const {
        if (node->leaf || !node->left || !node->right) {
            return node->value;
        }
        return row[node->feature] <= node->threshold ? predict_row(node->left.get(), row) : predict_row(node->right.get(), row);
    }

    VectorXd predict(const MatrixXd& x) const override {
        VectorXd out(x.rows());
        for (Eigen::Index i = 0; i < x.rows(); ++i) {
            out[i] = predict_row(root.get(), x.row(i));
        }
        return out;
    }
};

double subset_loss(const VectorXd& y, const std::vector<int>& subset) {
    double subset_mean = 0.0;
    for (int row : subset) {
        subset_mean += y[row];
    }
    subset_mean /= static_cast<double>(subset.size());
    double loss = 0.0;
    for (int row : subset) {
        const double diff = y[row] - subset_mean;
        loss += diff * diff;
    }
    return loss;
}

std::unique_ptr<RegressionTree::Node> build_tree(
    const MatrixXd& x,
    const VectorXd& y,
    const std::vector<int>& rows,
    int depth,
    int max_depth,
    int min_leaf,
    int max_features,
    std::mt19937& rng) {
    auto node = std::make_unique<RegressionTree::Node>();
    for (int row : rows) {
        node->value += y[row];
    }
    node->value /= static_cast<double>(rows.size());
    if (depth >= max_depth || static_cast<int>(rows.size()) <= 2 * min_leaf) {
        return node;
    }

    std::vector<int> features(x.cols());
    std::iota(features.begin(), features.end(), 0);
    std::shuffle(features.begin(), features.end(), rng);
    if (max_features > 0 && max_features < static_cast<int>(features.size())) {
        features.resize(static_cast<size_t>(max_features));
    }

    double best_loss = std::numeric_limits<double>::infinity();
    int best_feature = -1;
    double best_threshold = 0.0;
    std::vector<int> best_left;
    std::vector<int> best_right;
    for (int feature : features) {
        std::vector<double> values;
        values.reserve(rows.size());
        for (int row : rows) {
            values.push_back(x(row, feature));
        }
        std::sort(values.begin(), values.end());
        for (int q = 1; q <= 8; ++q) {
            const size_t idx = static_cast<size_t>(std::clamp(static_cast<int>(std::round(q * (values.size() - 1) / 9.0)), 0, static_cast<int>(values.size() - 1)));
            const double threshold = values[idx];
            std::vector<int> left;
            std::vector<int> right;
            for (int row : rows) {
                if (x(row, feature) <= threshold) {
                    left.push_back(row);
                } else {
                    right.push_back(row);
                }
            }
            if (static_cast<int>(left.size()) < min_leaf || static_cast<int>(right.size()) < min_leaf) {
                continue;
            }
            const double loss = subset_loss(y, left) + subset_loss(y, right);
            if (loss < best_loss) {
                best_loss = loss;
                best_feature = feature;
                best_threshold = threshold;
                best_left = left;
                best_right = right;
            }
        }
    }

    if (best_feature < 0) {
        return node;
    }
    node->leaf = false;
    node->feature = best_feature;
    node->threshold = best_threshold;
    node->left = build_tree(x, y, best_left, depth + 1, max_depth, min_leaf, max_features, rng);
    node->right = build_tree(x, y, best_right, depth + 1, max_depth, min_leaf, max_features, rng);
    return node;
}

std::unique_ptr<RegressionTree> fit_tree(
    const MatrixXd& x,
    const VectorXd& y,
    int max_depth,
    int min_leaf,
    int max_features,
    std::mt19937& rng,
    bool bootstrap) {
    std::vector<int> rows;
    rows.reserve(x.rows());
    if (bootstrap) {
        std::uniform_int_distribution<int> dist(0, static_cast<int>(x.rows()) - 1);
        for (Eigen::Index i = 0; i < x.rows(); ++i) {
            rows.push_back(dist(rng));
        }
    } else {
        for (Eigen::Index i = 0; i < x.rows(); ++i) {
            rows.push_back(static_cast<int>(i));
        }
    }
    auto tree = std::make_unique<RegressionTree>();
    tree->root = build_tree(x, y, rows, 0, max_depth, min_leaf, max_features, rng);
    return tree;
}

struct RandomForestModel final : RegressionModel {
    std::vector<std::unique_ptr<RegressionTree>> trees;
    VectorXd predict(const MatrixXd& x) const override {
        VectorXd out = VectorXd::Zero(x.rows());
        for (const auto& tree : trees) {
            out += tree->predict(x);
        }
        return trees.empty() ? out : out / static_cast<double>(trees.size());
    }
};

struct GradientBoostingModel final : RegressionModel {
    double bias{0.0};
    double learning_rate{0.05};
    std::vector<std::unique_ptr<RegressionTree>> trees;
    VectorXd predict(const MatrixXd& x) const override {
        VectorXd out = VectorXd::Constant(x.rows(), bias);
        for (const auto& tree : trees) {
            out.array() += learning_rate * tree->predict(x).array();
        }
        return out;
    }
};

std::unique_ptr<RandomForestModel> fit_random_forest(const MatrixXd& x, const VectorXd& y, const LearnerSpec& spec, int seed) {
    auto model = std::make_unique<RandomForestModel>();
    std::mt19937 rng(seed);
    const int max_features = spec.max_features > 0 ? spec.max_features : std::max(1, static_cast<int>(std::sqrt(static_cast<double>(x.cols()))));
    for (int t = 0; t < spec.n_estimators; ++t) {
        model->trees.push_back(fit_tree(x, y, spec.max_depth, spec.min_samples_leaf, max_features, rng, true));
    }
    return model;
}

std::unique_ptr<GradientBoostingModel> fit_hist_gradient_boosting(const MatrixXd& x, const VectorXd& y, const LearnerSpec& spec, int seed) {
    auto model = std::make_unique<GradientBoostingModel>();
    model->bias = y.mean();
    model->learning_rate = spec.learning_rate;
    std::mt19937 rng(seed);
    VectorXd pred = VectorXd::Constant(y.size(), model->bias);
    const int max_features = spec.max_features > 0 ? spec.max_features : x.cols();
    for (int iter = 0; iter < spec.max_iter; ++iter) {
        const VectorXd residual = y - pred;
        auto tree = fit_tree(x, residual, std::min(spec.max_depth, 3), spec.min_samples_leaf, max_features, rng, false);
        pred.array() += model->learning_rate * tree->predict(x).array();
        model->trees.push_back(std::move(tree));
    }
    return model;
}

std::unique_ptr<RegressionModel> fit_regressor(const MatrixXd& x, const VectorXd& y, const LearnerSpec& spec, int seed) {
    if (y.size() == 0) {
        return std::make_unique<MeanModel>(0.0);
    }
    if (x.rows() < 5) {
        return std::make_unique<MeanModel>(y.mean());
    }
    switch (spec.kind) {
        case LearnerKind::kLasso:
            return fit_elastic_net_internal(x, y, 1.0, spec.cv_folds, spec.max_iter);
        case LearnerKind::kElasticNet:
        {
            double best_loss = std::numeric_limits<double>::infinity();
            std::unique_ptr<LinearModel> best_model;
            for (double ratio : spec.l1_ratios) {
                auto candidate = fit_elastic_net_internal(x, y, ratio, spec.cv_folds, spec.max_iter);
                const VectorXd err = y - candidate->predict(x);
                const double loss = err.array().square().mean();
                if (loss < best_loss) {
                    best_loss = loss;
                    best_model = std::move(candidate);
                }
            }
            return best_model;
        }
        case LearnerKind::kRandomForest:
            return fit_random_forest(x, y, spec, seed);
        case LearnerKind::kHistGradientBoosting:
            return fit_hist_gradient_boosting(x, y, spec, seed);
    }
    throw std::invalid_argument("Unknown learner kind");
}

struct StageResult {
    double tau_hat{0.0};
    double intercept_hat{0.0};
    double std_error{0.0};
    double ci_lower{0.0};
    double ci_upper{0.0};
    int iterations{0};
    bool converged{false};
    double objective{0.0};
    double kept_fraction{1.0};
};

VectorXd tukey_weights(const VectorXd& u, double c) {
    VectorXd w = VectorXd::Zero(u.size());
    for (Eigen::Index i = 0; i < u.size(); ++i) {
        const double abs_u = std::abs(u[i]);
        if (abs_u < c) {
            const double t = 1.0 - (u[i] * u[i]) / (c * c);
            w[i] = t * t;
        }
    }
    return w;
}

VectorXd tukey_psi(const VectorXd& u, double c) {
    VectorXd psi = VectorXd::Zero(u.size());
    for (Eigen::Index i = 0; i < u.size(); ++i) {
        const double abs_u = std::abs(u[i]);
        if (abs_u < c) {
            const double ratio = u[i] / c;
            const double t = 1.0 - ratio * ratio;
            psi[i] = u[i] * t * t;
        }
    }
    return psi;
}

VectorXd tukey_psi_prime(const VectorXd& u, double c) {
    VectorXd out = VectorXd::Zero(u.size());
    for (Eigen::Index i = 0; i < u.size(); ++i) {
        const double abs_u = std::abs(u[i]);
        if (abs_u < c) {
            const double ratio2 = (u[i] * u[i]) / (c * c);
            out[i] = (1.0 - ratio2) * (1.0 - 5.0 * ratio2);
        }
    }
    return out;
}

VectorXd weighted_least_squares(const MatrixXd& x, const VectorXd& y, const VectorXd& weights) {
    VectorXd sqrt_w = weights.array().max(kEps).sqrt();
    MatrixXd xw = x;
    for (Eigen::Index i = 0; i < x.rows(); ++i) {
        xw.row(i) *= sqrt_w[i];
    }
    VectorXd yw = y.array() * sqrt_w.array();
    return xw.completeOrthogonalDecomposition().solve(yw);
}

StageResult fit_mm_stage(const VectorXd& y_res, const VectorXd& d_res, double c, int max_iter, double tolerance, double ci_level) {
    MatrixXd x(y_res.size(), 2);
    x.col(0) = VectorXd::Ones(y_res.size());
    x.col(1) = d_res;
    VectorXd beta = weighted_least_squares(x, y_res, VectorXd::Ones(y_res.size()));
    const VectorXd pilot_res = y_res - x * beta;
    std::vector<std::pair<double, int>> order;
    order.reserve(static_cast<size_t>(pilot_res.size()));
    for (Eigen::Index i = 0; i < pilot_res.size(); ++i) {
        order.emplace_back(std::abs(pilot_res[i]), static_cast<int>(i));
    }
    std::sort(order.begin(), order.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
    const int keep_n = std::max(2, static_cast<int>(std::round(0.75 * pilot_res.size())));
    std::vector<int> keep_idx;
    keep_idx.reserve(static_cast<size_t>(keep_n));
    for (int i = 0; i < keep_n; ++i) {
        keep_idx.push_back(order[static_cast<size_t>(i)].second);
    }
    beta = weighted_least_squares(select_rows(x, keep_idx), select_rows(y_res, keep_idx), VectorXd::Ones(keep_idx.size()));

    bool converged = false;
    VectorXd weights = VectorXd::Ones(y_res.size());
    int iterations = 0;
    for (iterations = 1; iterations <= max_iter; ++iterations) {
        const VectorXd residuals = y_res - x * beta;
        const double scale = mad_scale(residuals);
        const VectorXd u = residuals / scale;
        weights = tukey_weights(u, c);
        const VectorXd beta_next = weighted_least_squares(x, y_res, weights);
        if ((beta_next - beta).cwiseAbs().maxCoeff() < tolerance) {
            beta = beta_next;
            converged = true;
            break;
        }
        beta = beta_next;
    }

    const VectorXd residuals = y_res - x * beta;
    const double scale = mad_scale(residuals);
    const VectorXd u = residuals / scale;
    const VectorXd psi = tukey_psi(u, c);
    const VectorXd psi_prime = tukey_psi_prime(u, c);
    const VectorXd scaled_prime = psi_prime / std::max(scale, kEps);
    const MatrixXd a = x.transpose() * scaled_prime.asDiagonal() * x / static_cast<double>(x.rows());
    const MatrixXd b = x.transpose() * psi.array().square().matrix().asDiagonal() * x / static_cast<double>(x.rows());
    const MatrixXd a_inv = a.completeOrthogonalDecomposition().pseudoInverse();
    const MatrixXd cov = a_inv * b * a_inv.transpose() / static_cast<double>(x.rows());
    const double se = std::sqrt(std::max(cov(1, 1), kEps));
    const double z = normal_quantile(0.5 + ci_level / 2.0);
    return StageResult{
        beta[1],
        beta[0],
        se,
        beta[1] - z * se,
        beta[1] + z * se,
        iterations,
        converged,
        (weights.array() * residuals.array().square()).mean(),
        (weights.array() > 1e-8).cast<double>().mean(),
    };
}

}  // namespace

LearnerKind learner_kind_from_name(const std::string& name) {
    if (name == "lasso") {
        return LearnerKind::kLasso;
    }
    if (name == "elastic_net") {
        return LearnerKind::kElasticNet;
    }
    if (name == "random_forest") {
        return LearnerKind::kRandomForest;
    }
    if (name == "hist_gradient_boosting") {
        return LearnerKind::kHistGradientBoosting;
    }
    throw std::invalid_argument("Unknown learner_name=" + name);
}

MmDmlResult fit_mm_dml(const FitInput& input, const DmlConfig& config) {
    if (input.x.rows() != input.d.size() || input.x.rows() != input.y.size()) {
        throw std::invalid_argument("X, D, and Y must have the same number of rows");
    }
    const auto folds = make_folds(static_cast<int>(input.x.rows()), std::max(2, config.folds), config.seed);
    VectorXd y_hat = VectorXd::Zero(input.y.size());
    VectorXd d_hat = VectorXd::Zero(input.d.size());

    for (size_t fold_idx = 0; fold_idx < folds.size(); ++fold_idx) {
        const auto& test_idx = folds[fold_idx];
        std::vector<char> mask(static_cast<size_t>(input.x.rows()), 1);
        for (int idx : test_idx) {
            mask[static_cast<size_t>(idx)] = 0;
        }
        std::vector<int> train_idx;
        train_idx.reserve(static_cast<size_t>(input.x.rows()) - test_idx.size());
        for (int row = 0; row < input.x.rows(); ++row) {
            if (mask[static_cast<size_t>(row)]) {
                train_idx.push_back(row);
            }
        }
        const MatrixXd x_train = select_rows(input.x, train_idx);
        const VectorXd y_train = select_rows(input.y, train_idx);
        const VectorXd d_train = select_rows(input.d, train_idx);
        const MatrixXd x_test = select_rows(input.x, test_idx);

        auto y_model = fit_regressor(x_train, y_train, config.learner, config.seed + static_cast<int>(fold_idx) * 17 + 1);
        auto d_model = fit_regressor(x_train, d_train, config.learner, config.seed + static_cast<int>(fold_idx) * 17 + 2);
        const VectorXd y_pred = y_model->predict(x_test);
        const VectorXd d_pred = d_model->predict(x_test);
        for (size_t i = 0; i < test_idx.size(); ++i) {
            const int row = test_idx[i];
            y_hat[row] = y_pred[static_cast<Eigen::Index>(i)];
            d_hat[row] = d_pred[static_cast<Eigen::Index>(i)];
        }
    }

    const VectorXd y_res = input.y - y_hat;
    const VectorXd d_res = input.d - d_hat;
    const StageResult stage = fit_mm_stage(y_res, d_res, config.mm_c, config.max_iter, config.tolerance, config.ci_level);
    const VectorXd stage_residuals = y_res.array() - stage.intercept_hat - stage.tau_hat * d_res.array();
    return MmDmlResult{
        stage.tau_hat,
        stage.intercept_hat,
        stage.std_error,
        stage.ci_lower,
        stage.ci_upper,
        stage.iterations,
        stage.converged,
        stage.objective,
        stage.kept_fraction,
        standard_deviation(y_res),
        standard_deviation(d_res),
        std::vector<double>(y_res.data(), y_res.data() + y_res.size()),
        std::vector<double>(d_res.data(), d_res.data() + d_res.size()),
        std::vector<double>(stage_residuals.data(), stage_residuals.data() + stage_residuals.size()),
    };
}

}  // namespace mm_dml
