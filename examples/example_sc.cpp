#include "robustcause/robustcause.hpp"

#include <armadillo>

#include <iostream>

int main() {
  arma::mat outcomes(7, 4);
  outcomes.col(0) = arma::vec{1.0, 2.1, 9.0, 4.0, 5.1, 6.0, 7.0};
  outcomes.col(1) = arma::vec{1.1, 2.0, 3.0, 4.1, 5.0, 6.1, 7.2};
  outcomes.col(2) = arma::vec{0.8, 1.9, 3.2, 3.9, 5.2, 5.9, 6.9};
  outcomes.col(3) = arma::vec{1.3, 2.2, 2.8, 4.2, 4.9, 6.0, 7.1};

  robust::ScPanel panel;
  panel.outcomes = outcomes;
  panel.treated_index = 0;
  panel.treatment_start = 5;
  panel.predictor_treated = arma::vec{2.0, 4.0};
  panel.predictor_donors.set_size(2, 3);
  panel.predictor_donors.row(0) = arma::rowvec{2.0, 1.9, 2.1};
  panel.predictor_donors.row(1) = arma::rowvec{4.1, 3.9, 4.0};
  panel.predictor_weights = arma::vec{5.0, 5.0};

  const robust::ScData data = robust::prepare_sc_data(panel);
  const robust::ScResult sc_fit = robust::fit_sc(data);
  const robust::MmscResult mm_fit = robust::fit_mm_sc(data);
  const robust::ScPlaceboResult placebo = robust::fit_mm_sc_placebos(outcomes, 5);

  std::cout << "Standard SC weights:\n" << sc_fit.weights << "\n";
  std::cout << "MM-SC weights:\n" << mm_fit.weights << "\n";
  std::cout << "Robust time weights:\n" << mm_fit.robust_time_weights << "\n";
  std::cout << "Post-treatment gaps:\n" << mm_fit.post_gaps << "\n";
  std::cout << "Placebo RMSPE ratios:\n" << placebo.rmspe_ratio << "\n";
  return 0;
}
