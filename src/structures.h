/*
 * Author: Marcos Jiménez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 29/06/2025
 */

#include <tuple>

struct arguments_optim{

  // Manifolds and Estimators:
  int nmanifolds, nestimators, ntransforms;

  // Estimators stuff:
  std::vector<arma::uvec> indexes, target_indexes;
  double f = 0.00, f_null = 0.00, loglik = 0.00;

  // Optim stuff:
  double c1 = 10e-04, c2 = 0.5, eps = 1e-05, ss_fac = 2, ss_min = 0.1,
    ng = 1, inprod = 1, step_eps = 1e-09, df = 1000, df_eps = 1e-09;
  int M = 15L, step_maxit = 30L, iterations = 0L, maxit = 10000L, tcg_maxit = 5,
    rstarts = 1L, cores = 1L;
  double old_inprod = arma::datum::inf;
  int step_iteration = 0L;
  bool print = false;

  double ss = 1;
  arma::vec dir = {1};

  std::string search = "back";
  bool convergence = false;
  arma::vec parameters, dparameters, g, dg, rg, drg, dH;
  arma::vec transparameters, dtransparameters, grad, dgrad;
  arma::mat jacob, h, hess, B;
  arma::mat posterior;
  arma::mat latentloglik;
  arma::vec latentpars, loglatentpars;
  arma::vec weights;
  arma::vec se;
  std::vector<arma::mat> modhessian, dparam_dS;
  std::vector<arma::vec> classes; // P(X = c) // classes_hat
  std::vector<std::vector<arma::mat>> conditionals; // conditionals_hat

  // Checks:
  // Rcpp::Nullable<Rcpp::List> nullable_control = R_NilValue;
  std::string optimizer = "newton", std_error = "normal";
  arma::uvec lower, upper;

  // Manifolds:
  arma::mat X, dX, dL, dP, Phi, A, Phi_Target;
  arma::uvec oblq_indexes;

  // Output:
  // Rcpp::List lambda, phi, psi, Rhat, residuals, R;
  std::vector<double> fs;
  // int df = 0L, df_null = 0L, total_nobs = 0L;
  std::vector<int> nobs, p, q;
  // Rcpp::CharacterVector cor, estimator, projection;

  // Outcomes:
  // std::vector<std::vector<double>> doubles;
  // std::vector<std::vector<arma::vec>> vectors;
  // std::vector<std::vector<arma::mat>> matrices;
  // std::vector<std::vector<std::vector<arma::mat>>> list_matrices;

  // Manifolds:
  std::tuple<
    std::vector<std::vector<double>>,
    std::vector<std::vector<arma::vec>>,
    std::vector<std::vector<arma::mat>>,
    std::vector<std::vector<arma::cube>>,
    std::vector<std::vector<std::vector<arma::vec>>>,
    std::vector<std::vector<std::vector<arma::mat>>>
  > outputs_manifold;

  // Transformations:
  std::tuple<
    std::vector<std::vector<double>>,
    std::vector<std::vector<arma::vec>>,
    std::vector<std::vector<arma::mat>>,
    std::vector<std::vector<arma::cube>>,
    std::vector<std::vector<std::vector<arma::vec>>>,
    std::vector<std::vector<std::vector<arma::mat>>>
  > outputs_transform;

  // Estimators:
  std::tuple<
    std::vector<std::vector<double>>,
    std::vector<std::vector<arma::vec>>,
    std::vector<std::vector<arma::mat>>,
    std::vector<std::vector<arma::cube>>,
    std::vector<std::vector<std::vector<arma::vec>>>,
    std::vector<std::vector<std::vector<arma::mat>>>
  > outputs_estimator;

  int nparam, ntransparam, nrow_post, ncol_post;
  arma::uvec param2transparam, transparam2param;

};
