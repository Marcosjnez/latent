/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 24/08/2026
 */

#include <tuple>

struct arguments_optim{

  // Manifolds and Estimators:
  int nmanifolds, nestimators, ntransforms;
  int ncov_transform = 1L, ncov_params;
  arma::uvec vector_nparams;

  // Estimators stuff:
  std::vector<arma::uvec> indexes, target_indexes;
  double f = 0.00, f_null = 0.00, loglik = 0.00;

  // Optim stuff:
  double c1 = 10e-04, c2 = 0.5, eps = 1e-05, ss_fac = 2, ss_min = 0.1,
    ng = 1, max_rg = 1, inprod = 1, step_eps = 1e-09, df = 1000, df_eps = 1e-09;
  int M = 15L, step_maxit = 30L, iterations = 0L, maxit = 10000L, tcg_maxit = 5,
    rstarts = 1L, cores = 1L;
  double old_inprod = arma::datum::inf;
  int step_iteration = 0L;
  bool print = false;
  int print_interval = 50L;

  double ss = 1;
  arma::vec dir = {1};

  std::string search = "back";
  bool convergence = false;
  arma::vec parameters, dparameters, g, dg, rg, drg, dH;
  arma::vec transparameters, transparameters_init, dtransparameters, grad,
  dgrad, grad_init, dgrad_init;
  arma::mat h, v, B, inv_h, inv_hess;
  arma::sp_mat jacob, vcov, dconstr, d2constraints;
  arma::mat hess; // Returned as dgCMatrix class type from Matrix package
  arma::mat posterior;
  arma::mat freqs;
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

  // Manifolds:
  std::tuple<
    std::vector<std::vector<double>>,
    std::vector<std::vector<arma::vec>>,
    std::vector<std::vector<arma::mat>>,
    std::vector<std::vector<arma::cube>>,
    std::vector<std::vector<std::vector<arma::vec>>>,
    std::vector<std::vector<std::vector<arma::mat>>>,
    std::vector<std::vector<std::string>>,
    std::vector<std::vector<std::string>>,
    std::vector<std::vector<std::string>>,
    std::vector<std::vector<std::string>>,
    std::vector<std::vector<std::string>>,
    std::vector<std::vector<std::string>>
  > outputs_manifold;

  // Transformations:
  std::tuple<
    std::vector<std::vector<double>>,
    std::vector<std::vector<arma::vec>>,
    std::vector<std::vector<arma::mat>>,
    std::vector<std::vector<arma::cube>>,
    std::vector<std::vector<std::vector<arma::vec>>>,
    std::vector<std::vector<std::vector<arma::mat>>>,
    std::vector<std::vector<std::string>>,
    std::vector<std::vector<std::string>>,
    std::vector<std::vector<std::string>>,
    std::vector<std::vector<std::string>>,
    std::vector<std::vector<std::string>>,
    std::vector<std::vector<std::string>>
  > outputs_transform;

  // Estimators:
  std::tuple<
    std::vector<std::vector<double>>,
    std::vector<std::vector<arma::vec>>,
    std::vector<std::vector<arma::mat>>,
    std::vector<std::vector<arma::cube>>,
    std::vector<std::vector<std::vector<arma::vec>>>,
    std::vector<std::vector<std::vector<arma::mat>>>,
    std::vector<std::vector<std::string>>,
    std::vector<std::vector<std::string>>,
    std::vector<std::vector<std::string>>,
    std::vector<std::vector<std::string>>,
    std::vector<std::vector<std::string>>,
    std::vector<std::vector<std::string>>
  > outputs_estimator;

  int pick = 0L;

  int nparam, ntransparam, nrow_post, ncol_post;
  arma::uvec transparam2param;
  arma::uvec idx_transforms;
  int mstep_maxit = 100L;
  double mstep_eps = 1e-06;
  std::string mopt = "lbfgs";

};

inline void append_constraint_derivatives(
    arguments_optim& x,
    const arma::sp_mat& first_derivatives,
    const std::vector<arma::sp_mat>& second_derivatives) {

  const arma::uword ntrans = x.transparameters.n_elem;
  const arma::uword nnew = first_derivatives.n_cols;

  if(first_derivatives.n_rows != ntrans) {
    Rcpp::stop("Constraint first derivatives must have one row per transformed parameter.");
  }

  if(second_derivatives.size() != nnew) {
    Rcpp::stop("There must be one constraint Hessian for every constraint derivative column.");
  }

  if(x.dconstr.n_rows != ntrans) {
    Rcpp::stop("The incremental constraint derivative matrix has an invalid number of rows.");
  }

  const arma::uword nold = x.dconstr.n_cols;
  const arma::uword expected_old_dimension = nold*ntrans;

  if(x.d2constraints.n_rows != expected_old_dimension ||
     x.d2constraints.n_cols != expected_old_dimension) {
    Rcpp::stop("The incremental second-constraint derivative matrix is not aligned with dconstr.");
  }

  x.dconstr.resize(ntrans, nold+nnew);

  for(arma::sp_mat::const_iterator it = first_derivatives.begin();
      it != first_derivatives.end(); ++it) {
    x.dconstr(it.row(), nold+it.col()) = *it;
  }

  const arma::uword new_dimension = (nold+nnew)*ntrans;
  x.d2constraints.resize(new_dimension, new_dimension);

  for(arma::uword j=0L; j < nnew; ++j) {

    const arma::sp_mat& H = second_derivatives[j];

    if(H.n_rows != ntrans || H.n_cols != ntrans) {
      Rcpp::stop("Every constraint Hessian must be square with one row and column per transformed parameter.");
    }

    const arma::uword offset = (nold+j)*ntrans;

    for(arma::sp_mat::const_iterator it = H.begin();
        it != H.end(); ++it) {
      x.d2constraints(offset+it.row(),
                      offset+it.col()) = *it;
    }

  }

}
