/*
 * Author: Marcos Jiménez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 03/02/2025
 */

arma::vec grad_comp(arma::vec parameters,
                    Rcpp::List control_manifold,
                    Rcpp::List control_estimator,
                    Rcpp::List control_optimizer) {

  /*
   * parameters: vector of parameters
   */

  /* control_manifold: list of manifolds
   * Each list projects a set of parameters onto a manifold
   * NOTE: each list must contain nonoverlapping parameters
   * manifold_indices: list of indices that relate sets of parameters to manifolds
   */

  /* control_estimator: list of estimators
   * estimator_indices: list of indices that relate sets of parameters to estimators
   * estimator_target: list of indices that relate estimators to parameters
   */

  arguments_optim x;

  x.nmanifolds = control_manifold.size();
  x.nestimators = control_estimator.size();

  product_manifold* final_manifold;
  product_estimator* final_estimator;

  std::vector<manifolds*> xmanifolds(x.nmanifolds);
  std::vector<estimators*> xestimators(x.nestimators);

  optim* algorithm = choose_optim(x, control_optimizer);
  x.parameters = parameters;

  for(int i=0; i < x.nmanifolds; ++i) {
    xmanifolds[i] = choose_manifold(control_manifold[i], xmanifolds[i]);
  }

  for(int i=0; i < x.nestimators; ++i) {
    xestimators[i] = choose_estimator(control_estimator[i], xestimators[i]);
  }

  /*
   * Computations
   */

  Rcpp::List computations;

  final_estimator->param(x, xestimators);
  final_estimator->F(x, xestimators);
  final_estimator->G(x, xestimators);

  return x.g;

}
