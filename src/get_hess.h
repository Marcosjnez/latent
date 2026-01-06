/*
 * Author: Marcos Jiménez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 05/10/2025
 */

Rcpp::List get_hess(Rcpp::List control_manifold,
                    Rcpp::List control_transform,
                    Rcpp::List control_estimator,
                    Rcpp::List control_optimizer) {

  Rcpp::List result;
  arguments_optim x;

  x.nmanifolds = control_manifold.size();
  x.ntransforms = control_transform.size();
  x.nestimators = control_estimator.size();

  product_transform* final_transform;
  product_estimator* final_estimator;

  std::vector<manifolds*> xmanifolds(x.nmanifolds);
  std::vector<transformations*> xtransforms(x.ntransforms);
  std::vector<estimators*> xestimators(x.nestimators);

  optim* algorithm = choose_optim(x, control_optimizer);

  std::vector<arma::vec> parameters_list = control_optimizer["parameters"];
  std::vector<arma::vec> transparameters_list = control_optimizer["transparameters"];
  arma::vec parameters = parameters_list[0];
  arma::vec transparameters = transparameters_list[0];
  x.parameters = parameters;
  x.transparameters = transparameters;

  for(int i=0; i < x.nmanifolds; ++i) {
    xmanifolds[i] = choose_manifold(control_manifold[i]);
  }

  for(int i=0; i < x.ntransforms; ++i) {
    xtransforms[i] = choose_transform(control_transform[i]);
  }

  for(int i=0; i < x.nestimators; ++i) {
    xestimators[i] = choose_estimator(control_estimator[i]);
  }

  /*
   * Computations
   */

  final_transform->transform(x, xtransforms);
  final_estimator->param(x, xestimators);
  final_estimator->G(x, xestimators);
  final_transform->update_grad(x, xtransforms);

  int npar = x.parameters.n_elem;
  arma::mat h(npar, npar, arma::fill::zeros);

  for(int i=0; i<npar; ++i) {

    x.dparameters.zeros();
    x.dparameters(i) = 1.00;

    final_transform->dparam(x, xtransforms);
    final_estimator->dG(x, xestimators);
    final_transform->update_dgrad(x, xtransforms);
    h.col(i) = x.dg;

  }

  result["h"] = h;

  return result;

}
