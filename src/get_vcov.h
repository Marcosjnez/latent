/*
 * Author: Marcos Jiménez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 27/12/2025
 */

Rcpp::List get_vcov(Rcpp::List control_manifold,
                    Rcpp::List control_transform,
                    Rcpp::List control_estimator,
                    Rcpp::List control_optimizer,
                    arma::mat H,
                    int cores) {

  Rcpp::List result;
  arguments_optim x;

  x.nmanifolds = control_manifold.size();
  x.ntransforms = control_transform.size();
  x.nestimators = control_estimator.size();

  product_manifold* final_manifold;
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
  // final_transform->jacobian(x, xtransforms);

  x.h = H;
  if(x.h.is_empty()) {

    int npar = x.parameters.n_elem;
    x.h.set_size(npar, npar);

    for(int i=0; i < npar; ++i) {

      x.dparameters.zeros();
      x.dparameters(i) = 1.00;

      final_transform->dtransform(x, xtransforms);
      final_estimator->dG(x, xestimators);
      final_transform->update_dgrad(x, xtransforms);
      x.h.col(i) = x.dg;

    }

  }

  x.h = H;

//   if(x.h.is_empty()) {
//
//     int npar = x.parameters.n_elem;
//
//     arma::mat h(npar, npar, arma::fill::none);
//
// #pragma omp parallel num_threads(cores)
// {
//   // x.h is still empty here, so the thread-local copy does not
//   // duplicate the npar × npar Hessian matrix
//   arguments_optim x_local = x;
//
// #pragma omp for schedule(static)
//   for(int i=0; i < npar; ++i) {
//
//     x_local.dparameters.zeros();
//     x_local.dparameters(i) = 1.00;
//
//     final_transform->dtransform(x_local, xtransforms);
//     final_estimator->dG(x_local, xestimators);
//     final_transform->update_dgrad(x_local, xtransforms);
//
//     // Each iteration writes to a different column
//     std::copy_n(x_local.dg.memptr(), x_local.dg.n_elem,
//                 h.colptr(i));
//
//   }
// }
//
// x.h = std::move(h);
//
//   }


  final_transform->update_vcov(x, xtransforms);

  result["vcov"] = x.vcov;
  result["se"] = x.se;

  return result;

}
