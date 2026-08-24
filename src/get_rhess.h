/*
 * Author: Marcos Jiménez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 24/08/2026
 */

Rcpp::List get_rhess(Rcpp::List control_manifold,
                     Rcpp::List control_transform,
                     Rcpp::List control_estimator,
                     Rcpp::List control_optimizer,
                     int cores) {

  Rcpp::List result;
  arguments_optim x;

  x.nmanifolds = control_manifold.size();
  x.ntransforms = control_transform.size();
  x.nestimators = control_estimator.size();

  product_manifold final_manifold;
  product_transform final_transform;
  product_estimator final_estimator;

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

  final_transform.transform(x, xtransforms);
  final_estimator.param(x, xestimators);
  final_estimator.G(x, xestimators);
  final_transform.update_grad(x, xtransforms);
  final_manifold.param(x, xmanifolds);
  final_manifold.proj(x, xmanifolds);
  final_manifold.tangent_basis(x, xmanifolds);

  arma::mat T = final_manifold.T;

  if(T.n_rows != static_cast<arma::uword>(x.nparam)) {
    Rcpp::stop("The product tangent-space basis must have one row per free parameter.");
  }

  if(!T.is_finite()) {
    Rcpp::stop("The product tangent-space basis contains non-finite values.");
  }

  if(T.n_cols > 0L) {

    arma::mat gram = T.t()*T;
    arma::mat identity = arma::eye(T.n_cols, T.n_cols);

    if(!arma::approx_equal(gram, identity, "absdiff", 1e-08)) {
      Rcpp::stop("The product tangent-space basis is not orthonormal.");
    }

  }

  arma::uword ntangent = T.n_cols;
  arma::mat h(ntangent, ntangent, arma::fill::zeros);

  for(arma::uword i=0L; i < ntangent; ++i) {

    x.dparameters = T.col(i);

    final_transform.dtransform(x, xtransforms);
    final_estimator.dG(x, xestimators);
    final_transform.update_dgrad(x, xtransforms);
    final_manifold.hess(x, xmanifolds);

    h.col(i) = T.t()*x.dH;

  }

  result["h"] = h;
  result["T"] = T;

  (void)algorithm;
  (void)cores;

  return result;

}
