/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 13/07/2026
 */

// Line-search algorithm satisfying the Wolfe conditions:

void wolfe(arguments_optim& x,
           std::vector<transformations*>& xtransforms,
           std::vector<manifolds*>& xmanifolds,
           std::vector<estimators*>& xestimators) {

  product_manifold* final_manifold;
  product_transform* final_transform;
  product_estimator* final_estimator;

  // x.ss = std::max(x.ss_min, x.ss * x.ss_fac);
  // x.ss = x.ss*2;
  double f0 = x.f;
  arma::vec parameters = x.parameters;
  x.inprod = arma::dot(x.dir, x.rg);
  // x.inprod = arma::dot(-x.dir, x.rg);
  // x.inprod = arma::accu(x.dir % x.dir);

  x.step_iteration = 0L;

  do {

    ++x.step_iteration;
    x.parameters = parameters + x.ss*x.dir;

    // Projection onto the manifold
    final_manifold->param(x, xmanifolds);
    final_manifold->retr(x, xmanifolds);
    // final_manifold->param(x, xmanifolds); // Unnecessary
    // Parameterization
    final_transform->transform(x, xtransforms);
    final_estimator->param(x, xestimators);
    final_estimator->F(x, xestimators);
    final_estimator->G(x, xestimators);
    final_transform->update_grad(x, xtransforms);
    final_manifold->param(x, xmanifolds);
    final_manifold->proj(x, xmanifolds);
    x.df = x.f - f0;
    double inprod = arma::dot(x.dir, x.rg); // Armijo condition
    if (x.df > x.c1 * x.ss * x.inprod ||
        x.ss < x.step_eps) {
      x.ss *= x.c2;
    } else if (inprod < x.c2 * x.inprod) { // Curvature condition
      x.ss *= x.ss_fac;
    } else {
      break;
    }

  } while (x.step_iteration < x.step_maxit);

  if (x.ss < std::numeric_limits<double>::epsilon()) {
    x.ss = std::numeric_limits<double>::epsilon();
  }
  // Rcpp::Rcout << "Armijo =" << iteration << std::endl;

}
