/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 11/09/2026
 */

// Objective-only backtracking line search satisfying the Armijo condition:

void armijo(arguments_optim& x,
            std::vector<transformations*>& xtransforms,
            std::vector<manifolds*>& xmanifolds,
            std::vector<estimators*>& xestimators) {

  product_manifold final_manifold;
  product_transform final_transform;
  product_estimator final_estimator;

  double f0 = x.f;
  arma::vec parameters = x.parameters;
  x.inprod = arma::dot(x.dir, x.rg);
  x.step_iteration = 0L;

  // Initial step size:
  x.ss = std::max(x.ss_min, x.ss*x.ss_fac);
  double max_dir = arma::abs(x.dir).max();
  if(max_dir > 1.0) {
    x.ss = std::min(x.ss, 1.0/max_dir);
  }

  do {

    ++x.step_iteration;
    x.parameters = parameters + x.ss*x.dir;

    // Projection onto the manifold:
    final_manifold.param(x, xmanifolds);
    final_manifold.retr(x, xmanifolds);

    // Parameterization:
    final_transform.transform(x, xtransforms);
    final_estimator.param(x, xestimators);
    final_estimator.F(x, xestimators);

    x.df = x.f - f0;

    // Accept the step if the Armijo condition is satisfied:
    if(x.parameters.is_finite() && std::isfinite(x.f) &&
       x.df <= x.c1*x.ss*x.inprod) {
      break;
    }

    x.ss *= x.c2;

  } while(x.step_iteration < x.step_maxit);

  // Do not return a nonfinite final trial:
  if(!x.parameters.is_finite() || !std::isfinite(x.f)) {
    x.parameters = parameters;
    final_manifold.param(x, xmanifolds);
    final_transform.transform(x, xtransforms);
    final_estimator.param(x, xestimators);
    final_estimator.F(x, xestimators);
    x.df = 0.0;
  }

}
