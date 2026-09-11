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

  const double f0 = x.f;
  const arma::vec parameters = x.parameters;
  x.inprod = arma::dot(x.dir, x.rg);
  x.step_iteration = 0L;

  if(!std::isfinite(f0) || !x.dir.is_finite() ||
     !x.rg.is_finite() || !std::isfinite(x.inprod)) {
    Rf_error("Armijo requires a finite objective, gradient, and direction.");
  }

  // An exactly stationary starting point does not require a trial step.
  if(arma::dot(x.rg, x.rg) == 0.0) {
    x.df = 0.0;
    x.ss = 0.0;
    return;
  }

  if(x.inprod >= 0.0) {
    Rf_error("Armijo requires a descent direction.");
  }

  // Retain the existing Armijo initial-step convention.
  x.ss = std::max(x.ss_min, x.ss * x.ss_fac);
  double max_dir = arma::abs(x.dir).max();
  if(max_dir > 1.0) {
    x.ss = std::min(x.ss, 1.0/max_dir);
  }

  if(!std::isfinite(x.ss) || x.ss <= 0.0) {
    Rf_error("Armijo requires a finite positive initial step size.");
  }

  for(int i=0L; i < x.step_maxit; ++i) {

    ++x.step_iteration;
    x.parameters = parameters + x.ss*x.dir;

    final_manifold.param(x, xmanifolds);
    final_manifold.retr(x, xmanifolds);
    final_manifold.param(x, xmanifolds);

    // Every trial must refresh the complete transformation chain.
    final_transform.transform(x, xtransforms);
    final_estimator.param(x, xestimators);
    final_estimator.F(x, xestimators);
    x.df = x.f - f0;

    if(x.parameters.is_finite() && std::isfinite(x.f) &&
       x.df <= x.c1*x.ss*x.inprod) {
      return;
    }

    if(x.ss <= x.step_eps || i+1L == x.step_maxit) {
      break;
    }

    x.ss *= x.c2;

    if(x.ss == 0.0) {
      break;
    }

  }

  // Restore both the parameters and the shared component caches on failure.
  // A failed trial must not be returned as an accepted optimization step.
  x.parameters = parameters;
  final_manifold.param(x, xmanifolds);
  final_transform.transform(x, xtransforms);
  final_estimator.param(x, xestimators);
  final_estimator.F(x, xestimators);
  x.df = 0.0;

  Rf_error("Armijo could not find an acceptable step. "
             "Consider increasing control$step_maxit or reducing control$step_eps.");

}
