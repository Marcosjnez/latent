/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 13/07/2026
 */

// Line-search algorithm satisfying the armijo condition:

void armijo(arguments_optim& x,
            product_manifold* final_manifold, product_estimator* final_estimator,
            std::vector<manifolds*>& xmanifolds, std::vector<estimators*>& xestimators) {

  x.ss = std::max(x.ss_min, x.ss * x.ss_fac);
  // x.ss = x.ss*2;
  double f0 = x.f;
  arma::vec parameters = x.parameters;
  x.inprod = arma::accu(x.dir % x.rg);
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
    final_estimator->param(x, xestimators);
    final_estimator->F(x, xestimators);
    x.df = x.f - f0;
    // if (x.df < x.c1 * x.ss * x.inprod) break;
    if (x.df < x.c1 * x.ss * x.inprod || // armijo condition
        x.ss < x.step_eps) break;
    // if (x.df < 0) break;
    x.ss *= x.c2;

  } while (x.step_iteration <= x.step_maxit);

  if (x.ss < std::numeric_limits<double>::epsilon()) {
    x.ss = std::numeric_limits<double>::epsilon();
  }
  // Rcpp::Rcout << "Armijo =" << iteration << std::endl;

}
