/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 13/07/2026
 */

// Conjugate-gradient method to solve the Riemannian Newton equation:

void tcg(arguments_optim& x,
         std::vector<transformations*>& xtransforms,
         std::vector<manifolds*>& xmanifolds,
         std::vector<estimators*>& xestimators,
         bool& att_bnd, arma::vec c, double rad) {

  /*
   * Truncated conjugate gradient sub-solver for the trust-region sub-problem
   * From Liu (Algorithm 4; 2020)
   */

  product_manifold* final_manifold;
  product_transform* final_transform;
  product_estimator* final_estimator;

  x.dir.zeros();
  arma::vec dir0;

  double alpha, rr0, tau, beta, dHd;
  x.dparameters = -x.rg; // Initial search direction
  arma::vec r = x.dparameters; // Initial residual
  double rr = x.ng * x.ng;
  double tol = x.ng * std::min(pow(x.ng, c[0]), c[1]);

  int iter = 0;

  // In x, g should be already computed
  // final_estimator->param(x, xestimators); // Unnecessary
  // final_estimator->G(x, xestimators); // Unnecessary

  // final_manifold->param(x, xmanifolds); // Unnecessary
  // final_estimator->param(x, xestimators); // Unnecessary
  final_transform->dtransform(x, xtransforms);
  final_estimator->dG(x, xestimators);
  final_transform->update_dgrad(x, xtransforms);
  final_manifold->param(x, xmanifolds);
  final_manifold->hess(x, xmanifolds);

  do {

    dHd = arma::accu(x.dparameters % x.dH);

    if(dHd <= 0) {

      tau = root_quad(arma::accu(x.dparameters % x.dparameters),
                      2 * arma::accu(x.dir % x.dparameters),
                      arma::accu(x.dir % x.dir) - rad * rad); // Solve equation 39
      x.dir += tau * x.dparameters;
      att_bnd = true;

      break;

    }

    rr0 = rr;
    alpha = rr0 / dHd;
    dir0 = x.dir;
    x.dir += alpha * x.dparameters; // update proposal

    if (sqrt(arma::accu(x.dir % x.dir)) >= rad) {

      tau = root_quad(arma::accu(x.dparameters % x.dparameters),
                      2 * arma::accu(dir0 % x.dparameters),
                      arma::accu(dir0 % dir0) - rad * rad); // Solve equation 39
      x.dir = dir0 + tau * x.dparameters;
      att_bnd = true;

      break;

    }

    r = r - alpha * x.dH; // update gradient
    rr = arma::accu(r % r);

    if (sqrt(rr) < tol) {

      att_bnd = false;
      break;

    }

    beta = rr / rr0;
    x.dparameters = r + beta * x.dparameters;
    iter += 1;

    final_transform->dtransform(x, xtransforms);
    final_estimator->dG(x, xestimators);
    final_transform->update_dgrad(x, xtransforms);
    final_manifold->param(x, xmanifolds);
    final_manifold->hess(x, xmanifolds);

  } while (iter < x.tcg_maxit);

  x.dparameters = x.dir;

}
