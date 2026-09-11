/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 11/09/2026
 */

#ifndef LATENT_NEWTON_TCG_H
#define LATENT_NEWTON_TCG_H

// Liu's CG subsolver: a bounded trust-region step or an unbounded Newton direction.

void tcg(arguments_optim& x,
         std::vector<transformations*>& xtransforms,
         std::vector<manifolds*>& xmanifolds,
         std::vector<estimators*>& xestimators,
         bool& att_bnd, const arma::vec& c, double rad,
         bool constrained = true) {

  /*
   * Truncated conjugate gradient subsolver for both Newton variants
   * From Liu (Algorithm 4; 2020)
   */

  product_manifold final_manifold;
  product_transform final_transform;
  product_estimator final_estimator;

  x.dir.zeros();
  att_bnd = false;
  arma::vec dir0;

  double alpha, rr0, tau, beta, dHd;
  x.dparameters = -x.rg; // Initial search direction
  arma::vec r = x.dparameters; // Initial residual
  double rr = x.ng * x.ng;
  double tol = x.ng * std::min(pow(x.ng, c[0]), c[1]);

  int iter = 0;

  // Avoid the zero-gradient boundary equation (0/0) at a stationary point.
  if(x.ng == 0.0) {
    x.dparameters = x.dir;
    return;
  }

  // The gradient is already available at the current point.

  final_transform.dtransform(x, xtransforms);
  final_estimator.dG(x, xestimators);
  final_transform.update_dgrad(x, xtransforms);
  final_manifold.param(x, xmanifolds);
  final_manifold.hess(x, xmanifolds);

  do {

    dHd = arma::accu(x.dparameters % x.dH);

    if(!std::isfinite(dHd)) {
      if(!constrained && iter == 0L) x.dir = -x.rg;
      break;
    }

    if(dHd <= 0) {

      if(constrained) {

        tau = root_quad(arma::accu(x.dparameters % x.dparameters),
                        2 * arma::accu(x.dir % x.dparameters),
                        arma::accu(x.dir % x.dir) - rad * rad);
        x.dir += tau * x.dparameters;
        att_bnd = true;

      } else if(iter == 0L) {

        // Liu's unbounded negative-curvature fallback.
        x.dir = x.dparameters;

      }

      break;

    }

    rr0 = rr;
    alpha = rr0 / dHd;
    dir0 = x.dir;
    x.dir += alpha * x.dparameters; // update proposal

    if(constrained && sqrt(arma::accu(x.dir % x.dir)) >= rad) {

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

    final_transform.dtransform(x, xtransforms);
    final_estimator.dG(x, xestimators);
    final_transform.update_dgrad(x, xtransforms);
    final_manifold.param(x, xmanifolds);
    final_manifold.hess(x, xmanifolds);

  } while (iter < x.tcg_maxit);

  x.dparameters = x.dir;

}

// Unbounded CG for Newton with an Armijo or Wolfe line search.
void tcg(arguments_optim& x,
         std::vector<transformations*>& xtransforms,
         std::vector<manifolds*>& xmanifolds,
         std::vector<estimators*>& xestimators) {

  bool att_bnd = false;
  const arma::vec c = {1.0, 0.01};
  tcg(x, xtransforms, xmanifolds, xestimators, att_bnd, c, 0.0, false);

}

#endif
