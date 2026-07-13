/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 13/07/2026
 */

// Newton Trust-region algorithm:

optim_result ntr(arguments_optim x,
                 std::vector<transformations*>& xtransforms,
                 std::vector<manifolds*>& xmanifolds,
                 std::vector<estimators*>& xestimators) {

  product_manifold* final_manifold;
  product_transform* final_transform;
  product_estimator* final_estimator;

  // Ensure initial parameters are ok:
  final_manifold->param(x, xmanifolds);
  final_manifold->retr(x, xmanifolds);
  final_transform->transform(x, xtransforms);
  final_estimator->param(x, xestimators);

  // Objective
  final_estimator->F(x, xestimators);
  // Gradient
  final_estimator->G(x, xestimators);
  final_transform->update_grad(x, xtransforms);
  // Riemannian gradient
  final_manifold->proj(x, xmanifolds);

  /*
   * Riemannian newton trust-region algorithm
   * From Liu (Algorithm 2; 2020)
   */

  // Rcpp::Rcout << "x.f = " << x.f << std::endl;

  x.inprod = arma::accu(x.rg % x.rg);
  x.ng = std::sqrt(x.inprod);

  double max_rad = 10;

  arma::vec fac_rad(2);
  fac_rad[0] = 0.25;
  fac_rad[1] = 2;

  arma::vec crit_goa(3);
  crit_goa[0] = 0.20;
  crit_goa[1] = 0.25;
  crit_goa[2] = 0.75;

  arma::vec c(2);
  c[0] = 1;
  c[1] = 0.01;

  double rad = 1;
  bool att_bnd = false;

  x.iterations = 0;
  double goa, preddiff;

  x.convergence = false;

  do {

    // subsolver
    tcg(x, xtransforms, xmanifolds, xestimators, att_bnd, c, rad); // Update x.dparameters, x.dg, and x.dH

    final_transform->dtransform(x, xtransforms);
    final_estimator->dG(x, xestimators);
    final_transform->update_dgrad(x, xtransforms);
    final_manifold->param(x, xmanifolds);
    final_manifold->proj(x, xmanifolds); // REMOVE
    final_manifold->hess(x, xmanifolds);

    preddiff = - arma::accu(x.dparameters % ( x.rg + 0.5 * x.dH) );

    arguments_optim new_x = x;
    new_x.parameters += x.dir;
    // Projection onto the manifold
    final_manifold->param(new_x, xmanifolds);
    final_manifold->retr(new_x, xmanifolds);
    final_manifold->param(new_x, xmanifolds);
    // Parameterization
    final_transform->transform(new_x, xtransforms);
    final_estimator->param(new_x, xestimators);
    final_estimator->F(new_x, xestimators);

    x.df = x.f - new_x.f; // New stopping criteria

    if ( std::abs(preddiff) <= arma::datum::eps ) {

      goa = 1;

    } else {

      goa = x.df / preddiff;

    }

    if (goa < crit_goa[1]) {

      rad = fac_rad[0] * rad;

    } else if (goa > crit_goa[2] && att_bnd) {

      rad = std::min(fac_rad[1] * rad, max_rad);

    }

    // accepted iteration
    if (goa > crit_goa[0]) {

      x = new_x;

      // update gradient
      final_estimator->G(x, xestimators);
      final_transform->update_grad(x, xtransforms);
      // Riemannian gradient
      final_manifold->proj(x, xmanifolds);

      // x.inprod = arma::accu(arma::abs(x.dir % x.rg));
      // x.inprod = arma::accu(x.dir % x.dir);
      x.inprod = arma::accu(x.rg % x.rg);
      x.max_rg = arma::abs(x.rg).max();

      x.ng = std::sqrt(x.inprod);

    }

    ++x.iterations;

    // Print info:
    if((x.print & x.rstarts) == 1L) {

      if (x.iterations % x.print_interval == 0) {
        Rprintf("iter = %d  f = %.8f  ng = %.8f\r",
                x.iterations, x.f, x.ng);
        R_FlushConsole();
        R_ProcessEvents();
      }

    }

    // if ((x.ng/std::sqrt(x.rg.n_elem)) < x.eps && x.max_rg < x.eps) {
    if ((x.ng/std::sqrt(x.rg.n_elem)) < x.eps) {
      x.convergence = true;
      break;
    }

  } while (x.iterations < x.maxit);

  optim_result result = std::make_tuple(x.parameters,
                                        x.transparameters,
                                        x.f,
                                        x.iterations,
                                        x.convergence,
                                        x.ng, x.rg, x.g, x.dir);

  return result;

}
