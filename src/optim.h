/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 20/07/2025
 */

typedef std::tuple<arma::vec, arma::vec, double, int, bool, double, arma::mat, arma::mat> optim_result;

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

// Line-search algorithm satisfying the Wolfe conditions:

void wolfe(arguments_optim& x,
           product_transform* final_transform,
           product_manifold* final_manifold,
           product_estimator* final_estimator,
           std::vector<transformations*>& xtransforms,
           std::vector<manifolds*>& xmanifolds,
           std::vector<estimators*>& xestimators) {

  // x.ss = std::max(x.ss_min, x.ss * x.ss_fac);
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
    final_transform->transform(x, xtransforms);
    final_estimator->param(x, xestimators);
    final_estimator->F(x, xestimators);
    final_estimator->G(x, xestimators);
    final_transform->update_grad(x, xtransforms);
    final_manifold->param(x, xmanifolds);
    final_manifold->proj(x, xmanifolds);
    x.df = x.f - f0;
    double inprod = arma::accu(x.dir % x.rg); // Armijo condition
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

// Conjugate-gradient method to solve the Riemannian Newton equation:

void tcg(arguments_optim& x, product_manifold* final_manifold, product_estimator* final_estimator,
         std::vector<manifolds*>& xmanifolds, std::vector<estimators*>& xestimators,
         arma::vec& dir, bool& att_bnd, double ng, arma::vec c, double rad) {

  /*
   * Truncated conjugate gradient sub-solver for the trust-region sub-problem
   * From Liu (Algorithm 4; 2020)
   */

  dir.zeros();
  arma::vec dir0;

  double alpha, rr0, tau, beta, dHd;
  x.dparameters = -x.rg; // Initial search direction
  arma::vec r = x.dparameters; // Initial residual
  double rr = ng * ng;
  double tol = ng * std::min(pow(ng, c[0]), c[1]);

  int iter = 0;

  // In x, g should be already computed
  // final_estimator->param(x, xestimators); // Unnecessary
  // final_estimator->G(x, xestimators); // Unnecessary

  do {

    // final_estimator->param(x, xestimators); // Unnecessary
    final_estimator->dG(x, xestimators);
    final_manifold->param(x, xmanifolds);
    // final_manifold->proj(x, xmanifolds); // Unnecessary
    final_manifold->hess(x, xmanifolds);

    dHd = arma::accu(x.dparameters % x.dH);

    if(dHd <= 0) {

      tau = root_quad(arma::accu(x.dparameters % x.dparameters), 2 * arma::accu(dir % x.dparameters),
                      arma::accu(dir % dir) - rad * rad); // Solve equation 39
      dir = dir + tau * x.dparameters;
      att_bnd = true;

      break;

    }

    rr0 = rr;
    alpha = rr0 / dHd;
    dir0 = dir;
    dir = dir + alpha * x.dparameters; // update proposal

    if (sqrt(arma::accu(dir % dir)) >= rad) {

      tau = root_quad(arma::accu(x.dparameters % x.dparameters), 2 * arma::accu(dir0 % x.dparameters),
                      arma::accu(dir0 % dir0) - rad * rad); // Solve equation 39
      dir = dir0 + tau * x.dparameters;
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
    x.dir = dir;
    iter += 1;

  } while (iter < x.tcg_maxit);

}

// Gradient descent algorithm:

optim_result gd(arguments_optim x,
                std::vector<transformations*>& xtransforms,
                std::vector<manifolds*>& xmanifolds,
                std::vector<estimators*>& xestimators) {

  product_manifold* final_manifold;
  product_transform* final_transform;
  product_estimator* final_estimator;

  // Ensure initial parameters are ok:
  final_manifold->param(x, xmanifolds);
  final_manifold->retr(x, xmanifolds);
  final_manifold->param(x, xmanifolds);
  final_estimator->param(x, xestimators);

  // double ss_fac = 2, ss_min = 0.1;
  // x.c1 = 0.5; x.c2 = 0.5;
  // x.step_maxit = 10L;

  // Parameterization

  // final_estimator->param(x, xestimators);
  final_estimator->F(x, xestimators);
  // update gradient
  final_estimator->G(x, xestimators);
  // Riemannian gradient
  // final_manifold->param(x, xmanifolds);
  final_manifold->proj(x, xmanifolds);

  x.dir = -x.rg;
  x.inprod = arma::accu(-x.dir % x.rg);
  x.ng = sqrt(x.inprod);

  do {

    // x.ss *= 2;

    // armijo(x, final_manifold, final_estimator, xmanifolds, xestimators);
    wolfe(x, final_transform, final_manifold, final_estimator,
          xtransforms, xmanifolds, xestimators);

    // update gradient
    // final_estimator->param(x, xestimators); // param() was called in armijo
    final_estimator->G(x, xestimators);
    // Riemannian gradient
    final_manifold->param(x, xmanifolds);
    final_manifold->proj(x, xmanifolds);

    x.dir = -x.rg;
    x.inprod = arma::accu(-x.dir % x.rg);
    x.ng = sqrt(x.inprod);

    // if(x.old_inprod < x.inprod) {
    //   Rcpp::Rcout << "Continue 2" << std::endl;
    //   continue;
    // } else {
    //   x.old_inprod = x.inprod;
    // }

    ++x.iterations;

    if(x.print) {
      Rcpp::Rcout << "Iteration = " << x.iterations << std::endl;
      Rcpp::Rcout << "f = " << x.f << std::endl;
      Rcpp::Rcout << "step size = " << x.ss << std::endl;
      Rcpp::Rcout << "step iters = " << x.step_iteration << std::endl;
      Rcpp::Rcout << "ng = " << x.ng << std::endl;
      Rcpp::Rcout << "dif = " << std::sqrt(x.df*x.df) << std::endl;
      Rcpp::Rcout << "" << std::endl;
    }
    if (x.ng < x.eps | std::sqrt(x.df*x.df) < x.df_eps) {
      x.convergence = true;
      break;
    }

  } while (x.iterations < x.maxit);

  // final_estimator->outcomes(x, xestimators);

  optim_result result = std::make_tuple(x.parameters,
                                        x.transparameters,
                                        x.f,
                                        x.iterations,
                                        x.convergence,
                                        x.ng, x.rg, x.posterior);

  return result;

}

// L-BFGS algorithm:

optim_result lbfgs(arguments_optim x,
                   std::vector<transformations*>& xtransforms,
                   std::vector<manifolds*>& xmanifolds,
                   std::vector<estimators*>& xestimators) {

  product_manifold* final_manifold;
  product_transform* final_transform;
  product_estimator* final_estimator;

  // Ensure initial parameters are ok:
  // Rf_error("292");
  final_manifold->param(x, xmanifolds);
  final_manifold->retr(x, xmanifolds);
  final_manifold->param(x, xmanifolds);
  // Rf_error("296");
  final_transform->transform(x, xtransforms);
  // Rf_error("298");
  // Rf_error("299");
  final_estimator->param(x, xestimators);
  // Rf_error("300");
  x.iterations = 0L;

  // Parameterization
  // final_estimator->param(x, xestimators);
  final_estimator->F(x, xestimators);
  // Rprintf("307= %g \n", x.f);
  // update the gradient
  final_estimator->G(x, xestimators);
  // Rf_error("307");
  final_transform->update_grad(x, xtransforms);
  // Rf_error("306");
  // Riemannian gradient
  // final_manifold->param(x, xmanifolds);
  final_manifold->proj(x, xmanifolds);
  x.dir = -x.rg;
  x.inprod = arma::accu(-x.dir % x.rg);
  x.ng = sqrt(x.inprod);
  // x.ss = 1;
  int p1 = x.parameters.size();
  arma::mat B(p1, p1, arma::fill::eye);

  std::vector<arma::vec> s(x.maxit), y(x.maxit);
  std::vector<double> p(x.maxit), alpha(x.maxit), beta(x.maxit);

  x.convergence = false;

  do {

    // x.ss *= 2;

    int k = x.iterations;
    ++x.iterations;

    arma::uvec seq(2);
    seq[0] = x.M; seq[1] = k;
    int min = seq.min();
    arma::vec max(2);
    max[0] = min; max[1] = 0;
    int m = max.max();

    arma::vec old_parameters = x.parameters;
    arma::vec old_rg = x.rg;

    // for (arma::uword i = 0; i < x.parameters.n_elem; ++i) {
    //   Rprintf("%g ", x.parameters[i]);
    // }
    // Rprintf("\n");
    //
    // for (arma::uword i = 0; i < x.g.n_elem; ++i) {
    //   Rprintf("%g ", x.g[i]);
    // }
    // Rprintf("\n");
    //
    // for (arma::uword i = 0; i < x.rg.n_elem; ++i) {
    //   Rprintf("%g ", x.rg[i]);
    // }
    // Rprintf("\n");

    // armijo(x, final_manifold, final_estimator, xmanifolds, xestimators);
    wolfe(x, final_transform, final_manifold, final_estimator,
          xtransforms, xmanifolds, xestimators);
    // Rprintf("363= %g \n", x.f);

    // for (arma::uword i = 0; i < x.parameters.n_elem; ++i) {
    //   Rprintf("%g ", x.parameters[i]);
    // }
    // Rprintf("\n");

    // if(x.ss < arma::datum::eps) {
    //   x.convergence = false;
    //   break;
    // }

    // final_estimator->param(x, xestimators); // Necessary¿?
    final_estimator->G(x, xestimators);
    final_transform->update_grad(x, xtransforms);
    // After the retraction in armijo you need to param the manifolds:
    final_manifold->param(x, xmanifolds);
    final_manifold->proj(x, xmanifolds);

    arma::vec q = arma::vectorise(x.rg);
    s[k] = arma::vectorise(x.parameters - old_parameters);
    y[k] = arma::vectorise(x.rg - old_rg);
    p[k] = 1/arma::accu(y[k] % s[k]);

    for(int i=k; i > (k-m-1); --i) {

      alpha[i] = p[i]*arma::accu(s[i] % q);
      q -= alpha[i] * y[i];

    }

    double gamma = arma::accu(s[k] % y[k]) / arma::accu(y[k] % y[k]);
    if(gamma < arma::datum::eps) {
      // Rcpp::Rcout << "Continue" << std::endl;
      // x.dir = -x.rg;
      continue;
    }
    arma::mat H0 = gamma*B;
    arma::vec z = H0 * q;

    for(int i=(k-m); i < (k+1); ++i) {

      beta[i] = p[i]*arma::accu(y[i] % z);
      z += s[i] * (alpha[i] - beta[i]);

    }

    x.dir = -z;
    x.inprod = arma::accu(-x.dir % x.rg);
    // x.inprod = arma::accu(x.dir % x.dir);
    x.ng = sqrt(x.inprod);
    if (std::isnan(x.ng)) {
      x.convergence = false;
      break;
    }
    x.old_inprod = x.inprod;
    // if(x.old_inprod < x.inprod) {
    //   // Rcpp::Rcout << "Continue 2" << std::endl;
    //   continue;
    // } else {
    //   x.old_inprod = x.inprod;
    // }

    if(x.print) {
      Rcpp::Rcout << "Iteration = " << x.iterations << std::endl;
      Rcpp::Rcout << "f = " << x.f << std::endl;
      Rcpp::Rcout << "step size = " << x.ss << std::endl;
      Rcpp::Rcout << "step iters = " << x.step_iteration << std::endl;
      Rcpp::Rcout << "ng = " << x.ng << std::endl;
      Rcpp::Rcout << "dif = " << std::sqrt(x.df*x.df) << std::endl;
      Rcpp::Rcout << "" << std::endl;
    }

    // if (x.ng < x.eps | std::sqrt(x.df*x.df) < x.df_eps) {
    if (x.ng < x.eps) {
      x.convergence = true;
      break;
    }

  } while (x.iterations < x.maxit);

  // final_estimator->outcomes(x, xestimators);

  optim_result result = std::make_tuple(x.parameters,
                                        x.transparameters,
                                        x.f,
                                        x.iterations,
                                        x.convergence,
                                        x.ng, x.rg, x.posterior);

  return result;

}

// Newton Trust-region algorithm:

optim_result ntr(arguments_optim x,
                 std::vector<transformations*>& xtransforms,
                 std::vector<manifolds*>& xmanifolds,
                 std::vector<estimators*>& xestimators) {

  product_manifold* final_manifold;
  product_estimator* final_estimator;

  // Ensure initial parameters are ok:
  final_manifold->param(x, xmanifolds);
  final_manifold->retr(x, xmanifolds);
  final_manifold->param(x, xmanifolds);
  final_estimator->param(x, xestimators);

  /*
   * Riemannian trust-region algorithm
   * From Liu (Algorithm 2; 2020)
   */

  // Parameterization
  // final_manifold->param(x, xmanifolds);
  // final_manifold->retr(x, xmanifolds);
  // final_manifold->param(x, xmanifolds);
  // final_estimator->param(x, xestimators);

  // Objective
  final_estimator->F(x, xestimators);

  // Rcpp::Rcout << "x.f = " << x.f << std::endl;

  // Gradient
  final_estimator->G(x, xestimators);

  // Riemannian gradient
  final_manifold->proj(x, xmanifolds);

  x.ng = sqrt(arma::accu(x.rg % x.rg));

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

  arguments_optim new_x;
  arma::vec dir(x.parameters.n_elem);

  do {

    // subsolver
    tcg(x, final_manifold, final_estimator, xmanifolds, xestimators,
        dir, att_bnd, x.ng, c, rad); // Update x.dparameters, x.dg, and x.dH
    x.dparameters = dir;
    x.dir = dir;
    new_x = x;
    new_x.parameters += dir;

    final_estimator->dG(x, xestimators);
    final_manifold->param(x, xmanifolds);
    final_manifold->hess(x, xmanifolds);

    preddiff = - arma::accu(x.dparameters % ( x.rg + 0.5 * x.dH) );

    // Projection onto the manifold
    final_manifold->param(new_x, xmanifolds);
    final_manifold->retr(new_x, xmanifolds);
    // Parameterization
    final_estimator->param(new_x, xestimators);
    final_estimator->F(new_x, xestimators);

    x.df = x.f - new_x.f; // New stopping criteria

    if ( std::abs(preddiff) <= arma::datum::eps ) {

      goa = 1;

    } else {

      goa = (x.f - new_x.f) / preddiff;

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
      final_estimator->param(x, xestimators); // Unnecessary
      final_estimator->G(x, xestimators);
      // Riemannian gradient
      final_manifold->param(x, xmanifolds);
      final_manifold->proj(x, xmanifolds);

      x.ng = sqrt(arma::accu(x.rg % x.rg));

    }

    ++x.iterations;

    if(x.print) {
      Rcpp::Rcout << "Iteration = " << x.iterations << std::endl;
      Rcpp::Rcout << "f = " << x.f << std::endl;
      Rcpp::Rcout << "step size = " << x.ss << std::endl;
      Rcpp::Rcout << "step iters = " << x.step_iteration << std::endl;
      Rcpp::Rcout << "ng = " << x.ng << std::endl;
      Rcpp::Rcout << "dif = " << std::sqrt(x.df*x.df) << std::endl;
      Rcpp::Rcout << "" << std::endl;
    }

    if (x.ng < x.eps | std::sqrt(x.df*x.df) < x.df_eps) {
      x.convergence = true;
      break;
    }

  } while (x.iterations < x.maxit);

  // final_estimator->outcomes(x, xestimators);

  optim_result result = std::make_tuple(x.parameters,
                                        x.transparameters,
                                        x.f,
                                        x.iterations,
                                        x.convergence,
                                        x.ng, x.rg, x.posterior);

  return result;

}

// Newton Trust-region algorithm:

optim_result em(arguments_optim x,
                std::vector<transformations*>& xtransforms,
                std::vector<manifolds*>& xmanifolds,
                std::vector<estimators*>& xestimators) {

  product_estimator* final_estimator;

  x.convergence = false;
  arma::vec loglik(x.maxit); loglik[0] = x.loglik;
  x.iterations = 0L;

  do {

    ++x.iterations;
    final_estimator->E(x, xestimators);
    final_estimator->M(x, xestimators);
    loglik[x.iterations] = x.loglik;
    double diff = loglik[x.iterations]-loglik[x.iterations-1];
    x.ng = sqrt(arma::accu(diff*diff));
    if (x.ng < x.eps) {
      x.convergence = true;
      break;
    }
    // Rf_error("617");

  } while(x.iterations < x.maxit);

  // final_estimator->outcomes(x, xestimators);

  optim_result result = std::make_tuple(x.parameters,
                                        x.transparameters,
                                        x.loglik,
                                        x.iterations,
                                        x.convergence,
                                        x.ng, x.rg, x.posterior);

  return result;

}

// Optimization algorithms:

class optim {

public:

  virtual optim_result optimize(arguments_optim x,
                                std::vector<transformations*>& xtransforms,
                                std::vector<manifolds*>& xmanifolds,
                                std::vector<estimators*>& xestimators) = 0;

};

// Riemannian Gradient Descent:

class RGD:public optim {

public:

  optim_result optimize(arguments_optim x,
                        std::vector<transformations*>& xtransforms,
                        std::vector<manifolds*>& xmanifolds,
                        std::vector<estimators*>& xestimators) {

    return gd(x, xtransforms, xmanifolds, xestimators);

  }

};

// L-BFGS:

class LBFGS:public optim {

public:

  optim_result optimize(arguments_optim x,
                        std::vector<transformations*>& xtransforms,
                        std::vector<manifolds*>& xmanifolds,
                        std::vector<estimators*>& xestimators) {

    return lbfgs(x, xtransforms, xmanifolds, xestimators);

  }

};

// Riemannian Newton Trust-Region:

class RNTR:public optim {

public:

  optim_result optimize(arguments_optim x,
                        std::vector<transformations*>& xtransforms,
                        std::vector<manifolds*>& xmanifolds,
                        std::vector<estimators*>& xestimators) {

    return ntr(x, xtransforms, xmanifolds, xestimators);

  }

};

class EM:public optim {

public:

  optim_result optimize(arguments_optim x,
                        std::vector<transformations*>& xtransforms,
                        std::vector<manifolds*>& xmanifolds,
                        std::vector<estimators*>& xestimators) {

    return em(x, xtransforms, xmanifolds, xestimators);

  }

};

optim* choose_optim(arguments_optim& x, Rcpp::List control_optimizer) {

  // if(control_optimizer.containsElementNamed("parameters")) {
    std::vector<arma::vec> params = control_optimizer["parameters"];
    x.nparam = params[0].n_elem;
    x.dir.set_size(x.nparam); x.dir.zeros();
    x.parameters.set_size(x.nparam);
  // }

  // if(control_optimizer.containsElementNamed("transparameters")) {
    std::vector<arma::vec> transparameters = control_optimizer["transparameters"];
    x.ntransparam = transparameters[0].n_elem;
    x.transparameters.set_size(x.ntransparam);
  // }

  arma::uvec param2transparam = control_optimizer["param2transparam"];
  arma::uvec transparam2param = control_optimizer["transparam2param"];
  x.param2transparam = param2transparam;
  x.transparam2param = transparam2param;

  if(control_optimizer.containsElementNamed("posterior")) {
    std::vector<arma::mat> post = control_optimizer["posterior"];
    x.nrow_post = post[0].n_rows;
    x.ncol_post = post[0].n_cols;
    x.posterior.set_size(x.nrow_post, x.ncol_post);
    x.latentloglik.set_size(x.nrow_post, x.ncol_post);
  }

  if(control_optimizer.containsElementNamed("nlatent")) {
    int S = control_optimizer["S"];
    int nlatent = control_optimizer["nlatent"];
    x.latentloglik.set_size(S, nlatent);
    x.latentpars.set_size(nlatent);
    x.loglatentpars.set_size(nlatent);
  }

  if(control_optimizer.containsElementNamed("opt")) {
    std::string opt = control_optimizer["opt"];
    x.optimizer = opt;
  }
  if(control_optimizer.containsElementNamed("maxit")) {
    int maxit = control_optimizer["maxit"];
    x.maxit = maxit;
  }
  if(control_optimizer.containsElementNamed("step_maxit")) {
    int step_maxit = control_optimizer["step_maxit"];
    x.step_maxit = step_maxit;
  }
  if(control_optimizer.containsElementNamed("step_eps")) {
    double step_eps = control_optimizer["step_eps"];
    x.step_eps = step_eps;
  }
  if(control_optimizer.containsElementNamed("df_eps")) {
    double df_eps = control_optimizer["df_eps"];
    x.df_eps = df_eps;
  }
  if(control_optimizer.containsElementNamed("M")) {
    int M = control_optimizer["M"];
    x.M = M;
  }
  if(control_optimizer.containsElementNamed("rstarts")) {
    int rstarts = control_optimizer["rstarts"];
    x.rstarts = rstarts;
  }
  if(control_optimizer.containsElementNamed("cores")) {
    int cores = control_optimizer["cores"];
    x.cores = cores;
  }
  if(control_optimizer.containsElementNamed("c2")) {
    double c2 = control_optimizer["c2"];
    x.c2 = c2;
  }
  if(control_optimizer.containsElementNamed("ss")) {
    double ss = control_optimizer["ss"];
    x.ss = ss;
  }
  if(control_optimizer.containsElementNamed("ss_fac")) {
    double ss_fac = control_optimizer["ss_fac"];
    x.ss_fac = ss_fac;
  }
  if(control_optimizer.containsElementNamed("ss_min")) {
    double ss_min = control_optimizer["ss_min"];
    x.ss_min = ss_min;
  }
  if(control_optimizer.containsElementNamed("eps")) {
    double eps = control_optimizer["eps"];
    x.eps = eps;
  }
  if(control_optimizer.containsElementNamed("tcg_maxit")) {
    double tcg_maxit = control_optimizer["tcg_maxit"];
    x.tcg_maxit = tcg_maxit;
  }
  if(control_optimizer.containsElementNamed("print")) {
    bool print = control_optimizer["print"];
    x.print = print;
  }

  // Select the optimization algorithm:

  optim* algorithm;
  if(x.optimizer == "grad") {

    if(control_optimizer.containsElementNamed("c1")) {
      double c1 = control_optimizer["c1"];
      x.c1 = c1;
    } else {
      x.c1 = 0.5;
    }

    algorithm = new RGD();

  } else if(x.optimizer == "lbfgs") {

    if(control_optimizer.containsElementNamed("c1")) {
      double c1 = control_optimizer["c1"];
      x.c1 = c1;
    } else {
      // x.c1 = 10e-04;
      x.c1 = 0.5;
    }

    algorithm = new LBFGS();

  } else if(x.optimizer == "newton") {
    algorithm = new RNTR();
  } else if(x.optimizer == "em") {

    algorithm = new EM();

  } else {

    Rcpp::stop("Available optimization rutines for factor extraction: \n grad, lbfgs, newton, em");

  }

  return algorithm;

}

