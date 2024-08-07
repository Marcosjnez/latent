/*
 * Author: Marcos Jimenez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 18/03/2022
 *
 */

// #include "structures.h"
// #include "auxiliary_manifolds.h"
// #include "manifold.h"
// #include "auxiliary_criteria.h"
// #include "criteria.h"

// Line-search satisfying the armijo condition:

void armijo(arguments_rotate& x, rotation_manifold *manifold,
            rotation_criterion *criterion,
            double ss_fac, double ss_min, double max_iter,
            double c1, double c2, double eps) {

  x.ss = std::max(ss_min, x.ss * ss_fac);
  // x.ss = x.ss*2;
  double f0 = x.f;
  int iteration = 0;
  arma::mat X = x.T;
  x.inprod = arma::accu(x.dir % x.rg);

  do{

    ++iteration;
    x.T = X + x.ss*x.dir;
    // Projection onto the manifold
    manifold->retr(x); // update x.T
    // Parameterization
    manifold->param(x); // update x.L, x.Phi and x.Inv_T
    criterion->F(x);
    double df = x.f - f0;
    if (df < c1 * x.ss * x.inprod || x.ss < 1e-09) // armijo condition
      break;
    x.ss *= c2;

  } while (iteration <= max_iter);

  bool convergence = true;
  if(iteration > max_iter) {

    convergence = false;

  }

}

// Conjugate-gradient method to solve the Riemannian Newton equation in the Trust-Region algorithm:

void tcg(arguments_rotate x, rotation_manifold *manifold, rotation_criterion *criterion,
         arma::mat& dir, bool& att_bnd, double ng, arma::vec c, double rad) {

  /*
   * Truncated conjugate gradient sub-solver for the trust-region sub-problem
   * From Liu (Algorithm 4; 2020)
   */

  dir.zeros();
  arma::mat dir0;

  double alpha, rr0, tau, beta, dHd;
  x.dT = -x.rg; // Initial search direction
  arma::mat r = x.dT; // Initial residual
  double rr = ng * ng;
  double tol = ng * std::min(pow(ng, c[0]), c[1]);

  int iter = 0;

  do{

    // Differential of L and P
    manifold->dLP(x);

    // Differential of the gradient of L and P
    criterion->dgLP(x);

    // Differential of g
    manifold->dgrad(x);

    // Riemannian hessian
    manifold->hess(x);

    dHd = arma::accu(x.dT % x.dH);

    if(dHd <= 0) {

      tau = root_quad(arma::accu(x.dT % x.dT), 2 * arma::accu(dir % x.dT),
                      arma::accu(dir % dir) - rad * rad); // Solve equation 39
      dir = dir + tau * x.dT;
      att_bnd = true;

      break;

    }

    rr0 = rr;
    alpha = rr0 / dHd;
    dir0 = dir;
    dir = dir + alpha * x.dT; // update proposal

    if (sqrt(arma::accu(dir % dir)) >= rad) {

      tau = root_quad(arma::accu(x.dT % x.dT), 2 * arma::accu(dir0 % x.dT),
                      arma::accu(dir0 % dir0) - rad * rad); // Solve equation 39
      dir = dir0 + tau * x.dT;
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
    x.dT = r + beta * x.dT;
    iter = iter + 1;

  } while (iter < 5);

}

// Newton Trust-region algorithm:

NTR ntr(arguments_rotate x, rotation_manifold *manifold, rotation_criterion *criterion) {

  /*
   * Riemannian trust-region algorithm
   * From Liu (Algorithm 2; 2020)
   */

  // Parameterization
  manifold->param(x); // update x.L, x.Phi and x.Inv_T

  // Objective
  criterion->F(x); // update x.f, x.L2, x.IgCL2N, x.term, x.f1 and x.f2

  // Rcpp::Rcout << "x.f = " << x.f << std::endl;

  // Gradient wrt L and P
  criterion->gLP(x); // update x.gL, x.gP, x.f1, x.f2 and x.LoL2

  // Rcpp::Rcout << "x.gL = " << x.gL << std::endl;

  // Gradient wtr T
  manifold->grad(x); // update x.g

  // Riemannian gradient
  manifold->proj(x); // update x.rg and x.A

  // Differential of the gradient of L and P
  // criterion->dgLP(x); // update dgL and dgP

  // Rcpp::Rcout << "x.dgL = " << x.dgL << std::endl;

  x.ng = sqrt(arma::accu(x.rg % x.rg));

  double max_rad = 10;

  arma::vec fac_rad(2);
  fac_rad[0] = 0.25;
  fac_rad[1] = 2;

  arma::vec crit_goa(3);
  crit_goa[0] = 0.2;
  crit_goa[1] = 0.25;
  crit_goa[2] = 0.75;

  arma::vec c(2);
  c[0] = 1;
  c[1] = 0.01;

  double rad = 1;
  bool att_bnd = false;

  x.iteration = 0;
  double goa, preddiff;

  arguments_rotate new_x;
  arma::mat dir(x.q, x.q);

  x.convergence = false;

  do{

    // subsolver
    tcg(x, manifold, criterion, dir, att_bnd, x.ng, c, rad);
    x.dT = dir;
    x.dir = dir;
    new_x = x;
    new_x.T += dir;

    // Projection onto the manifold
    manifold->retr(new_x); // update x.T

    // Differential of L and P
    manifold->dLP(x); // update x.dL, x.dP and Inv_T_dt

    // Differential of the gradient of L and P
    criterion->dgLP(x); // update dgL and dgP

    // Differential of g
    manifold->dgrad(x); // update dg

    // Riemannian hessian
    manifold->hess(x); // update dH

    preddiff = - arma::accu(x.dT % ( x.rg + 0.5 * x.dH) );

    // Parameterization
    manifold->param(new_x);

    // objective
    criterion->F(new_x);

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

      // update the gradient
      criterion->gLP(x);
      manifold->grad(x);

      // Riemannian gradient
      manifold->proj(x);

      x.ng = sqrt(arma::accu(x.rg % x.rg));

    }

    ++x.iteration;
    if (x.ng < x.eps) {
      x.convergence = true;
      break;
    }

  } while (x.iteration <= x.maxit);

  NTR result = std::make_tuple(x.L, x.Phi, x.T, x.f, x.iteration, x.convergence, x.dir);

  return result;

}

// Gradient descent algorithm:

NTR gd(arguments_rotate x, rotation_manifold *manifold, rotation_criterion *criterion) {

  x.iteration = 0;
  double ss_fac = 2, ss_min = 0.1, c1 = 0.5, c2 = 0.5;

  // Parameterization
  manifold->param(x); // update x.L, x.Phi and x.Inv_T
  criterion->F(x);
  // update gradient
  criterion->gLP(x);
  manifold->grad(x);
  // Riemannian gradient
  manifold->proj(x);
  x.dir = -x.rg;
  x.inprod = arma::accu(-x.dir % x.rg);
  x.ng = sqrt(x.inprod);
  // x.ss = 1;

  x.convergence = false;

  do{

    // x.ss *= 2;

    armijo(x, manifold, criterion, ss_fac, ss_min,
           10, c1, c2, x.eps);

    // update gradient
    criterion->gLP(x);
    manifold->grad(x);
    // Riemannian gradient
    manifold->proj(x);
    x.dir = -x.rg;
    x.inprod = arma::accu(-x.dir % x.rg);
    x.ng = sqrt(x.inprod);

    ++x.iteration;
    if (x.ng < x.eps) {
      x.convergence = true;
      break;
    }

  } while (x.iteration <= x.maxit);

  NTR result = std::make_tuple(x.L, x.Phi, x.T, x.f, x.iteration, x.convergence, x.dir);

  return result;

}

// BFGS algorithm:

NTR bfgs(arguments_rotate x, rotation_manifold *manifold, rotation_criterion *criterion) {

  x.iteration = 0;
  double ss_fac = 2, ss_min = 0.1, c1 = 10e-04, c2 = 0.5;

  // Parameterization
  manifold->param(x); // update x.L, x.Phi and x.Inv_T
  criterion->F(x);
  // update the gradient
  criterion->gLP(x);
  manifold->grad(x);
  // Riemannian gradient
  manifold->proj(x);
  x.dir = -x.rg;
  x.inprod = arma::accu(-x.dir % x.rg);
  x.ng = sqrt(x.inprod);
  // x.ss = 1;
  int p1 = x.T.n_rows;
  int p2 = x.T.n_cols;
  arma::mat B(p1*p2, p1*p2, arma::fill::eye);

  x.convergence = false;

  do{

    // x.ss *= 2;

    arma::mat old_T = x.T;
    arma::mat old_rg = x.rg;

    armijo(x, manifold, criterion, ss_fac, ss_min,
           30, c1, c2, x.eps);

    // update gradient
    criterion->gLP(x);
    manifold->grad(x);
    // Riemannian gradient
    manifold->proj(x);

    arma::vec y = arma::vectorise(x.rg - old_rg);
    arma::vec s = arma::vectorise(x.T - old_T);
    double sy = arma::accu(s%y);
    B += (sy + y.t() * B * y) % (s * s.t()) / (sy*sy) -
      (B * y * s.t() + s * y.t() * B) / sy;
    arma::mat dir = B * arma::vectorise(x.rg);
    dir.reshape(p1, p2);

    x.dir = -dir;
    x.inprod = arma::accu(-x.dir % x.rg);
    x.ng = sqrt(x.inprod);

    ++x.iteration;
    if (x.ng < x.eps) {
      x.convergence = true;
      break;
    }

  } while (x.iteration <= x.maxit);

  NTR result = std::make_tuple(x.L, x.Phi, x.T, x.f, x.iteration, x.convergence, x.dir);

  return result;

}

// L-BFGS algorithm:

NTR lbfgs(arguments_rotate x, rotation_manifold *manifold, rotation_criterion *criterion) {

  x.iteration = 0;
  double ss_fac = 2, ss_min = 0.1, c1 = 10e-04, c2 = 0.5;

  // Parameterization
  manifold->param(x); // update x.L, x.Phi and x.Inv_T
  criterion->F(x);    // Compute the objective with x.L and x.Phi
  // update the gradient
  criterion->gLP(x);  // Update the gradient wrt x.L and x.Phi
  manifold->grad(x);  // Update the gradient wrt x.T
  // Riemannian gradient
  manifold->proj(x);  // Update the Riemannian gradient x.rg
  x.dir = -x.rg;
  x.inprod = arma::accu(-x.dir % x.rg);
  x.ng = sqrt(x.inprod);
  // x.ss = 1;
  int p1 = x.T.n_rows;
  int p2 = x.T.n_cols;
  arma::mat B(p1*p2, p1*p2, arma::fill::eye);

  int M = 15;
  std::vector<arma::vec> s(x.maxit), y(x.maxit);
  std::vector<double> p(x.maxit), alpha(x.maxit), beta(x.maxit);

  x.convergence = false;

  do{

    // x.ss *= 2;

    int k = x.iteration;
    arma::uvec seq(2);
    seq[0] = M; seq[1] = k;
    int min = seq.min();
    arma::vec max(2);
    max[0] = min; max[1] = 0;
    int m = max.max();

    arma::mat old_T = x.T;
    arma::mat old_rg = x.rg;

    // Update x.ss, x.T, x.L, x.Phi and x.Inv_T and x.f
    armijo(x, manifold, criterion, ss_fac, ss_min,
           30, c1, c2, x.eps);

    // update gradient
    criterion->gLP(x);
    manifold->grad(x);
    // Riemannian gradient
    manifold->proj(x);

    arma::vec q = arma::vectorise(x.rg);
    s[k] = arma::vectorise(x.T - old_T);
    y[k] = arma::vectorise(x.rg - old_rg);
    p[k] = 1/arma::accu(y[k] % s[k]);

    for(int i=k; i > (k-m-1); --i) {

      alpha[i] = p[i]*arma::accu(s[i] % q);
      q -= alpha[i] * y[i];

    }

    double gamma = arma::accu(s[k] % y[k]) / arma::accu(y[k] % y[k]);
    arma::mat H0 = gamma*B;
    arma::mat z = H0 * q;

    for(int i=(k-m); i < (k+1); ++i) {

      beta[i] = p[i]*arma::accu(y[i] % z);
      z += s[i] * (alpha[i] - beta[i]);

    }

    z.reshape(p1, p2);

    x.dir = -z;
    x.inprod = arma::accu(-x.dir % x.rg);
    x.ng = sqrt(x.inprod);

    ++x.iteration;
    if (x.ng < x.eps) {
      x.convergence = true;
      break;
    }

  } while (x.iteration <= x.maxit);

  NTR result = std::make_tuple(x.L, x.Phi, x.T, x.f, x.iteration, x.convergence, x.dir);

  return result;

}
