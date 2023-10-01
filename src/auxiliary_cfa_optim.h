/*
 * Author: Marcos Jimenez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 17/09/2023
 *
 */

// Line-search satisfying the armijo condition:

void armijo(arguments_optim& x, std::vector<arguments_cfa>& structs,
            double ss_fac, double ss_min, double max_iter,
            double c1, double c2, double eps) {

  // int nblocks = structs.size();
  cfa_manifold2* manifold = new ultimate_manifold();
  cfa_criterion2* criterion = new ultimate_criterion();

  x.ss = std::max(ss_min, x.ss * ss_fac);
  // x.ss = x.ss*2;
  double f0 = x.f;
  int iteration = 0;
  arma::vec parameters = x.parameters;
  x.inprod = arma::accu(x.dir % x.rg);

  do{

    ++iteration;
    x.parameters = parameters + x.ss*x.dir;
    // for(int i=0; i < nblocks; ++i) structs[i].parameters = x.parameters;
    // Projection onto the manifold
    manifold->retr(x, structs);
    // Rcpp::Rcout<<structs.phi[0]<<std::endl;
    // Rcpp::Rcout<<structs.T[0]<<std::endl;
    // Parameterization
    manifold->param(x, structs);
    criterion->F(x, structs);
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

void tcg(arguments_optim& x, std::vector<arguments_cfa>& structs,
         arma::vec& dir, bool& att_bnd, double ng, arma::vec c, double rad) {

  /*
   * Truncated conjugate gradient sub-solver for the trust-region sub-problem
   * From Liu (Algorithm 4; 2020)
   */

  cfa_manifold2* manifold = new ultimate_manifold();
  cfa_criterion2* criterion = new ultimate_criterion();

  dir.zeros();
  arma::vec dir0;

  double alpha, rr0, tau, beta, dHd;
  x.dparameters = -x.rg; // Initial search direction
  arma::vec r = x.dparameters; // Initial residual
  double rr = ng * ng;
  double tol = ng * std::min(pow(ng, c[0]), c[1]);

  int iter = 0;

  do{

    // Differential of transformed parameters
    manifold->dparam(x, structs);

    // Differential of the euclidean gradient
    criterion->dG(x, structs);

    // Differential of g
    manifold->dgrad(x, structs);

    // Riemannian hessian
    manifold->hess(x, structs);

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
    iter = iter + 1;

  } while (iter < 5);

}

// Newton Trust-region algorithm:

void ntr(arguments_optim& x, std::vector<arguments_cfa>& structs) {

  /*
   * Riemannian trust-region algorithm
   * From Liu (Algorithm 2; 2020)
   */

  cfa_manifold2* manifold = new ultimate_manifold();
  cfa_criterion2* criterion = new ultimate_criterion();

  // Parameterization
  manifold->param(x, structs);

  // Objective
  criterion->F(x, structs);

  // Gradient wrt transformed parameters
  criterion->G(x, structs);

  // Gradient
  manifold->grad(x, structs);

  // Riemannian gradient
  manifold->proj(x, structs);

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

  arguments_optim new_x;
  arma::vec dir(x.parameters.size());

  x.convergence = false;

  do{

    // subsolver
    tcg(x, structs, dir, att_bnd, x.ng, c, rad);
    x.dparameters = dir;
    new_x = x;
    new_x.parameters += dir;

    // Projection onto the manifold
    manifold->retr(new_x, structs);

    // Differential of transformed parameters
    manifold->dparam(x, structs);

    // Differential of the gradient wrt transformed parameters
    criterion->dG(x, structs);

    // Differential of g
    manifold->dgrad(x, structs);

    // Riemannian hessian
    manifold->hess(x, structs);

    preddiff = - arma::accu(x.dparameters % ( x.rg + 0.5 * x.dH) );

    // Parameterization
    manifold->param(new_x, structs);

    // objective
    criterion->F(new_x, structs);

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
      criterion->G(x, structs);
      manifold->grad(x, structs);

      // Riemannian gradient
      manifold->proj(x, structs);

      x.ng = sqrt(arma::accu(x.rg % x.rg));

    }

    ++x.iteration;
    if (x.ng < x.eps) {
      x.convergence = true;
      break;
    }

  } while (x.iteration < x.maxit);

}

// Gradient descent algorithm:

void gd(arguments_optim& x, std::vector<arguments_cfa>& structs) {

  cfa_manifold2* manifold = new ultimate_manifold();
  cfa_criterion2* criterion = new ultimate_criterion();

  x.iteration = 0;
  double ss_fac = 2, ss_min = 0.1, c1 = 0.5, c2 = 0.5;

  // Parameterization
  manifold->param(x, structs);
  criterion->F(x, structs);

  // update gradient
  criterion->G(x, structs);
  manifold->grad(x, structs);

  // Riemannian gradient
  manifold->proj(x, structs);
  x.dir = -x.rg;
  x.inprod = arma::accu(-x.dir % x.rg);
  x.ng = sqrt(x.inprod);
  // x.ss = 1;

  x.convergence = false;

  do{

    // x.ss *= 2;

    armijo(x, structs, ss_fac, ss_min,
           10, c1, c2, x.eps);

    // Update the gradient
    criterion->G(x, structs);
    manifold->grad(x, structs);

    // Riemannian gradient
    manifold->proj(x, structs);
    x.dir = -x.rg;
    x.inprod = arma::accu(-x.dir % x.rg);
    x.ng = sqrt(x.inprod);

    ++x.iteration;
    if (x.ng < x.eps) {
      x.convergence = true;
      break;
    }

  } while (x.iteration < x.maxit);

}

// BFGS algorithm:

void bfgs(arguments_optim& x, std::vector<arguments_cfa>& structs) {

  cfa_manifold2* manifold = new ultimate_manifold();
  cfa_criterion2* criterion = new ultimate_criterion();

  x.iteration = 0;
  double ss_fac = 2, ss_min = 0.1, c1 = 10e-04, c2 = 0.5;

  // Parameterization
  manifold->param(x, structs); // update x.L, x.Phi and x.Inv_T
  criterion->F(x, structs);
  // update the gradient
  criterion->G(x, structs);
  manifold->grad(x, structs);
  // Riemannian gradient
  manifold->proj(x, structs);
  x.dir = -x.rg;
  x.inprod = arma::accu(-x.dir % x.rg);
  x.ng = sqrt(x.inprod);
  // x.ss = 1;
  int p = x.parameters.size();
  arma::mat B(p, p, arma::fill::eye);

  x.convergence = false;

  do{

    // x.ss *= 2;

    arma::mat old_parameters = x.parameters;
    arma::mat old_rg = x.rg;

    armijo(x, structs, ss_fac, ss_min, 30, c1, c2, x.eps);

    // update gradient
    criterion->G(x, structs);
    manifold->grad(x, structs);
    // Riemannian gradient
    manifold->proj(x, structs);

    arma::vec y = x.rg - old_rg;
    arma::vec s = x.parameters - old_parameters;
    double sy = arma::accu(s%y);
    B += (sy + y.t() * B * y) % (s * s.t()) / (sy*sy) -
      (B * y * s.t() + s * y.t() * B) / sy;
    arma::vec dir = B * x.rg;

    x.dir = -dir;
    x.inprod = arma::accu(-x.dir % x.rg);
    x.ng = sqrt(x.inprod);

    ++x.iteration;
    if (x.ng < x.eps) {
      x.convergence = true;
      break;
    }

  } while (x.iteration < x.maxit);

  x.B = B;

}

// L-BFGS algorithm:

void lbfgs(arguments_optim& x, std::vector<arguments_cfa>& structs) {

  cfa_manifold2* manifold = new ultimate_manifold();
  cfa_criterion2* criterion = new ultimate_criterion();

  x.iteration = 0;
  double ss_fac = 2, ss_min = 0.1, c1 = 10e-04, c2 = 0.5;

  // Parameterization
  // Rcpp::Rcout<<"param"<<std::endl;
  manifold->param(x, structs); // update x.L, x.Phi and x.Inv_T
  // Rcpp::Rcout<<"F"<<std::endl;
  criterion->F(x, structs);    // Compute the objective with x.L and x.Phi
  // update the gradient
  // Rcpp::Rcout<<"G"<<std::endl;
  criterion->G(x, structs);  // Update the gradient wrt x.L and x.Phi
  // Rcpp::Rcout<<"grad"<<std::endl;
  manifold->grad(x, structs);  // Update the gradient wrt x.T
  // Rcpp::Rcout<<"proj"<<std::endl;
  // Riemannian gradient
  manifold->proj(x, structs);  // Update the Riemannian gradient x.rg
  // Rcpp::Rcout<<"after proj"<<std::endl;
  x.dir = -x.rg;
  x.inprod = arma::accu(-x.dir % x.rg);
  x.ng = std::sqrt(x.inprod);
  // x.ss = 1;
  int pp = x.parameters.size();
  arma::mat B(pp, pp, arma::fill::eye);

  int M = 15;
  std::vector<arma::vec> s(x.maxit), y(x.maxit);
  std::vector<double> p(x.maxit), alpha(x.maxit), beta(x.maxit);

  x.convergence = false;
  // Rcpp::Rcout<<"before do"<<std::endl;

  do{

    // x.ss *= 2;

    int k = x.iteration;
    arma::uvec seq(2);
    seq[0] = M; seq[1] = k;
    int min = seq.min();
    arma::vec max(2);
    max[0] = min; max[1] = 0;
    int m = max.max();

    arma::vec old_parameters = x.parameters;
    arma::vec old_rg = x.rg;

    // Update x.ss, x.T, x.L, x.Phi and x.Inv_T and x.f
    armijo(x, structs, ss_fac, ss_min, 30, c1, c2, x.eps);

    // update gradient
    criterion->G(x, structs);
    manifold->grad(x, structs);
    // Riemannian gradient
    manifold->proj(x, structs);

    arma::vec q = x.rg;
    s[k] = x.parameters - old_parameters;
    y[k] = x.rg - old_rg;
    p[k] = 1/arma::accu(y[k] % s[k]);

    for(int i=k; i > (k-m-1); --i) {

      alpha[i] = p[i]*arma::accu(s[i] % q);
      q -= alpha[i] * y[i];

    }

    double gamma = arma::accu(s[k] % y[k]) / arma::accu(y[k] % y[k]);
    arma::mat H0 = gamma*B;
    arma::vec z = H0 * q;

    for(int i=(k-m); i < (k+1); ++i) {

      beta[i] = p[i]*arma::accu(y[i] % z);
      z += s[i] * (alpha[i] - beta[i]);

    }

    x.dir = -z;
    x.inprod = arma::accu(-x.dir % x.rg);
    x.ng = sqrt(x.inprod);

    ++x.iteration;
    if (x.ng < x.eps) {
      x.convergence = true;
      break;
    }

  } while (x.iteration < x.maxit);
  // Rcpp::Rcout<<"end_optim"<<std::endl;

  // cfa_NTR result = std::make_tuple(x.parameters, x.f, x.iteration, x.convergence);

  // return result;

}

