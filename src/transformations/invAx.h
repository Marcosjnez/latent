/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 12/09/2026
 */

#ifndef LATENT_INVAX_H
#define LATENT_INVAX_H

#include <unordered_map>

// y = inv(A)*x, with A = I-B.
// Inputs are B (p x p) and x (length p), in that order.
// The output has length p. B need not be symmetric; I-B must be nonsingular.
// One pivoted LU factorization is reused for the value and derivatives.
// No full inverse is formed unless its entries are requested in jacobian().
// Derivative updates require transform() at the current parameter point.
// update_dgrad() also requires update_grad() and dtransform() there.

class invAx: public transformations {

public:

  int p;
  arma::uvec indices_B, indices_x, indices_invAx;
  arma::mat A, dB, grad_in_B;
  arma::vec v, y, dv, dy, grad_out, grad_in_x, output_weights;

  void transform(arguments_optim& x) {

    evaluate(x);
    x.transparameters(indices_invAx) = y;

  }

  void update_grad(arguments_optim& x) {

    grad_out = x.grad(indices_invAx)%output_weights;
    grad_in_x = solve_A(grad_out, true);
    grad_in_B = grad_in_x*y.t();

    // Explicit accumulation also handles repeated/shared input indices.
    for(arma::uword k = 0; k < indices_B.n_elem; ++k) {
      x.grad[indices_B[k]] += grad_in_B[k];
    }
    for(arma::uword k = 0; k < indices_x.n_elem; ++k) {
      x.grad[indices_x[k]] += grad_in_x[k];
    }

  }

  void dtransform(arguments_optim& x) {

    dB = arma::reshape(x.dtransparameters(indices_B), p, p);
    dv = x.dtransparameters(indices_x);

    // A*dy = dx+dB*y, since dA = -dB.
    dy = solve_A(dv+dB*y);
    x.dtransparameters(indices_invAx) = dy;

  }

  void update_dgrad(arguments_optim& x) {

    arma::vec dgrad_out = x.dgrad(indices_invAx)%output_weights;
    arma::vec dgrad_in_x = solve_A(dgrad_out+dB.t()*grad_in_x, true);
    arma::mat dgrad_in_B = dgrad_in_x*y.t()+grad_in_x*dy.t();

    for(arma::uword k = 0; k < indices_B.n_elem; ++k) {
      x.dgrad[indices_B[k]] += dgrad_in_B[k];
    }
    for(arma::uword k = 0; k < indices_x.n_elem; ++k) {
      x.dgrad[indices_x[k]] += dgrad_in_x[k];
    }

  }

  void jacobian(arguments_optim& x) {

    // Refresh caches without overwriting any parameter or derivative vector.
    evaluate(x);
    const arma::uword n = static_cast<arma::uword>(p);
    const arma::uword n2 = n*n;
    jacob.set_size(n, n2+n);

    // Column order: vec(B), x. J = [kron(y.t(), inv(A)), inv(A)].
    // Obtain each required inverse column by solving with a unit vector.
    arma::vec e(n, arma::fill::zeros);
    for(arma::uword i = 0; i < n; ++i) {
      e[i] = 1.00;
      arma::vec u = solve_A(e);
      e[i] = 0.00;

      jacob.col(n2+i) = u;
      for(arma::uword j = 0; j < n; ++j) {
        jacob.col(i+j*n) = y[j]*u;
      }
    }

  }

  void outcomes(arguments_optim& x) {

    (void)x;
    matrices.resize(1);
    matrices[0] = jacob;
    names_matrices.resize(1);
    names_matrices[0] = "jacobian";

  }

private:

  arma::mat L, U, Lt, Ut;
  arma::uvec permutation, inverse_permutation;

  void evaluate(arguments_optim& x) {

    A = -arma::reshape(x.transparameters(indices_B), p, p);
    A.diag() += 1.00;
    v = x.transparameters(indices_x); // Input vector x.

    if(!A.is_finite() || !v.is_finite()) {
      Rcpp::stop("invAx requires finite B and x inputs.");
    }

    arma::mat P;
    if(!arma::lu(L, U, P, A) || arma::any(U.diag() == 0.00) ||
       !L.is_finite() || !U.is_finite()) {
      Rcpp::stop("invAx: I-B is singular or its LU factorization failed.");
    }

    // Armadillo uses P*A = L*U. Apply P and P.t() by indexing.
    permutation = arma::index_max(P, 1);
    inverse_permutation = arma::sort_index(permutation);
    Lt = L.t();
    Ut = U.t();
    y = solve_A(v);

  }

  arma::vec solve_A(const arma::vec& rhs, bool transpose = false) const {

    arma::vec intermediate, result;
    const auto opts = arma::solve_opts::fast+arma::solve_opts::no_approx;
    bool ok;

    if(!transpose) {
      ok = arma::solve(intermediate, arma::trimatl(L), rhs.elem(permutation), opts);
      if(ok) {
        ok = arma::solve(result, arma::trimatu(U), intermediate, opts);
      }
    } else {
      ok = arma::solve(intermediate, arma::trimatl(Ut), rhs, opts);
      if(ok) {
        ok = arma::solve(result, arma::trimatu(Lt), intermediate, opts);
        if(ok) {
          result = result.elem(inverse_permutation);
        }
      }
    }

    if(!ok || !result.is_finite()) {
      Rcpp::stop("invAx: triangular solve failed or returned nonfinite values.");
    }

    return result;

  }

};

invAx* choose_invAx(const Rcpp::List& trans_setup) {

  std::vector<arma::uvec> indices_in = trans_setup["indices_in"];
  std::vector<arma::uvec> indices_out = trans_setup["indices_out"];
  int p = trans_setup["p"];

  if(p < 1) {
    Rcpp::stop("invAx requires a positive p dimension.");
  }
  if(indices_in.size() != 2L || indices_out.size() != 1L) {
    Rcpp::stop("invAx requires B and x inputs and one output.");
  }

  const arma::uword n = static_cast<arma::uword>(p);
  if(indices_in[0].n_elem != n*n || indices_in[1].n_elem != n ||
     indices_out[0].n_elem != n) {
    Rcpp::stop("The invAx parameter indices have incompatible dimensions.");
  }

  // Usually all output indices are distinct and all weights are one.
  // Shared output indices must represent structurally equal entries.
  std::unordered_map<arma::uword, arma::uword> counts;
  for(arma::uword index : indices_out[0]) {
    ++counts[index];
  }
  arma::vec output_weights(n);
  for(arma::uword k = 0; k < n; ++k) {
    output_weights[k] = 1.00/static_cast<double>(counts[indices_out[0][k]]);
  }

  invAx* mytrans = new invAx();
  mytrans->indices_in = arma::join_cols(indices_in[0], indices_in[1]);
  mytrans->indices_out = indices_out[0];
  mytrans->indices_B = indices_in[0];
  mytrans->indices_x = indices_in[1];
  mytrans->indices_invAx = indices_out[0];
  mytrans->output_weights = output_weights;
  mytrans->p = p;

  return mytrans;

}

#endif
