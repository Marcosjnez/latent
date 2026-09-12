/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 12/09/2026
 */

#ifndef LATENT_INVA_X_INVAT_H
#define LATENT_INVA_X_INVAT_H

#include <unordered_map>

// S = inv(A)*X*inv(A).t(), with A = I-B.
// Inputs are B and X, in that order; both and the output are p x p.
// Neither input needs to be symmetric. A must be nonsingular.
// One pivoted LU factorization is reused; no full inverse is formed.
// Derivative updates require transform() at the current parameter point.
// update_dgrad() also requires update_grad() and dtransform() there.

class invA_X_invAt: public transformations {

public:

  int p;
  arma::uvec indices_B, indices_X, indices_invA_X_invAt;
  arma::mat A, X, S, dB, dX, dS, grad_out, grad_in_B, grad_in_X;
  arma::mat output_weights;

  void transform(arguments_optim& x) {

    evaluate(x);
    x.transparameters(indices_invA_X_invAt) = arma::vectorise(S);

  }

  void update_grad(arguments_optim& x) {

    grad_out = arma::reshape(x.grad(indices_invA_X_invAt), p, p);
    grad_out %= output_weights;

    grad_in_B = solve_A(grad_out*S.t()+grad_out.t()*S, true);
    grad_in_X = sandwich(grad_out, true);

    // Explicit accumulation also handles repeated/shared input indices.
    for(arma::uword k = 0; k < indices_B.n_elem; ++k) {
      x.grad[indices_B[k]] += grad_in_B[k];
      x.grad[indices_X[k]] += grad_in_X[k];
    }

  }

  void dtransform(arguments_optim& x) {

    dB = arma::reshape(x.dtransparameters(indices_B), p, p);
    dX = arma::reshape(x.dtransparameters(indices_X), p, p);

    arma::mat K = solve_A(dB);
    dS = sandwich(dX)+K*S+S*K.t();
    x.dtransparameters(indices_invA_X_invAt) = arma::vectorise(dS);

  }

  void update_dgrad(arguments_optim& x) {

    arma::mat dgrad_out = arma::reshape(x.dgrad(indices_invA_X_invAt), p, p);
    dgrad_out %= output_weights;

    arma::mat dgrad_in_B = solve_A(dB.t()*grad_in_B+
      dgrad_out*S.t()+grad_out*dS.t()+dgrad_out.t()*S+grad_out.t()*dS, true);

    arma::mat K = solve_A(dB.t(), true).t(); // dB*inv(A)
    arma::mat dgrad_in_X = sandwich(dgrad_out, true)+
      K.t()*grad_in_X+grad_in_X*K;

    for(arma::uword k = 0; k < indices_B.n_elem; ++k) {
      x.dgrad[indices_B[k]] += dgrad_in_B[k];
      x.dgrad[indices_X[k]] += dgrad_in_X[k];
    }

  }

  void jacobian(arguments_optim& x) {

    // Refresh caches without overwriting any parameter or derivative vector.
    evaluate(x);
    const arma::uword n = static_cast<arma::uword>(p);
    const arma::uword n2 = n*n;
    jacob.zeros(n2, 2*n2);

    // Column order: vec(B), vec(X). Solve for individual columns rather
    // than building a full inverse. This dense Jacobian has O(p^4) entries.
    arma::vec e(n, arma::fill::zeros);
    for(arma::uword i = 0; i < n; ++i) {
      e.zeros();
      e[i] = 1.00;
      arma::mat u = solve_A(e);

      for(arma::uword j = 0; j < n; ++j) {
        jacob.col(i+j*n) = arma::vectorise(u*S.row(j)+S.col(j)*u.t());
      }

      for(arma::uword j = i; j < n; ++j) {
        e.zeros();
        e[j] = 1.00;
        arma::mat v = (j == i) ? u : solve_A(e);
        arma::mat column = u*v.t();
        jacob.col(n2+i+j*n) = arma::vectorise(column);
        if(j != i) {
          jacob.col(n2+j+i*n) = arma::vectorise(column.t());
        }
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
    X = arma::reshape(x.transparameters(indices_X), p, p);

    if(!A.is_finite() || !X.is_finite()) {
      Rcpp::stop("invA_X_invAt requires finite B and X inputs.");
    }

    arma::mat P;
    if(!arma::lu(L, U, P, A) || arma::any(U.diag() == 0.00) ||
       !L.is_finite() || !U.is_finite()) {
      Rcpp::stop("invA_X_invAt: I-B is singular or its LU factorization failed.");
    }

    // Armadillo uses P*A = L*U. Row permutations avoid dense P products.
    permutation = arma::index_max(P, 1);
    inverse_permutation = arma::sort_index(permutation);
    Lt = L.t();
    Ut = U.t();
    S = sandwich(X);

  }

  arma::mat solve_A(const arma::mat& rhs, bool transpose = false) const {

    arma::mat intermediate, result;
    const auto opts = arma::solve_opts::fast+arma::solve_opts::no_approx;
    bool ok;

    if(!transpose) {
      ok = arma::solve(intermediate, arma::trimatl(L), rhs.rows(permutation), opts);
      if(ok) {
        ok = arma::solve(result, arma::trimatu(U), intermediate, opts);
      }
    } else {
      ok = arma::solve(intermediate, arma::trimatl(Ut), rhs, opts);
      if(ok) {
        ok = arma::solve(result, arma::trimatu(Lt), intermediate, opts);
        if(ok) {
          result = result.rows(inverse_permutation);
        }
      }
    }

    if(!ok || !result.is_finite()) {
      Rcpp::stop("invA_X_invAt: triangular solve failed or returned nonfinite values.");
    }

    return result;

  }

  arma::mat sandwich(const arma::mat& rhs, bool transpose = false) const {

    arma::mat left = solve_A(rhs, transpose);
    arma::mat result = solve_A(left.t(), transpose).t();
    return result;

  }

};

invA_X_invAt* choose_invA_X_invAt(const Rcpp::List& trans_setup) {

  std::vector<arma::uvec> indices_in = trans_setup["indices_in"];
  std::vector<arma::uvec> indices_out = trans_setup["indices_out"];
  int p = trans_setup["p"];

  if(p < 1) {
    Rcpp::stop("invA_X_invAt requires a positive p dimension.");
  }
  if(indices_in.size() != 2L || indices_out.size() != 1L) {
    Rcpp::stop("invA_X_invAt requires B and X inputs and one output.");
  }

  const arma::uword n2 = static_cast<arma::uword>(p)*p;
  if(indices_in[0].n_elem != n2 || indices_in[1].n_elem != n2 ||
     indices_out[0].n_elem != n2) {
    Rcpp::stop("The invA_X_invAt parameter indices have incompatible dimensions.");
  }

  // Split adjoints only for output entries sharing a parameter index.
  // For symmetric storage this halves the off-diagonal contributions,
  // as in XtYX; for arbitrary full-matrix storage all weights are one.
  // Shared output indices must represent equal entries (e.g., symmetry).
  std::unordered_map<arma::uword, arma::uword> counts;
  for(arma::uword index : indices_out[0]) {
    ++counts[index];
  }
  arma::mat output_weights(p, p);
  for(arma::uword k = 0; k < n2; ++k) {
    output_weights[k] = 1.00/static_cast<double>(counts[indices_out[0][k]]);
  }

  invA_X_invAt* mytrans = new invA_X_invAt();
  mytrans->indices_in = arma::join_cols(indices_in[0], indices_in[1]);
  mytrans->indices_out = indices_out[0];
  mytrans->indices_B = indices_in[0];
  mytrans->indices_X = indices_in[1];
  mytrans->indices_invA_X_invAt = indices_out[0];
  mytrans->output_weights = output_weights;
  mytrans->p = p;

  return mytrans;

}

#endif
