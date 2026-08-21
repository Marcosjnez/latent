/*
 * Author: Marcos Jimenez
 * Modification date: 21/08/2026
 *
 * Numerical helpers for symmetric positive-definite matrices.
 */

#ifndef LATENT_SYMMETRIC_MATRIX_UTILS_H
#define LATENT_SYMMETRIC_MATRIX_UTILS_H

inline arma::mat make_sympd_latent(const arma::mat& X,
                                   double minimum_eigenvalue = 1e-05) {

  if(X.n_rows != X.n_cols ||
     X.n_rows == 0L ||
     !X.is_finite() ||
     !std::isfinite(minimum_eigenvalue) ||
     minimum_eigenvalue <= 0.0) {
    throw std::runtime_error("Cannot regularize an invalid covariance matrix.");
  }

  arma::mat result = 0.5*(X+X.t());

  if(result.is_sympd()) {
    return result;
  }

  arma::vec eigenvalues;
  arma::mat eigenvectors;

  bool success =
    arma::eig_sym(eigenvalues, eigenvectors, result);

  if(!success ||
     eigenvalues.n_elem != result.n_rows ||
     !eigenvalues.is_finite() ||
     !eigenvectors.is_finite()) {
    throw std::runtime_error("The covariance eigendecomposition failed.");
  }

  eigenvalues.transform(
    [minimum_eigenvalue](double value) {
      return std::max(value, minimum_eigenvalue);
    }
  );

  result =
    eigenvectors*arma::diagmat(eigenvalues)*eigenvectors.t();
  result = 0.5*(result+result.t());

  if(!result.is_sympd()) {
    throw std::runtime_error("The covariance matrix could not be made positive definite.");
  }

  return result;

}

#endif
