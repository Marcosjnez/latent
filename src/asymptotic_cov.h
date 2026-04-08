/*
 * Author: Marcos Jimenez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 08/04/2026
 *
 */

arma::vec diagcov(arma::mat X) {

  const size_t numCols = X.n_cols;
  arma::vec diag_cov(numCols);

  for (size_t i = 0; i < numCols; ++i) {
    arma::vec col = X.col(i);
    arma::uvec validIndices = arma::find_finite(col);
    arma::vec validCol = col(validIndices);
    validCol -= arma::mean(validCol);
    double d = arma::accu(validCol % validCol) / (validIndices.n_elem - 1);
    diag_cov(i) = d;
  }

  return diag_cov;
}

arma::mat center_finite(arma::mat X) {

  // Center each column ignoring NaNs

  for (arma::uword j = 0; j < X.n_cols; ++j) {
    arma::vec col = X.col(j);
    arma::uvec idx = arma::find_finite(col);
    if (idx.n_elem > 0) {
      double mu = arma::mean(col(idx));
      col(idx) -= mu;
      X.col(j) = col;
    }
  }
  return X;
}

arma::mat asymptotic_general(arma::mat X, bool cov) {

  /*
   * Browne and Shapiro (Equation 3.2; 1986)
   *
   * cov = false : asymptotic covariance matrix of item correlations
   *               (off-diagonal only)
   * cov = true  : asymptotic covariance matrix of item covariances
   *               (including variances)
   */

  int q = X.n_cols;
  int qq = q * q;

  // Center variables first (important for fourth moments)
  arma::mat Xc = center_finite(X);

  arma::vec d;
  arma::mat S, P;

  if (Xc.has_nan()) {
    d = arma::sqrt(diagcov(Xc));
    P = pairwise_cor(Xc);                 // assumed available in your codebase
    arma::mat D = arma::diagmat(d);
    S = D * P * D;
  } else {
    S = (Xc.t() * Xc) / Xc.n_rows;        // asymptotic scaling; correlation unaffected
    d = arma::sqrt(arma::diagvec(S));
    arma::mat Dinv = arma::diagmat(1.0 / d);
    P = Dinv * S * Dinv;
  }

  arma::mat Theta(qq, qq, arma::fill::zeros);

  int ij = 0;
  int kh = 0;

  for (int j = 0; j < q; ++j) {
    for (int i = 0; i < q; ++i) {

      kh = 0;

      for (int h = 0; h < q; ++h) {
        for (int k = 0; k < q; ++k) {

          arma::vec m = Xc.col(i) % Xc.col(j) % Xc.col(k) % Xc.col(h);
          arma::uvec validIndices = arma::find_finite(m);
          arma::vec v = m(validIndices);

          double val = arma::mean(v);

          if (cov) {
            Theta(ij, kh) = val;
          } else {
            Theta(ij, kh) = val / (d[i] * d[j] * d[k] * d[h]);
          }

          ++kh;
        }
      }

      ++ij;
    }
  }

  if (cov) {
    arma::vec s = arma::vectorise(S);
    arma::mat Gamma = Theta - s * s.t();

    // lower triangle INCLUDING diagonal
    arma::uvec lower_indices = arma::trimatl_ind(arma::size(S), 0);
    return Gamma(lower_indices, lower_indices);
  }

  arma::vec p = arma::vectorise(P);
  arma::mat Gamma = Theta - p * p.t();

  arma::mat Ms = dxt(q, q) * 0.5;
  Ms.diag() += 0.5;

  arma::mat I(q, q, arma::fill::eye);
  arma::mat Kd(qq, q, arma::fill::zeros);

  for (int i = 0; i < q; ++i) {
    int ii = i * q + i;
    Kd(ii, i) = 1.0;
  }

  arma::mat A = Ms * arma::kron(I, P) * Kd;
  arma::mat B = Gamma * Kd;
  arma::mat G = Kd.t() * Gamma * Kd;

  arma::mat asymptotic = Gamma - A * B.t() - B * A.t() + A * G * A.t();

  // lower triangle EXCLUDING diagonal
  arma::uvec lower_indices = arma::trimatl_ind(arma::size(P), 0);
  return asymptotic(lower_indices, lower_indices);
}

arma::mat asymptotic_normal(const arma::mat& S, bool cov) {

  /*
   * Browne and Shapiro (Equation 4.1; 1986)
   *
   * Input is now the covariance matrix S, not the correlation matrix.
   *
   * cov = false : asymptotic covariance matrix of item correlations
   *               (off-diagonal only)
   * cov = true  : asymptotic covariance matrix of item covariances
   *               (including variances)
   */

  int q = S.n_rows;
  int qq = q * q;

  arma::mat Ms = dxt(q, q) * 0.5;
  Ms.diag() += 0.5;

  if (cov) {
    arma::mat Gamma = 2.0 * Ms * arma::kron(S, S);

    // lower triangle INCLUDING diagonal
    arma::uvec lower_indices = arma::trimatl_ind(arma::size(S), 0);
    return Gamma(lower_indices, lower_indices);
  }

  arma::vec d = arma::sqrt(arma::diagvec(S));
  arma::mat Dinv = arma::diagmat(1.0 / d);
  arma::mat P = Dinv * S * Dinv;

  arma::mat I(q, q, arma::fill::eye);
  arma::mat Kd(qq, q, arma::fill::zeros);

  for (int i = 0; i < q; ++i) {
    int ii = i * q + i;
    Kd(ii, i) = 1.0;
  }

  arma::mat A = Ms * arma::kron(I, P) * Kd;
  arma::mat Gamma = 2.0 * Ms * arma::kron(P, P);
  arma::mat B = Gamma * Kd;
  arma::mat G = 2.0 * P % P;   // cheaper than Kd.t() * Gamma * Kd

  arma::mat asymptotic = Gamma - A * B.t() - B * A.t() + A * G * A.t();

  // lower triangle EXCLUDING diagonal
  arma::uvec lower_indices = arma::trimatl_ind(arma::size(P), 0);
  return asymptotic(lower_indices, lower_indices);
}

arma::mat asymptotic_elliptical(const arma::mat& S, double eta, bool cov) {

  /*
   * Browne and Shapiro style elliptical correction.
   *
   * Input is now the covariance matrix S, not the correlation matrix.
   *
   * cov = false : asymptotic covariance matrix of item correlations
   *               (off-diagonal only)
   * cov = true  : asymptotic covariance matrix of item covariances
   *               (including variances)
   *
   * Here eta = 1 + kappa.
   */

  if (!cov) {
    return eta * asymptotic_normal(S, false);
  }

  int q = S.n_rows;
  int qq = q * q;

  arma::mat Ms = dxt(q, q) * 0.5;
  Ms.diag() += 0.5;

  arma::mat Gamma_normal = 2.0 * Ms * arma::kron(S, S);
  arma::vec s = arma::vectorise(S);

  arma::mat Gamma =
    eta * Gamma_normal +
    (eta - 1.0) * (s * s.t());

  // lower triangle INCLUDING diagonal
  arma::uvec lower_indices = arma::trimatl_ind(arma::size(S), 0);
  return Gamma(lower_indices, lower_indices);
}
