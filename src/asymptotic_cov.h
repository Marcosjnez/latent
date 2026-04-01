/*
 * Author: Marcos Jimenez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 01/04/2026
 *
 */

arma::vec diagcov(arma::mat X) {

  // Get the variance of variables with missing data
  const size_t numCols = X.n_cols;

  // Initialize the correlation matrix
  arma::vec diag_cov(numCols);

  // Loop over all pairs of columns
  for (size_t i = 0; i < numCols; ++i) {
    // Get the columns for the pair (i, j)
    arma::vec col = X.col(i);

    // Find indices where both columns have non-NaN values
    arma::uvec validIndices = arma::find_finite(col);

    // Extract non-NaN values from both columns
    arma::vec validCol = col(validIndices);
    validCol -= arma::mean(validCol);

    // Calculate the correlation between the two columns
    double d = arma::accu(validCol % validCol) / (validIndices.n_elem-1);

    // Assign the correlation value to the correlation matrix
    diag_cov(i) = d;
  }

  return diag_cov;

}

arma::mat asymptotic_general(arma::mat X) {

  /*
   * Browne and Shapiro (Equation 3.2; 1986)
   */

  // X is the raw data of scores

  arma::vec d;
  arma::mat P;

  // Compute the standard deviations and correlation matrix of X:
  if(X.has_nan()) {
    d = arma::sqrt(diagcov(X));
    P = pairwise_cor(X);
  } else {
    arma::mat colmeans = arma::mean(X, 0);
    X.each_row() -= colmeans; // Centered matrix
    arma::mat S = arma::cov(X, 1); // Covariance matrix
    d = arma::sqrt(arma::diagvec(S));
    arma::mat diag_d_inv = arma::diagmat(1/d);
    P = diag_d_inv * S * diag_d_inv; // Correlation matrix
  }

  arma::vec p = arma::vectorise(P);
  int q = X.n_cols;
  int qq = q*q;
  arma::mat Theta(qq, qq); // Fourth-order moments

  int ij = 0;
  int kh;

  for(int j=0; j < q; ++j) {
    for(int i=0; i < q; ++i) {
      kh = 0;
      for(int h=0; h < q; ++h) {
        for(int k=0; k < q; ++k) {
          arma::vec m = X(arma::span::all, i) % X(arma::span::all, j) %
          X(arma::span::all, k) % X(arma::span::all, h);
          // Find indices where the vector has non-NaN values:
          arma::uvec validIndices = arma::find_finite(m);
          // Extract non-NaN values from the vector:
          arma::vec v = m(validIndices);
          // m.replace(arma::datum::nan, 0);  // replace each NaN with 0
          Theta(ij, kh) = arma::mean(v) / (d[i]*d[j]*d[k]*d[h]);
          ++kh;
        }
      }
      ++ij;
    }
  }

  arma::mat Gamma = Theta - p * p.t();

  arma::mat Ms = dxt(q, q)*0.5;
  Ms.diag() += 0.5;
  arma::mat I(q, q, arma::fill::eye);
  arma::mat Kd(qq, q, arma::fill::zeros);
  for(int i=0; i < q; ++i) {
    int ii = i * q + i;
    Kd(ii, i) = 1;
  }
  arma::mat A = Ms * arma::kron(I, P) * Kd;
  arma::mat B = Gamma * Kd;
  arma::mat G = Kd.t() * Gamma * Kd;

  arma::mat asymptotic = Gamma - A*B.t() - B*A.t() + A*G*A.t();

  arma::uvec lower_indices = arma::trimatl_ind(arma::size(P), -1);

  return asymptotic(lower_indices, lower_indices);

}

arma::mat asymptotic_normal(arma::mat P) {

  /*
   * Browne and Shapiro (Equation 4.1; 1986)
   */

  // P is the correlation matrix

  int q = P.n_rows;
  int qq = q*q;

  arma::mat Ms = dxt(q, q)*0.5;
  Ms.diag() += 0.5;
  arma::mat I(q, q, arma::fill::eye);
  arma::mat Kd(qq, q, arma::fill::zeros);
  for(int i=0; i < q; ++i) {
    int ii = i * q + i;
    Kd(ii, i) = 1;
  }
  arma::mat A = Ms * arma::kron(I, P) * Kd;
  arma::mat Gamma = 2*Ms * arma::kron(P, P);
  arma::mat B = Gamma * Kd;
  arma::mat G = 2*P % P; // Cheaper than Kd.t() * Gamma * Kd

  arma::mat asymptotic = Gamma - A*B.t() - B*A.t() + A*G*A.t();

  arma::uvec lower_indices = arma::trimatl_ind(arma::size(P), -1);

  return asymptotic(lower_indices, lower_indices);

}

arma::mat asymptotic_elliptical(arma::mat P, double eta) {

  /*
   * Browne and Shapiro (Equation 4.2; 1986)
   */

  arma::mat asymptotic = eta * asymptotic_normal(P);

  return asymptotic;

}
