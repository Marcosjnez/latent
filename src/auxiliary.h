/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 27/10/2025
 */

double root_quad(double a, double b, double c) {

  // Solution to quadratic interpolation:

  double res = 0.5 * (- b + sqrt(b * b - 4 * a * c) ) / a;

  return res;

}

// [[Rcpp::export]]
arma::mat dxt(int p, int q) {

  /*
   * derivative of a matrix wrt its transpose
   * also a commutation matrix
   */

  int pq = p*q;

  arma::mat res(pq, pq);
  arma::mat temp(p, q);

  for(int i=0; i < pq; ++i) {
    temp.zeros();
    temp(i) = 1;
    res.col(i) = arma::vectorise(temp.t(), 0);
  }

  return res;

}

// [[Rcpp::export]]
arma::mat commutation(const arma::uword m, const arma::uword n) {
  arma::mat K(m * n, m * n, arma::fill::zeros);

  for (arma::uword i = 0; i < m; ++i) {
    for (arma::uword j = 0; j < n; ++j) {
      arma::uword col = j * m + i;      // column index
      arma::uword row = j + n * i;      // row index
      K(row, col) = 1.0;
    }
  }

  return K;
}

arma::mat skew(arma::mat X) {

  // Skew-symmetric matrix

  return 0.5 * (X - X.t());

}

arma::mat symm(arma::mat X) {

  // Symmetric matrix

  return 0.5 * (X + X.t());

}

arma::vec orthogonalize(arma::mat X, arma::vec x) {

  // Make x orthogonal to every column of X (X must be orthogonal)

  int K = X.n_cols;
  for(int i=0; i < K; ++i) {

    // x -= arma::accu(X.col(i) % x) / arma::accu(X.col(i) % X.col(i)) * X.col(i);
    x -= arma::accu(X.col(i) % x) * X.col(i);

  }

  // x /= sqrt(arma::accu(x % x));

  return x;

}

arma::uvec consecutive(int lower, int upper) {

  // Generate a sequence of integers from lower to upper

  int size = upper - lower + 1;
  arma::uvec ivec(size);
  std::iota(ivec.begin(), ivec.end(), lower);

  return ivec;
}

// Derivatives wrt model correlation

arma::mat gLRhat(arma::mat lambda, arma::mat psi) {

  int p = lambda.n_rows;
  int q = lambda.n_cols;
  arma::mat I(p, p, arma::fill::eye);
  arma::mat LP = lambda * psi;
  arma::mat g1 = arma::kron(LP, I);
  arma::mat g21 = arma::kron(I, LP);
  arma::mat g2 = g21 * dxt(p, q);
  arma::mat J = g1 + g2;

  return J;

}

arma::mat gPRhat(arma::mat lambda, int q) {

  arma::mat g1 = arma::kron(lambda, lambda);
  arma::mat g2 = g1 * dxt(q, q);
  arma::mat J = g1 + g2;
  arma::uvec indices_diag_q(q);
  for(int i=0; i < q; ++i) indices_diag_q[i] = i * q + i;
  J.cols(indices_diag_q) *= 0.5;

  return J;
}

arma::mat gURhat(int p) {

  arma::mat J = dxt(p, p);
  J.diag().ones();

  return J;

}

arma::mat cbind_diag(arma::mat X) { // For geomin

  /*
   * Transform every column into a diagonal matrix and bind
   */

  int p = X.n_rows;
  int q = X.n_cols;
  arma::mat res(p, 0);

  for(int i=0; i < q; ++i) {
    res = arma::join_rows(res, arma::diagmat(X.col(i)));
  }

  return res;

}

arma::mat kdiag(arma::mat X) { // For geomin

  /*
   * Transform every column into a diagonal matrix and bind the results
   */

  int pq = X.n_rows;
  int q = X.n_cols;
  int p = pq/q;

  arma::mat res2(pq, 0);

  for(int j=0; j < q; ++j) {

    arma::mat res1(0, p);

    for(int i=0; i < q; ++i) {
      int index_1 = i*p;
      int index_2 = index_1 + (p-1);
      arma::mat temp = arma::diagmat(X(arma::span(index_1, index_2), j));
      res1 = arma::join_cols(res1, temp);
    }

    res2 = arma::join_rows(res2, res1);

  }

  return res2;

}

// [[Rcpp::export]]
arma::vec soft(arma::vec x, double a) {

  arma::vec x_trunc_exp = arma::trunc_exp(x);
  arma::vec probs = x_trunc_exp / arma::accu(x_trunc_exp) * a;
  return probs;

}

const double SQRT2M_PI = std::sqrt(2 * M_PI);
const double logSQRT2M_PI = std::log(SQRT2M_PI);

arma::vec ddnorm2(double x, double mu, double sd) {
  // Derivative of normal density wrt mu and s, where sd = standard deviation
  x -= mu;
  double sd2 = sd*sd;
  double sd3 = sd2*sd;
  double d = arma::trunc_exp(-0.5*x*x/sd2) / (SQRT2M_PI*sd);
  double dmu = x*d/sd2;
  double dsd = (x*x-sd2)*d/sd3;
  arma::vec res = {dmu, dsd};
  res.replace(arma::datum::nan, 0); // Do not remove, it's useful
  return res;
}
arma::vec ddnorm(double x, double mu, double s) {
  // Derivative of lognormal wrt mu and s, where s = log(standard deviation)
  x -= mu;
  double sd = arma::trunc_exp(s);
  double sd2 = sd*sd;
  double d = arma::trunc_exp(-0.5*x*x/sd2) / (SQRT2M_PI*sd);
  double dmu = x*d/sd2;
  double dsd = (x*x-sd2)*d/sd2;
  arma::vec res = {dmu, dsd};
  res.replace(arma::datum::nan, 0); // Do not remove, it's useful
  return res;
}
double logdnorm2(double x, double mu, double sd) {
  x -= mu;
  double sd2 = sd*sd;
  double res = -0.5*x*x/sd2 - arma::trunc_log(SQRT2M_PI*sd);
  if (std::isnan(res)) {
    res = 0.0;
  }
  return res;
}
arma::vec logdnorm(arma::vec x, double mu, double s) {
  x -= mu;
  double sd = arma::trunc_exp(s);
  double sd2 = sd*sd;
  arma::vec res = -0.5* x % x/sd2 - arma::trunc_log(SQRT2M_PI*sd);
  // if (std::isnan(res)) {
  //   res = 0.0;
  // }
  return res;
}
double Dnorm(double x, double mu, double s2) {
  x -= mu;
  return arma::trunc_exp(-0.5*x*x/s2) / (SQRT2M_PI*sqrt(s2));
}

double logDnorm(double x, double mu, double s2) {
  x -= mu;
  return -0.5*x*x/s2 - arma::trunc_log(SQRT2M_PI*sqrt(s2));
}

// [[Rcpp::export]]
arma::mat duplication(int p, bool halflower = true) {

  // Duplication matrix

  int pp = p*p;

  arma::mat res(pp, pp);
  arma::mat null(p, p, arma::fill::zeros);

  for(int i=0; i < pp; ++i) {
    null.zeros();
    null(i) = 1;
    null = arma::symmatl(null);
    res.col(i) = arma::vectorise(null);
  }

  arma::uvec lower_indices = arma::trimatl_ind( arma::size(null) );

  if(halflower) {
    arma::uvec lower = arma::trimatl_ind( arma::size(null), -1 );
    res.cols(lower) *= 0.5;

  }

  return res.cols(lower_indices);

}

// [[Rcpp::export]]
arma::uvec mytest(int p) {

  arma::uvec diag_ind   = arma::regspace<arma::uvec>(0, p - 1) * p
  + arma::regspace<arma::uvec>(0, p - 1);
  return diag_ind;
}

arma::mat pairwise_cor(arma::mat X) {

  const size_t numCols = X.n_cols;

  // Initialize the correlation matrix
  arma::mat corrMatrix(numCols, numCols, arma::fill::eye);

  // Loop over all pairs of columns
  for (size_t i = 0; i < (numCols-1); ++i) {
    for (size_t j = (i+1); j < numCols; ++j) {
      // Get the columns for the pair (i, j)
      arma::vec col1 = X.col(i);
      arma::vec col2 = X.col(j);

      // Find indices where both columns have non-NaN values
      arma::uvec validIndices = arma::find_finite(col1 % col2);

      // Extract non-NaN values from both columns
      arma::vec validCol1 = col1(validIndices);
      arma::vec validCol2 = col2(validIndices);

      // Calculate the correlation between the two columns
      double correlation = as_scalar(arma::cor(validCol1, validCol2));

      // Assign the correlation value to the correlation matrix
      corrMatrix(i, j) = corrMatrix(j, i) = correlation;
    }
  }

  return corrMatrix;
}

// arma::mat center_mat(arma::mat X) {
//   return X.each_row() - arma::mean(X, 0);
// }

// arma::mat bspline(arma::vec x, arma::vec knots, int degree, arma::vec boundaries,
//                   bool center = false, bool intercept = false) {
//
//   int n = x.size();
//   unsigned int p = knots.size() + degree;
//   double z0, z1, output;
//
//   arma::vec lower_boundary(degree);
//   arma::vec upper_boundary(degree);
//   for(int i = 0; i < degree; i++) {
//     lower_boundary(i) = boundaries(0);
//   }
//   for(int i = 0; i < degree; i++) {
//     upper_boundary(i) = boundaries(1);
//   }
//
//   knots = arma::join_cols(lower_boundary, knots, upper_boundary);
//
//   arma::mat X(n, p);
//
//   for(int j = 0; j < n; j++) {
//     for(int k = 0; k < p; k++) {
//       if((x[j] <= knots[k+1]) && (x[j] > knots[k])) {
//         X(j, k) = 1;
//       } else {
//         X(j, k) = 0;
//       }
//     }
//   }
//
//   X.insert_cols(p, 1);
//
//   for(int m = 1; m < (degree+1); m++) {
//     for(int j = 0; j < n; j++) {
//       for(int k = 0; k < p; k++) {
//
//         if(knots[k+m] == knots[k]) {
//           z0 = 0;
//         } else {
//           z0 = (x[j] - knots[k]) / (knots[k+m] - knots[k]);
//         }
//         if(knots[k+m+1] == knots[k+1]) {
//           z1 = 0;
//         } else {
//           z1 = (knots[k+m+1] - x[j]) / (knots[k+m+1] - knots[k+1]);
//         }
//
//         X(j, k) = z0 * X(j, k) + z1 * X(j, k+1);
//       }
//     }
//   }
//
//   arma::uvec indices = {p};
//   X.shed_cols(indices);
//
//   if(center) {
//     X = center_mat(X, n, p);
//   }
//
//   if(intercept) {
//     X.insert_cols(0, arma::ones(n));
//   }
//
//   return X;
// }

// Cox–de Boor recursion for one (x, k, d)
// knots: full knot sequence, k: basis index, d: degree, x: input value
// double bspline_basis(double x, const arma::vec &knots, int k, int d) {
//   if (d == 0) {
//     return (x >= knots[k] && x < knots[k+1]) ? 1.0 : 0.0;
//   } else {
//     double denom1 = knots[k+d] - knots[k];
//     double denom2 = knots[k+d+1] - knots[k+1];
//     double term1 = 0.0;
//     double term2 = 0.0;
//     if (denom1 > 0) {
//       term1 = (x - knots[k]) / denom1 * bspline_basis(x, knots, k, d-1);
//     }
//     if (denom2 > 0) {
//       term2 = (knots[k+d+1] - x) / denom2 * bspline_basis(x, knots, k+1, d-1);
//     }
//     return term1 + term2;
//   }
// }
//
// // [[Rcpp::export]]
// arma::mat bspline(arma::vec x, arma::vec internal_knots, int degree,
//                   arma::vec boundaries, bool center = false, bool intercept = false) {
//
//   int n = x.size();
//   // Knot sequence: degree repeated at each boundary, plus internal knots
//   arma::vec knots = arma::join_cols(
//     arma::vec(degree, arma::fill::value(boundaries[0])),
//     internal_knots,
//     arma::vec(degree, arma::fill::value(boundaries[1]))
//   );
//
//   int K = knots.size() - degree - 1; // Number of basis functions
//
//   arma::mat X(n, K);
//   for (int i = 0; i < n; ++i) {
//     for (int j = 0; j < K; ++j) {
//       X(i, j) = bspline_basis(x[i], knots, j, degree);
//     }
//   }
//
//   if (center) {
//     X = center_mat(X);
//   }
//
//   if (intercept) {
//     X.insert_cols(0, arma::ones(n));
//   }
//
//   return X;
// }
//
