/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 03/02/2025
 */

double root_quad(double a, double b, double c) {

  // Solution to quadratic interpolation:

  double res = 0.5 * (- b + sqrt(b * b - 4 * a * c) ) / a;

  return res;

}

arma::mat dxt(int p, int q) {

  /*
   * derivative of a matrix wrt its transpose
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

arma::mat lyap_sym(arma::mat Y, arma::mat Q) {

  // Solve the lyapunov equation YX + XY = Q with symmetric Q and X:

  int q = Y.n_cols;
  arma::vec I(q, arma::fill::ones);

  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, Y);

  arma::mat M = eigvec.t() * Q * eigvec;
  arma::mat W1 = I * eigval.t();
  arma::mat W = W1 + W1.t();
  arma::mat YY = M / W;
  arma::mat A = eigvec * YY * eigvec.t();

  return A;

}

arma::uvec consecutive(int lower, int upper) {

  // Generate a sequence of integers from lower to upper

  int size = upper - lower + 1;
  arma::uvec ivec(size);
  std::iota(ivec.begin(), ivec.end(), lower);

  return ivec;
}

// Derivatives wrt model correlation

arma::mat gLRhat(arma::mat Lambda, arma::mat Phi) {

  int p = Lambda.n_rows;
  int q = Lambda.n_cols;
  arma::mat I(p, p, arma::fill::eye);
  arma::mat LP = Lambda * Phi;
  arma::mat g1 = arma::kron(LP, I);
  arma::mat g21 = arma::kron(I, LP);
  arma::mat g2 = g21 * dxt(p, q);
  arma::mat g = g1 + g2;

  return g;

}

arma::mat gPRhat(arma::mat Lambda, arma::mat Phi) {

  int q = Phi.n_cols;
  arma::mat g1 = arma::kron(Lambda, Lambda);
  arma::mat g2 = g1 * dxt(q, q);
  arma::mat g_temp = g1 + g2;
  arma::uvec indexes_diag_q(q);
  for(int i=0; i < q; ++i) indexes_diag_q[i] = i * q + i;
  g_temp.cols(indexes_diag_q) *= 0.5;

  return g_temp;
}

arma::mat gURhat(int p) {

  arma::mat gPsi = dxt(p, p);
  gPsi.diag().ones();

  return gPsi;

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

arma::vec soft(arma::vec x, double a) {

  arma::vec x_trunc_exp = arma::trunc_exp(x);
  arma::vec probs = x_trunc_exp / arma::accu(x_trunc_exp) * a;
  return probs;

}

const double SQRT2M_PI = std::sqrt(2 * M_PI);
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

// [[Rcpp::export]]
double logDnorm(double x, double mu, double s2) {
  x -= mu;
  return -0.5*x*x/s2 - arma::trunc_log(SQRT2M_PI*sqrt(s2));
}

