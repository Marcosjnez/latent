/*
 * Author: Marcos Jimenez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 18/03/2022
 *
 */

// #include "auxiliary_manifolds.h"

// [[Rcpp::export]]
arma::mat asymptotic_general(arma::mat X) {

  /*
   * Browne and Shapiro (Equation 3.2; 1986)
   */

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

  return asymptotic;

}

// [[Rcpp::export]]
arma::mat asymptotic_normal(arma::mat P) {

  /*
   * Browne and Shapiro (Equation 4.1; 1986)
   */

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

  return asymptotic;

}

// [[Rcpp::export]]
arma::mat asymptotic_elliptical(arma::mat P, double eta) {

  /*
   * Browne and Shapiro (Equation 4.2; 1986)
   */

  arma::mat asymptotic = eta * asymptotic_normal(P);

  return asymptotic;

}

// [[Rcpp::export]]
arma::mat asymptotic_poly(const arma::mat X, const arma::mat R, const int cores) {

  /*
   * Function to estimate the Asymptotic covariance matrix of the polychoric correlations
   */

  // Rcpp::Timer timer;
  // Rcpp::List result;

  const int n = X.n_rows;
  const int q = X.n_cols;

  std::vector<std::vector<int>> cols(q);
  std::vector<int> maxs(q);
  std::vector<std::vector<double>> taus(q);
  std::vector<std::vector<double>> mvphi(q);

  omp_set_num_threads(cores);
#pragma omp parallel for
  for(size_t i = 0; i < q; ++i) {
    cols[i] = arma::conv_to<std::vector<int>>::from(X.col(i));
    maxs[i] = *max_element(cols[i].begin(), cols[i].end());
    std::vector<int> frequencies = count(cols[i], n, maxs[i]-1L);
    mvphi[i] = cumsum(frequencies);
    taus[i] = mvphi[i]; // Cumulative frequencies
    for (size_t j = 0; j < maxs[i]; ++j) {
      mvphi[i][j] /= n;
      taus[i][j] = Qnorm(mvphi[i][j]);
    }
    mvphi[i].push_back(1.0);
    mvphi[i].insert(mvphi[i].begin(), 0.0);
    taus[i].push_back(pos_inf);
    taus[i].insert(taus[i].begin(), neg_inf);
  }

  // timer.step("Thresholds");
  // result["taus"] = taus;

  int dq = 0.5*q*(q-1);
  int k = 0;
  std::vector<std::vector<int>> indexes(dq, std::vector<int>(2));
  for(int i=0; i < (q-1); ++i) {
    for(int j=i+1L; j < q; ++j) {
      indexes[k][0] = i;
      indexes[k][1] = j;
      ++k;
    }
  }

  arma::mat ACOV(dq, dq);
  double f = 0.00;

  for(int i=0; i < dq; ++i) {
    int indexes1 = indexes[i][0];
    int indexes2 = indexes[i][1];
    int s = taus[indexes1].size()+1L;
    int r = taus[indexes2].size()+1L;
    double rho1 = R(indexes1, indexes2);
    Rcpp::List deriv = poly_derivatives(rho1, taus[indexes1], taus[indexes2],
                                        mvphi[indexes1], mvphi[indexes2]);
    arma::mat ppi = deriv["ppi"];
    arma::mat dppidp = deriv["dppidp"];
    arma::mat dppidtau1 = deriv["dppidtau1"];
    arma::mat dppidtau2 = deriv["dppidtau2"];
    // result["ppi"] = ppi;
    // result["dppidp"] = dppidp;
    // result["dppidtau1"] = dppidtau1;
    // result["dppidtau2"] = dppidtau2;
    // return result;
    Rcpp::List x1 = COV(rho1, taus[indexes1], taus[indexes2],
                        mvphi[indexes1], mvphi[indexes2], ppi, dppidp,
                        dppidtau1, dppidtau2);
    // return x1;
    arma::mat Gamma1 = x1["Gamma"];
    double omega1 = x1["omega"];

    for(int j=0; j < dq; ++j) {
      int indexes3 = indexes[j][0];
      int indexes4 = indexes[j][1];
      int y = taus[indexes3].size()+1L;
      int w = taus[indexes4].size()+1L;
      double rho2 = R(indexes3, indexes4);
      Rcpp::List deriv = poly_derivatives(rho2, taus[indexes3], taus[indexes4],
                                          mvphi[indexes3], mvphi[indexes4]);
      arma::mat ppi = deriv["ppi"];
      arma::mat dppidp = deriv["dppidp"];
      arma::mat dppidtau3 = deriv["dppidtau1"];
      arma::mat dppidtau4 = deriv["dppidtau2"];
      Rcpp::List x2 = COV(rho2, taus[indexes3], taus[indexes4],
                          mvphi[indexes3], mvphi[indexes4], ppi, dppidp,
                          dppidtau3, dppidtau4);
      arma::mat Gamma2 = x2["Gamma"];
      double omega2 = x2["omega"];
      double omega_prod = omega1*omega2;
      for(int a=0; a < s; ++a) {
        for(int b=0; b < r; ++b) {
          for(int c=0; c < y; ++c) {
            for(int d=0; d < w; ++d) {
              for(int l=0; l < n; ++l) {
                if(X(l, indexes1) == a && X(l, indexes2) == b &&
                   X(l, indexes3) == c && X(l, indexes4) == d) {
                  f += Gamma1(a, b) * Gamma2(d, c) - omega_prod;
                }
              }
            }
          }
        }
      }
      ACOV(i, j) = ACOV(j, i) = f/n;
      f = 0.00;
    }
  }

  // timer.step("ACOV");

  // result["ACOV"] = ACOV;
  // result["elapsed"] = timer;

  return ACOV;

}

arma::mat asymp_cov(arma::mat S,
                    Rcpp::Nullable<arma::mat> nullable_X,
                    double eta, std::string type) {

  arma::mat asymptotic_cov, X;

  if(type == "normal") {
    asymptotic_cov = asymptotic_normal(S);
  } else if(type == "elliptical") {
    asymptotic_cov = asymptotic_elliptical(S, eta);
  } else if(type == "general") {
    if(nullable_X.isNotNull()) {
      X = Rcpp::as<arma::mat>(nullable_X);
    } else {
      Rcpp::stop("The asymptotic covariance matrix of a general correlation matrix requires the raw data");
    }
    asymptotic_cov = asymptotic_general(X);
  } else if(type == "poly") {
    if(nullable_X.isNotNull()) {
      X = Rcpp::as<arma::mat>(nullable_X);
    } else {
      Rcpp::stop("The asymptotic covariance matrix of a general correlation matrix requires the raw data");
    }
    asymptotic_cov = asymptotic_poly(X, S, 1L);
  } else{
    Rcpp::stop("The available asymptotic covariance estimators are 'normal', 'elliptical', 'general', and 'poly'");
  }

  return asymptotic_cov;

}

