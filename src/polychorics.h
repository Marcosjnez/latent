/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 03/09/2025
 */

// #include <armadillo>
#include <iostream>
#include <array>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <limits>
#include <unordered_map>

// [[Rcpp::export]]
arma::mat real_sqrtmat(arma::mat R) {
  arma::mat X = arma::real(arma::sqrtmat(R));
  return X;
}

// [[Rcpp::export]]
std::vector<int> count(const std::vector<int>& X, const int n, const int max_X) {

  // Allocate memory space for table
  std::vector<int> frequency(max_X);
  // std::vector<int> frequency;
  // frequency.reserve(max_X);

  // Populate table
  for (int i = 0; i < n; i++) {
    if((!std::isnan(X[i])) && (X[i] < max_X)) frequency[X[i]]++;
    // if(X[i] < max_X) frequency[X[i]]++;
  }

  // Return joint frequency table
  return frequency;
}

// [[Rcpp::export]]
std::vector<std::vector<int>> joint_frequency_table(const std::vector<int>& X, const int n, const int max_X,
                                                    const std::vector<int>& Y, const int max_Y) {

  // Allocate memory space for table
  std::vector<std::vector<int>> joint_frequency(max_X + 1L, std::vector<int>(max_Y + 1L, 0));

  // Populate table
  for (int i = 0; i < n; i++) {
    // if(!std::isnan(X[i]) & !std::isnan(Y[i]) ) joint_frequency[X[i]][Y[i]]++;
    joint_frequency[X[i]][Y[i]]++;
  }

  return joint_frequency;
}

// Constants for drezner:
const int NX = 5L;
const std::vector<double> X2 = {.04691008, .23076534, .5, .76923466, .95308992};
const std::vector<double> W2 = {.018854042, .038088059, .0452707394, .038088059, .018854042};

double drezner(double h1, double hk, const double pnorm_tau_h, const double pnorm_tau_k, double r) {

  double bv = 0;
  double r1, r2, rr, rr2, r3, h3, h5, h6, h7, aa, ab, h11;
  double cor_max = 0.7;
  double bv_fac1 = 0.13298076;
  double bv_fac2 = 0.053051647;

  // computation
  double h2 = hk;
  double h12 = (h1*h1+h2*h2)/2;
  double r_abs = std::abs(r);
  if (r_abs > cor_max){
    r2 = 1.0 - r*r;
    r3 = std::sqrt(r2);
    if (r<0){
      h2 = -h2;
    }
    h3 = h1*h2;
    h7 = std::exp( -h3 / 2.0);
    if ( r_abs < 1){
      h6 = std::abs(h1-h2);
      h5 = h6*h6 / 2.0;
      h6 = h6 / r3;
      aa = 0.5 - h3 / 8.0;
      ab = 3.0 - 2.0 * aa * h5;
      bv = bv_fac1*h6*ab*(1-Pnorm(h6))-std::exp(-h5/r2)*(ab + aa*r2)*bv_fac2;
      for (int ii=0; ii<NX; ii++){
        r1 = r3*X2[ii];
        rr = r1*r1;
        r2 = std::sqrt( 1.0 - rr);
        bv += - W2[ii]*std::exp(- h5/rr)*(std::exp(-h3/(1.0+r2))/r2/h7 - 1.0 - aa*rr);
      }
    }
    h11 = std::min(h1, h2);
    bv = bv*r3*h7 + Pnorm(h11);
    // if(h1 < h2) {
    //   bv = bv*r3*h7 + pnorm_tau_h;
    // } else {
    //   bv = bv*r3*h7 + pnorm_tau_k;
    // }
    if (r < 0){
      bv = pnorm_tau_h - bv;
    }

  } else {
    h3=h1*h2;
    for (int ii=0; ii<NX; ii++){
      r1 = r*X2[ii];
      rr2 = 1.0 - r1*r1;
      bv += W2[ii] * std::exp(( r1*h3 - h12)/rr2)/ std::sqrt(rr2);
    }
    bv = pnorm_tau_h*pnorm_tau_k + r*bv;
  }
  //--- OUTPUT
  return bv;
}

// [[Rcpp::export]]
double fpoly(double p, std::vector<double> tau1, std::vector<double> tau2,
             std::vector<double> pnorm_tau1, std::vector<double> pnorm_tau2,
             std::vector<std::vector<int>> n) {

  // tau1 = Vector of thresholds for the first variable (It must start at -Infinite and end at Infinite)
  // tau2 = Vector of thresholds for the second variable (It must start at -Infinite and end at Infinite)
  // n =  Contingency table for the variables
  // pnorm_tau1 = pnorm of tau1
  // pnorm_tau2 = pnorm of tau2

  const size_t s1 = tau1.size()-1L;
  const size_t s2 = tau2.size()-1L;

  double f = 0.0;
  for (size_t i = 0; i < s1; ++i) {
    for (size_t j = 0; j < s2; ++j) {
      f -= n[i][j] * arma::trunc_log(pbinorm(p,
                                             tau1[i], tau2[j], tau1[i + 1], tau2[j + 1],
                                             pnorm_tau1[i], pnorm_tau2[j], pnorm_tau1[i+1], pnorm_tau2[j+1]));
    }
  }

  return f;

}

// [[Rcpp::export]]
double fpoly2(arma::mat R, const std::vector<arma::vec>& tau,
              const std::vector<arma::vec>& pnorm_tau,
              const std::vector<std::vector<std::vector<int>>>& n) {


  const size_t q = tau.size(); // Number of items
  int m = 0L;

  double f = 0.0;
  for(size_t l=0; l < (q-1L); ++l) {
    const size_t s1 = tau[l].size()-1L;
    for(int k=(l+1L); k < q; ++k) {
      const size_t s2 = tau[k].size()-1L;
      for (size_t i = 0; i < s1; ++i) {
        for (size_t j = 0; j < s2; ++j) {
          f -= n[m][i][j] * arma::trunc_log(pbinorm(R(l, k),
                                            tau[l][i], tau[k][j],
                                            tau[l][i+1], tau[k][j+1],
                                            pnorm_tau[l][i], pnorm_tau[k][j],
                                            pnorm_tau[l][i+1], pnorm_tau[k][j+1]));
        }
      }
      ++m;
    }
  }

  return f;

}

const double GOLDEN_RATIO = (3.0 - std::sqrt(5.0)) / 2.0;
constexpr double ZEPS = 1.0e-10;

std::vector<double> optimize2(const std::vector<double>& tau1, const std::vector<double>& tau2, const std::vector<std::vector<int>>& n,
                              const size_t s1, const size_t s2, const std::vector<double>& pnorm_tau1, const std::vector<double>& pnorm_tau2,
                              const int nobs, const double cor) {

  // tau1 = Vector of thresholds for the first variable (It must start at -Infinite and end at Infinite)
  // tau2 = Vector of thresholds for the second variable (It must start at -Infinite and end at Infinite)
  // n =  Contingency table for the variables
  // s1 = Length of tau1 - 1L
  // s2 = Length of tau2 - 1L
  // pnorm_tau1 = pnorm of tau1
  // pnorm_tau2 = pnorm of tau2
  // nobs =  Sample size
  // cor = Initial value for the correlation

  double asin_p = std::asin(cor);
  double p = cor; // Parameters to be estimated
  double cos_asin_p = std::cos(asin_p);
  double iteration = 1;

  // Start the iterative algorithm
  for(int i=0; i < 20L; ++i) {
    // double f = 0.0;  // Objective value (no needed)
    double g = 0.0;     // Gradient
    double h = 0.0; // Approximated Hessian (asymptotic formula)

    for (size_t i = 0; i < s1; ++i) {
      for (size_t j = 0; j < s2; ++j) {
        // CDF of the bivariate normal:
        double prop = pbinorm(p, tau1[i], tau2[j], tau1[i + 1], tau2[j + 1],
                              pnorm_tau1[i], pnorm_tau2[j], pnorm_tau1[i+1], pnorm_tau2[j+1]);
        // PDF of the Bivariate normal:
        double gij = dbinorm(p, tau1[i+1], tau2[j+1]) -
          dbinorm(p, tau1[i], tau2[j+1]) -
          dbinorm(p, tau1[i+1], tau2[j]) +
          dbinorm(p, tau1[i], tau2[j]);
        // Derivative of the PDF of the Bivariate normal:
        double hij = ddbinorm(p, tau1[i+1], tau2[j+1]) -
          ddbinorm(p, tau1[i], tau2[j+1]) -
          ddbinorm(p, tau1[i+1], tau2[j]) +
          ddbinorm(p, tau1[i], tau2[j]);
        // f -= n[i][j] * std::log(prop) / nobs; // No need to compute the objective value
        if(prop < 1e-09) prop = 1e-09; // Avoid division by zero
        double gij_cos = gij*cos_asin_p;
        g -= n[i][j] / prop * gij_cos / nobs; // Update Gradient
        double term = hij*cos_asin_p*cos_asin_p - gij*p;
        h += n[i][j]*(gij_cos*gij_cos - prop*term)/(prop*prop) / nobs; // Update Hessian
      }
    }
    double dir = g/h; // Approximated Newton's Descent direction
    asin_p -= dir;             // Update parameter (no need for step-size)
    p = std::sin(asin_p);
    if((g*g) < 1e-09) break; // Tolerance criteria
    ++ iteration;
    cos_asin_p = std::cos(asin_p);
  }

  return {p, iteration};

}

std::vector<double> optimize(const std::vector<double>& tau1, const std::vector<double>& tau2, const std::vector<std::vector<int>>& n,
                             const size_t s1, const size_t s2, const std::vector<double>& pnorm_tau1, const std::vector<double>& pnorm_tau2,
                             const int nobs, const double cor) {

  // tau1 = Vector of thresholds for the first variable (It must start at -Infinite and end at Infinite)
  // tau2 = Vector of thresholds for the second variable (It must start at -Infinite and end at Infinite)
  // n =  Contingency table for the variables
  // s1 = Length of tau1 - 1L
  // s2 = Length of tau2 - 1L
  // pnorm_tau1 = pnorm of tau1
  // pnorm_tau2 = pnorm of tau2
  // nobs =  Sample size
  // cor = Initial value for the correlation

  double p = cor; // Parameters to be estimated
  double iteration = 1;

  // Start the iterative algorithm
  for(int i=0; i < 20L; ++i) {
    // double f = 0.0;  // Objective value (no needed)
    double g = 0.0;     // Gradient
    double score = 0.0; // Approximated Hessian (asymptotic formula)

    for (size_t i = 0; i < s1; ++i) {
      for (size_t j = 0; j < s2; ++j) {
        // CDF of the bivariate normal:
        double prop = pbinorm(p, tau1[i], tau2[j], tau1[i + 1], tau2[j + 1],
                              pnorm_tau1[i], pnorm_tau2[j], pnorm_tau1[i+1], pnorm_tau2[j+1]);
        // PDF of the Bivariate normal:
        double gij = dbinorm(p, tau1[i+1], tau2[j+1]) -
          dbinorm(p, tau1[i], tau2[j+1]) -
          dbinorm(p, tau1[i+1], tau2[j]) +
          dbinorm(p, tau1[i], tau2[j]);
        // No need to compute the likelihood
        // f -= n[i][j] * std::log(prop) / nobs;
        // A Taylor expansion results in a Fisher scoring method to update the correlation
        if(prop < 1e-09) prop = 1e-09; // Avoid division by zero
        g -= n[i][j] / prop * gij / nobs; // Update Gradient
        score += gij*gij / prop;          // Update Hessian
      }
    }
    double dir = g/score; // Approximated Newton's Descent direction
    p -= dir;             // Update parameter (no need for step-size)
    if(p > 1 || p < -1) {
      return optimize2(tau1, tau2, n, s1, s2, pnorm_tau1, pnorm_tau2, nobs, cor);
    }
    if((g*g) < 1e-09) break; // Tolerance criteria
    ++ iteration;
  }

  return {p, iteration};

}

// [[Rcpp::export]]
Rcpp::List poly_deriv(double rho, std::vector<double> tau1, std::vector<double> tau2,
                      std::vector<double> pnorm_tau1, std::vector<double> pnorm_tau2,
                      std::vector<std::vector<int>> n) {

  // tau1 = Vector of thresholds for the first variable (It must start at -Infinite and end at Infinite)
  // tau2 = Vector of thresholds for the second variable (It must start at -Infinite and end at Infinite)
  // pnorm_tau1 = pnorm of tau1
  // pnorm_tau2 = pnorm of tau2

  int s = tau1.size()-1L;
  int r = tau2.size()-1L;
  arma::mat ppi(s, r);
  arma::mat dppidp(s, r);
  double denominator = std::sqrt(1-rho*rho);
  double dldp = 0.0;
  for (size_t i = 0; i < s; ++i) {
    for (size_t j = 0; j < r; ++j) {
      // CDF of the bivariate normal:
      ppi(i, j) = pbinorm(rho, tau1[i], tau2[j], tau1[i+1], tau2[j+1],
          pnorm_tau1[i], pnorm_tau2[j], pnorm_tau1[i+1], pnorm_tau2[j+1]);
      // PDF of the Bivariate normal:
      dppidp(i, j) = dbinorm(rho, tau1[i+1], tau2[j+1]) -
        dbinorm(rho, tau1[i], tau2[j+1]) -
        dbinorm(rho, tau1[i+1], tau2[j]) +
        dbinorm(rho, tau1[i], tau2[j]);
      dldp -= n[i][j] * dppidp(i,j) / ppi(i,j);
    }
  }

  // arma::mat dppidtau1(s*r, s-1, arma::fill::zeros);
  arma::vec dldtau1(s-1L, arma::fill::zeros);
  for(int k=0; k < (s-1); ++k) {
    for(int j=0; j < r; ++j) {
      double numerator1 = tau2[j+1]-rho*tau1[k+1];
      double numerator2 = tau2[j]-rho*tau1[k+1];
      // dppidtau1(j*s+k, k) = Dnorm(tau1[k+1])*(Pnorm(numerator1/denominator) -
      //   Pnorm(numerator2/denominator));
      // dppidtau1(j*s+k+1, k) = -dppidtau1(j*s+k, k);
      double dppidtau1 = Dnorm(tau1[k+1])*(Pnorm(numerator1/denominator) -
                               Pnorm(numerator2/denominator));
      dldtau1(k) += dppidtau1 * (n[k][j]/ppi(k,j)-n[k+1][j]/ppi(k+1,j));
    }
  }

  // arma::mat dppidtau2(r*s, r-1, arma::fill::zeros);
  arma::vec dldtau2(r-1L, arma::fill::zeros);
  for(int m=0; m < (r-1); ++m) {
    for(int i=0; i < s; ++i) {
      double numerator1 = tau1[i+1]-rho*tau2[m+1];
      double numerator2 = tau1[i]-rho*tau2[m+1];
      // dppidtau2(m*s+i, m) = Dnorm(tau2[m+1])*(Pnorm(numerator1/denominator) -
      //   Pnorm(numerator2/denominator));
      // dppidtau2(m*s+i+s, m) = -dppidtau2(m*s+i, m);
      double dppidtau2 = Dnorm(tau2[m+1])*(Pnorm(numerator1/denominator) -
                               Pnorm(numerator2/denominator));
      dldtau2(m) += dppidtau2 * (n[i][m]/ppi(i,m)-n[i][m+1]/ppi(i,m+1));
    }
  }

  Rcpp::List result;
  result["ppi"] = ppi;
  result["dppidp"] = dppidp;
  // result["dppidtau1"] = dppidtau1;
  // result["dppidtau2"] = dppidtau2;
  result["dldp"] = dldp;
  result["dldtau1"] = -dldtau1;
  result["dldtau2"] = -dldtau2;

  return result;
}

Rcpp::List poly_derivatives(double rho, std::vector<double> tau1, std::vector<double> tau2,
                            std::vector<double> pnorm_tau1, std::vector<double> pnorm_tau2) {

  // tau1 = Vector of thresholds for the first variable (It must start at -Infinite and end at Infinite)
  // tau2 = Vector of thresholds for the second variable (It must start at -Infinite and end at Infinite)
  // pnorm_tau1 = pnorm of tau1
  // pnorm_tau2 = pnorm of tau2

  int s = tau1.size()-1L;
  int r = tau2.size()-1L;
  arma::mat ppi(s, r);
  arma::mat dppidp(s, r);
  double denominator = std::sqrt(1-rho*rho);
  for (size_t i = 0; i < s; ++i) {
    for (size_t j = 0; j < r; ++j) {
      // CDF of the bivariate normal:
      ppi(i, j) = pbinorm(rho, tau1[i], tau2[j], tau1[i + 1], tau2[j + 1],
          pnorm_tau1[i], pnorm_tau2[j], pnorm_tau1[i+1], pnorm_tau2[j+1]);
      // PDF of the Bivariate normal:
      dppidp(i, j) = dbinorm(rho, tau1[i+1], tau2[j+1]) -
        dbinorm(rho, tau1[i], tau2[j+1]) -
        dbinorm(rho, tau1[i+1], tau2[j]) +
        dbinorm(rho, tau1[i], tau2[j]);
    }
  }

  arma::mat dppidtau1(s*r, s-1, arma::fill::zeros);
  for(int k=0; k < (s-1); ++k) {
    for(int j=0; j < r; ++j) {
      double numerator1 = tau2[j+1]-rho*tau1[k+1];
      double numerator2 = tau2[j]-rho*tau1[k+1];
      dppidtau1(j*s+k, k) = Dnorm(tau1[k+1])*(Pnorm(numerator1/denominator) -
        Pnorm(numerator2/denominator));
      dppidtau1(j*s+k+1, k) = -dppidtau1(j*s+k, k);
    }
  }
  arma::mat dppidtau2(r*s, r-1, arma::fill::zeros);
  for(int m=0; m < (r-1); ++m) {
    for(int i=0; i < s; ++i) {
      double numerator1 = tau1[i+1]-rho*tau2[m+1];
      double numerator2 = tau1[i]-rho*tau2[m+1];
      dppidtau2(m*s+i, m) = Dnorm(tau2[m+1])*(Pnorm(numerator1/denominator) -
        Pnorm(numerator2/denominator));
      dppidtau2(m*s+i+s, m) = -dppidtau2(m*s+i, m);
    }
  }

  Rcpp::List result;
  result["ppi"] = ppi;
  result["dppidp"] = dppidp;
  result["dppidtau1"] = dppidtau1;
  result["dppidtau2"] = dppidtau2;

  return result;
}

Rcpp::List COV(double rho, std::vector<double> tau1, std::vector<double> tau2,
               std::vector<double> pnorm_tau1, std::vector<double> pnorm_tau2,
               arma::mat ppi, arma::mat dppidp, arma::mat dppidtau1, arma::mat dppidtau2) {

  tau1.erase(tau1.begin());  // remove the first element
  tau1.erase(tau1.end() - 1);  // remove the last element
  tau2.erase(tau2.begin());  // remove the first element
  tau2.erase(tau2.end() - 1);  // remove the last element

  int s = tau1.size() + 1L;
  arma::mat Ag(s, s-1L, arma::fill::zeros);
  arma::vec dnorm_tau1(s-1L);
  for(int i=0; i < (s-1L); ++i) dnorm_tau1(i) = Dnorm(tau1[i]);
  Ag.diag() = dnorm_tau1;
  Ag.diag(-1) = -dnorm_tau1;
  arma::mat Dg = arma::diagmat(1/arma::sum(ppi, 1));
  arma::mat Bg = (arma::inv(Ag.t() * Dg * Ag) * Ag.t() * Dg).t();

  int r = tau2.size() + 1L;
  arma::mat Ah(r, r-1L, arma::fill::zeros);
  arma::vec dnorm_tau2(r-1L);
  for(int i=0; i < (r-1L); ++i) dnorm_tau2(i) = Dnorm(tau2[i]);
  Ah.diag() = dnorm_tau2;
  Ah.diag(-1) = -dnorm_tau2;
  arma::mat Dh = arma::diagmat(1/arma::sum(ppi, 0));
  arma::mat Bh = (arma::inv(Ah.t() * Dh * Ah) * Ah.t() * Dh).t();

  double D = arma::accu(dppidp % dppidp / ppi);
  arma::mat alpha = (1/D) * (1/ppi) % dppidp;

  // Set alpha cells to zero?

  arma::vec Betag(s-1L);
  for(int i=0; i < (s-1L); ++i) {
    Betag[i] = (1/D) * arma::accu(1/arma::vectorise(ppi) % arma::vectorise(dppidp) % dppidtau1.col(i));
  }
  arma::vec Betah(r-1L);
  for(int i=0; i < (r-1L); ++i) {
    Betah[i] = (1/D) * arma::accu(1/arma::vectorise(ppi) % arma::vectorise(dppidp) % dppidtau2.col(i));
  }
  arma::rowvec ones_r(r, arma::fill::ones);
  arma::vec ones_s(s, arma::fill::ones);
  // Rcpp::List result;
  // result["ppi"] = ppi;
  // result["dppidp"] = dppidp;
  // result["D"] = D;
  // result["tau1"] = tau1;
  // result["tau2"] = tau2;
  // result["dnorm_tau1"] = dnorm_tau1;
  // result["dnorm_tau2"] = dnorm_tau2;
  // result["alpha"] = alpha;
  // result["Bg"] = Bg;
  // result["Betag"] = Betag;
  // result["ones_r"] = ones_r;
  // result["dppidtau1"] = dppidtau1;
  // result["dppidtau2"] = dppidtau2;
  // result["alpha"] = alpha;
  // return result;
  arma::mat Gamma = alpha + Bg * Betag * ones_r + ones_s * Betah.t() * Bh.t();
  double omega = arma::accu(Gamma % ppi);

  Rcpp::List result;
  result["Gamma"] = Gamma;
  result["omega"] = omega;

  return result;

}

Rcpp::List COV2(double rho,
                std::vector<double> tau1, std::vector<double> tau2,
                std::vector<double> pnorm_tau1, std::vector<double> pnorm_tau2) {

  // Compute the asymptotic variance of the polychoric correlations

  Rcpp::List deriv = poly_derivatives(rho, tau1, tau2, pnorm_tau1, pnorm_tau2);
  arma::mat ppi = deriv["ppi"];
  arma::mat dppidp = deriv["dppidp"];
  arma::mat dppidtau1 = deriv["dppidtau1"];
  arma::mat dppidtau2 = deriv["dppidtau2"];

  tau1.erase(tau1.begin());  // remove the first element
  tau1.erase(tau1.end() - 1);  // remove the last element
  tau2.erase(tau2.begin());  // remove the first element
  tau2.erase(tau2.end() - 1);  // remove the last element

  int s = tau1.size() + 1L;
  arma::mat Ag(s, s-1L, arma::fill::zeros);
  arma::vec dnorm_tau1(s-1L);
  for(int i=0; i < (s-1L); ++i) dnorm_tau1(i) = Dnorm(tau1[i]);
  Ag.diag() = dnorm_tau1;
  Ag.diag(-1) = -dnorm_tau1;
  arma::mat Dg = arma::diagmat(1/arma::sum(ppi, 1));
  arma::mat Bg = (arma::inv(Ag.t() * Dg * Ag) * Ag.t() * Dg).t();

  int r = tau2.size() + 1L;
  arma::mat Ah(r, r-1L, arma::fill::zeros);
  arma::vec dnorm_tau2(r-1L);
  for(int i=0; i < (r-1L); ++i) dnorm_tau2(i) = Dnorm(tau2[i]);
  Ah.diag() = dnorm_tau2;
  Ah.diag(-1) = -dnorm_tau2;
  arma::mat Dh = arma::diagmat(1/arma::sum(ppi, 0));
  arma::mat Bh = (arma::inv(Ah.t() * Dh * Ah) * Ah.t() * Dh).t();

  double D = arma::accu(dppidp % dppidp / ppi);
  arma::mat alpha = (1/D) * (1/ppi) % dppidp;

  // Set alpha cells to zero?

  arma::vec Betag(s-1L);
  for(int i=0; i < (s-1L); ++i) {
    Betag[i] = (1/D) * arma::accu(1/arma::vectorise(ppi) % arma::vectorise(dppidp) % dppidtau1.col(i));
  }
  arma::vec Betah(r-1L);
  for(int i=0; i < (r-1L); ++i) {
    Betah[i] = (1/D) * arma::accu(1/arma::vectorise(ppi) % arma::vectorise(dppidp) % dppidtau2.col(i));
  }
  arma::rowvec ones_r(r, arma::fill::ones);
  arma::vec ones_s(s, arma::fill::ones);

  arma::mat Gamma = alpha + Bg * Betag * ones_r + ones_s * Betah.t() * Bh.t();
  double omega = arma::accu(Gamma % ppi);

  Rcpp::List result;
  result["Gamma"] = Gamma;
  result["omega"] = omega;

  return result;

}

arma::mat std_2_matrix(std::vector<std::vector<int>> tabs, int n) {

  int p = tabs.size();
  int q = tabs[0].size();
  arma::mat matrix(p, q);

  for(int i=0; i < p; ++i) { // Fill rows
    for(int j=0; j < q; ++j) {
      matrix(i, j) = (tabs[i][j] + 0.0) / (n + 0.0);
    }
  }

  return matrix;

}

arma::mat DACOV2(int n, arma::mat poly,
                 std::vector<std::vector<std::vector<int>>> tabs,
                 std::vector<std::vector<double>> taus,
                 std::vector<std::vector<double>> pnorm_taus) {

  // Compute the asymptotic variance of the polychoric correlations

  int p = taus.size();
  arma::mat dacov(p, p);
  int k = 0;
  for(int i=0; i < (p-1); ++i) {
    for(int j=(i+1); j < p; ++j) {
      Rcpp::List x = COV2(poly(i, j), taus[i], taus[j], pnorm_taus[i], pnorm_taus[j]);
      arma::mat Gamma = x["Gamma"];
      double omega = x["omega"];
      arma::mat table = std_2_matrix(tabs[k], n);
      dacov(i, j) = arma::accu(Gamma % table % Gamma.t() - omega*omega); // Use Gamma.t() or Gamma???
      dacov(j, i) = dacov(i, j);
      ++k;
    }
  }

  return dacov;

}
