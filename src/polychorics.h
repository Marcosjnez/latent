#include <armadillo>
#include <iostream>
#include <array>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <limits>
#include <unordered_map>

// [[Rcpp::export]]
std::vector<int> count(const std::vector<int>& X, const int n, const int max_X) {

  // Allocate memory space for table
  std::vector<int> frequency(max_X + 1);
  // std::vector<int> frequency;
  // frequency.reserve(max_X + 1);

  // Populate table
  for (int i = 0; i < n; i++) {
    // if(!std::isnan(X[i])) frequency[X[i]]++;
    frequency[X[i]]++;
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

  // Return joint frequency table
  return joint_frequency;
  // By row: joint_frequency[0] is row 1 of the frequency table
}

double erf_inverse (double a) {
  double p, r, t;
  t = std::fmaf (a, 0.0f - a, 1.0f);
  t = std::log(t);
  if (fabsf(t) > 6.125f) { // maximum ulp error = 2.35793
    p =              3.03697567e-10f; //  0x1.4deb44p-32
    p = std::fmaf (p, t,  2.93243101e-8f); //  0x1.f7c9aep-26
    p = std::fmaf (p, t,  1.22150334e-6f); //  0x1.47e512p-20
    p = std::fmaf (p, t,  2.84108955e-5f); //  0x1.dca7dep-16
    p = std::fmaf (p, t,  3.93552968e-4f); //  0x1.9cab92p-12
    p = std::fmaf (p, t,  3.02698812e-3f); //  0x1.8cc0dep-9
    p = std::fmaf (p, t,  4.83185798e-3f); //  0x1.3ca920p-8
    p = std::fmaf (p, t, -2.64646143e-1f); // -0x1.0eff66p-2
    p = std::fmaf (p, t,  8.40016484e-1f); //  0x1.ae16a4p-1
  } else { // maximum ulp error = 2.35002
    p =              5.43877832e-9f;  //  0x1.75c000p-28
    p = std::fmaf (p, t,  1.43285448e-7f); //  0x1.33b402p-23
    p = std::fmaf (p, t,  1.22774793e-6f); //  0x1.499232p-20
    p = std::fmaf (p, t,  1.12963626e-7f); //  0x1.e52cd2p-24
    p = std::fmaf (p, t, -5.61530760e-5f); // -0x1.d70bd0p-15
    p = std::fmaf (p, t, -1.47697632e-4f); // -0x1.35be90p-13
    p = std::fmaf (p, t,  2.31468678e-3f); //  0x1.2f6400p-9
    p = std::fmaf (p, t,  1.15392581e-2f); //  0x1.7a1e50p-7
    p = std::fmaf (p, t, -2.32015476e-1f); // -0x1.db2aeep-3
    p = std::fmaf (p, t,  8.86226892e-1f); //  0x1.c5bf88p-1
  }
  r = a * p;
  return r;
}

const double SQRT2M_PI = std::sqrt(2 * M_PI);
double Dnorm(double x) {

  return std::exp(-0.5*x*x) / SQRT2M_PI;

}

double Qnorm(double p) {
  return M_SQRT2 * erf_inverse(2 * p - 1);
}

double MVPHI(double z) {
  return 0.5 * std::erfc(-z * M_SQRT1_2);
}

const double TWOPI = 6.283185307179586;
const std::array<std::array<double, 3>, 10> X = {{
  {{-0.9324695142031522, -0.9815606342467191, -0.9931285991850949}},
  {{-0.6612093864662647, -0.9041172563704750, -0.9639719272779138}},
  {{-0.2386191860831970, -0.7699026741943050, -0.9122344282513259}},
  {{ 0.0,              -0.5873179542866171, -0.8391169718222188}},
  {{ 0.0,              -0.3678314989981802, -0.7463319064601508}},
  {{ 0.0,              -0.1252334085114692, -0.6360536807265150}},
  {{ 0.0,               0.0,              -0.5108670019508271}},
  {{ 0.0,               0.0,              -0.3737060887154196}},
  {{ 0.0,               0.0,              -0.2277858511416451}},
  {{ 0.0,               0.0,              -0.07652652113349733}}
  }};
const std::array<std::array<double, 3>, 10> W = {{
  {{ 0.1713244923791705,  0.04717533638651177, 0.01761400713915212}},
  {{ 0.3607615730481384,  0.1069393259953183,  0.04060142980038694}},
  {{ 0.4679139345726904,  0.1600783285433464,  0.06267204833410906}},
  {{ 0.0,                0.2031674267230659,  0.08327674157670475}},
  {{ 0.0,                0.2334925365383547,  0.1019301198172404}},
  {{ 0.0,                0.2491470458134029,  0.1181945319615184}},
  {{ 0.0,                0.0,                0.1316886384491766}},
  {{ 0.0,                0.0,                0.1420961093183821}},
  {{ 0.0,                0.0,                0.1491729864726037}},
  {{ 0.0,                0.0,                0.1527533871307259}}
  }};

double genz(const double sh, const double sk, const double mvphi_h, const double mvphi_k, const double r) {

  int lg, ng;
  if (std::abs(r) < 0.3) {
    ng = 0;
    lg = 3;
  } else if (std::abs(r) < 0.75) {
    ng = 1;
    lg = 6;
  } else {
    ng = 2;
    lg = 10;
  }

  double h = sh;
  double k = sk;
  double hk = h * k;
  double bvn = 0.0;

  if (std::abs(r) < 0.925) {
    double hs = (h * h + k * k) / 2;
    double asr = std::asin(r);
    for (int i = 0; i < lg; ++i) {
      double sn = std::sin(asr * (X[i][ng] + 1) / 2);
      bvn += W[i][ng] * std::exp((sn * hk - hs) / (1 - sn * sn));
      sn = std::sin(asr * (-X[i][ng] + 1) / 2);
      bvn += W[i][ng] * std::exp((sn * hk - hs) / (1 - sn * sn));
    }
    bvn = bvn * asr / (2 * TWOPI) + (1-mvphi_h) * (1-mvphi_k);
  } else {
    if (r < 0) {
      k = -k;
      hk = -hk;
    }
    if (std::abs(r) < 1) {
      double as = (1 - r) * (1 + r);
      double a = std::sqrt(as);
      double bs = (h - k) * (h - k);
      double c = (4 - hk) / 8;
      double d = (12 - hk) / 16;
      bvn = a * std::exp(-(bs / as + hk) / 2)
        * (1 - c * (bs - as) * (1 - d * bs / 5) / 3 + c * d * as * as / 5);
      if (hk > -160) {
        double b = std::sqrt(bs);
        bvn -= std::exp(-hk / 2) * std::sqrt(TWOPI) * MVPHI(-b / a) * b
          * (1 - c * bs * (1 - d * bs / 5) / 3);
      }
      a = a / 2;
      for (int i = 0; i < lg; ++i) {
        double xs = std::pow(a * (X[i][ng] + 1), 2);
        double rs = std::sqrt(1 - xs);
        bvn += a * W[i][ng]
        * (std::exp(-bs / (2 * xs) - hk / (1 + rs)) / rs
             - std::exp(-(bs / xs + hk) / 2) * (1 + c * xs * (1 + d * xs)));
             xs = as * std::pow(-X[i][ng] + 1, 2) / 4;
             rs = std::sqrt(1 - xs);
             bvn += a * W[i][ng] * std::exp(-(bs / xs + hk) / 2)
               * (std::exp(-hk * (1 - rs) / (2 * (1 + rs))) / rs
                    - (1 + c * xs * (1 + d * xs)));
      }
      bvn = -bvn / TWOPI;
    }
    if (r > 0) {
      if(h > k) {
        bvn += 1-mvphi_h;
      } else {
        bvn += 1-mvphi_k;
      }
    }
    if (r < 0) {
      bvn = -bvn + std::max(0.0, (1-mvphi_h) - mvphi_k);
    }
  }

  return bvn;

}

const int NX = 5L;
const std::vector<double> X2 = {.04691008, .23076534, .5, .76923466, .95308992};
const std::vector<double> W2 = {.018854042, .038088059, .0452707394, .038088059, .018854042};

double drezner(double h1, double hk, const double mvphi_h, const double mvphi_k, double r) {

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
      bv = bv_fac1*h6*ab*(1-MVPHI(h6))-std::exp(-h5/r2)*(ab + aa*r2)*bv_fac2;
      for (int ii=0; ii<NX; ii++){
        r1 = r3*X2[ii];
        rr = r1*r1;
        r2 = std::sqrt( 1.0 - rr);
        bv += - W2[ii]*std::exp(- h5/rr)*(std::exp(-h3/(1.0+r2))/r2/h7 - 1.0 - aa*rr);
      }
    }
    h11 = std::min(h1, h2);
    bv = bv*r3*h7 + MVPHI(h11);
    // if(h1 < h2) {
    //   bv = bv*r3*h7 + mvphi_h;
    // } else {
    //   bv = bv*r3*h7 + mvphi_k;
    // }
    if (r < 0){
      bv = mvphi_h - bv;
    }

  } else {
    h3=h1*h2;
    for (int ii=0; ii<NX; ii++){
      r1 = r*X2[ii];
      rr2 = 1.0 - r1*r1;
      bv += W2[ii] * std::exp(( r1*h3 - h12)/rr2)/ std::sqrt(rr2);
    }
    bv = mvphi_h*mvphi_k + r*bv;
  }
  //--- OUTPUT
  return bv;
}

const double neg_inf = -std::numeric_limits<double>::infinity();
const double pos_inf = std::numeric_limits<double>::infinity();

// [[Rcpp::export]]
double dbinorm(double p, double x, double y) {

  /*
   * Function for the bivariate normal density
   */

  if(!std::isfinite(x) | !std::isfinite(y)) {
    return 0;
  }
  double z1 = x*x + y*y - 2*p*x*y;
  double p2 = p*p;
  double z2 = std::exp(-z1/2/(1-p2));
  double pd = 0.5*z2/M_PI/sqrt(1-p2);

  return pd;
}

double pbinorm(const double lower0, const double lower1, const double upper0, const double upper1, const double rho,
               const double mvphi0, const double mvphi1, const double mvphi2, const double mvphi3) {

  bool ll1 = lower0 == pos_inf;
  bool ll2 = lower1 == pos_inf;
  bool uu1 = upper0 == neg_inf;
  bool uu2 = upper1 == neg_inf;

  if(lower0 > upper0 || lower1 > upper1 || ll1 || ll2 || uu1 || uu2) {
    return 1.0;
  }

  bool l1 = lower0 == neg_inf;
  bool l2 = lower1 == neg_inf;
  bool u1 = upper0 == pos_inf;
  bool u2 = upper1 == pos_inf;

  if(l1) {
    if(l2) {
      if(u1) {
        // return std::normal_distribution<>{0.0, 1.0}(upper1);
        return mvphi3;
      }
      if(u2) {
        // return std::normal_distribution<>{0.0, 1.0}(upper0);
        return mvphi2;
      }
      return genz(-upper0, -upper1, 1.00-mvphi2, 1.00-mvphi3, rho);
    }
    if(u1) {
      if(u2) {
        // return std::normal_distribution<>{0.0, 1.0}(-lower1);
        return 1.00 - mvphi1;
      }
      // return std::normal_distribution<>{0.0, 1.0}(upper1) - std::normal_distribution<>{0.0, 1.0}(lower1);
      return mvphi3 - mvphi1;
    }
    if(u2) {
      return genz(-upper0, lower1, 1.00-mvphi2, mvphi1, -rho);
    }
    return genz(-upper0, -upper1, 1.00-mvphi2, 1.00-mvphi3, rho) -
      genz(-upper0, -lower1, 1.00-mvphi2, 1.00-mvphi1, rho);
  }

  if(u1) {
    if(u2) {
      if(l2) {
        // return std::normal_distribution<>{0.0, 1.0}(-lower0);
        return 1.00 - mvphi0;
      }
      return genz(lower0, lower1, mvphi0, mvphi1, rho);
    }
    if(l2) {
      return genz(-upper1, lower0, 1.00-mvphi3, mvphi0, -rho);
    }
    return genz(lower0, lower1, mvphi0, mvphi1, rho) -
      genz(lower0, upper1, mvphi0, mvphi3, rho);
  }

  if(l2) {
    if(u2) {
      // return std::normal_distribution<>{0.0, 1.0}(upper0) - std::normal_distribution<>{0.0, 1.0}(lower0);
      return mvphi2 - mvphi0;
    }
    return genz(-upper0, -upper1, 1.00-mvphi2, 1.00-mvphi3, rho) -
      genz(-lower0, -upper1, 1.00-mvphi0, 1.00-mvphi3, rho);
  }

  if(u2) {
    return genz(lower0, lower1, mvphi0, mvphi1, rho) -
      genz(upper0, lower1, mvphi2, mvphi1, rho);
  }

  return genz(upper0, upper1, mvphi2, mvphi3, rho) -
    genz(lower0, upper1, mvphi0, mvphi3, rho) -
    genz(upper0, lower1, mvphi2, mvphi1, rho) +
    genz(lower0, lower1, mvphi0, mvphi1, rho);
}

double ddbinorm(const double p, const double x, const double y) {

  /*
   * Function for the derivative of the bivariate normal density
   */

  if(!std::isfinite(x) | !std::isfinite(y)) {
    return 0;
  }

  const double z1 = x*x + y*y - 2*p*x*y;
  const double p2 = p*p;
  const double C = 1-p2;
  const double z2 = std::exp(-0.5*z1/C);
  const double dz1 = -2*x*y;
  const double dp2 = 2*p;
  const double dz2 = -z2 * 0.5*(dz1*C + 2*p*z1)/(C*C);
  const double denom = sqrt(C);
  const double ddenom = 1/(2*denom);
  const double dpd = (1/(2*M_PI)) * (dz2*denom + 2*ddenom*p*z2)/C;

  return dpd;
}

double fpoly(double p, const std::vector<double>& a, const std::vector<double>& b, const std::vector<std::vector<int>>& n,
             const size_t s1, const size_t s2, const std::vector<double>& mvphi1, const std::vector<double>& mvphi2) {

  double f = 0.0;
  for (size_t i = 0; i < s1; ++i) {
    for (size_t j = 0; j < s2; ++j) {
      f -= n[i][j] * std::log(pbinorm(a[i], b[j], a[i + 1], b[j + 1], p,
                                      mvphi1[i], mvphi2[j], mvphi1[i+1], mvphi2[j+1]));
    }
  }

  return f;

}

const double GOLDEN_RATIO = (3.0 - std::sqrt(5.0)) / 2.0;
constexpr double ZEPS = 1.0e-10;

std::vector<double> optimize2(const std::vector<double>& tau1, const std::vector<double>& tau2, const std::vector<std::vector<int>>& n,
                              const size_t s1, const size_t s2, const std::vector<double>& mvphi1, const std::vector<double>& mvphi2,
                              const int nobs, const double cor) {

  // tau1 = Vector of thresholds for the first variable (It must start at -Infinite and end at Infinite)
  // tau2 = Vector of thresholds for the second variable (It must start at -Infinite and end at Infinite)
  // n =  Contingency table for the variables
  // s1 = Length of tau1 - 1L
  // s2 = Length of tau2 - 1L
  // mvphi1 = pnorm of tau1
  // mvphi1 = pnorm of tau2
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
        double prop = pbinorm(tau1[i], tau2[j], tau1[i + 1], tau2[j + 1], p,
                              mvphi1[i], mvphi2[j], mvphi1[i+1], mvphi2[j+1]);
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
                             const size_t s1, const size_t s2, const std::vector<double>& mvphi1, const std::vector<double>& mvphi2,
                             const int nobs, const double cor) {

  // tau1 = Vector of thresholds for the first variable (It must start at -Infinite and end at Infinite)
  // tau2 = Vector of thresholds for the second variable (It must start at -Infinite and end at Infinite)
  // n =  Contingency table for the variables
  // s1 = Length of tau1 - 1L
  // s2 = Length of tau2 - 1L
  // mvphi1 = pnorm of tau1
  // mvphi1 = pnorm of tau2
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
        double prop = pbinorm(tau1[i], tau2[j], tau1[i + 1], tau2[j + 1], p,
                              mvphi1[i], mvphi2[j], mvphi1[i+1], mvphi2[j+1]);
        // PDF of the Bivariate normal:
        double gij = dbinorm(p, tau1[i+1], tau2[j+1]) -
          dbinorm(p, tau1[i], tau2[j+1]) -
          dbinorm(p, tau1[i+1], tau2[j]) +
          dbinorm(p, tau1[i], tau2[j]);
        // f -= n[i][j] * std::log(prop) / nobs; // No need to compute the objective value
        if(prop < 1e-09) prop = 1e-09; // Avoid division by zero
        g -= n[i][j] / prop * gij / nobs; // Update Gradient
        score += gij*gij / prop;          // Update Hessian
      }
    }
    double dir = g/score; // Approximated Newton's Descent direction
    p -= dir;             // Update parameter (no need for step-size)
    if(p > 1 || p < -1) {
      return optimize2(tau1, tau2, n, s1, s2, mvphi1, mvphi2, nobs, cor);
    }
    if((g*g) < 1e-09) break; // Tolerance criteria
    ++ iteration;
  }

  return {p, iteration};

}

std::vector<double> cumsum(const std::vector<int> input) {
  std::vector<double> output(input.size());
  std::partial_sum(input.begin(), input.end(), output.begin());
  return output;
}

Rcpp::List poly_derivatives(double rho, std::vector<double> tau1, std::vector<double> tau2,
                            std::vector<double> mvphi1, std::vector<double> mvphi2) {

  int s = tau1.size()-1L;
  int r = tau2.size()-1L;
  arma::mat ppi(s, r);
  arma::mat dppidp(s, r);
  double denominator = std::sqrt(1-rho*rho);
  for (size_t i = 0; i < s; ++i) {
    for (size_t j = 0; j < r; ++j) {
      // CDF of the bivariate normal:
      ppi(i, j) = pbinorm(tau1[i], tau2[j], tau1[i + 1], tau2[j + 1], rho,
          mvphi1[i], mvphi2[j], mvphi1[i+1], mvphi2[j+1]);
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
      dppidtau1(j*s+k, k) = Dnorm(tau1[k+1])*(MVPHI(numerator1/denominator) -
        MVPHI(numerator2/denominator));
      dppidtau1(j*s+k+1, k) = -dppidtau1(j*s+k, k);
    }
  }
  arma::mat dppidtau2(r*s, r-1, arma::fill::zeros);
  for(int m=0; m < (r-1); ++m) {
    for(int i=0; i < s; ++i) {
      double numerator1 = tau1[i+1]-rho*tau2[m+1];
      double numerator2 = tau1[i]-rho*tau2[m+1];
      dppidtau2(m*s+i, m) = Dnorm(tau2[m+1])*(MVPHI(numerator1/denominator) -
        MVPHI(numerator2/denominator));
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
               std::vector<double> mvphi1, std::vector<double> mvphi2,
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
                std::vector<double> mvphi1, std::vector<double> mvphi2) {

  // Compute the asymptotic variance of the polychoric correlations

  Rcpp::List deriv = poly_derivatives(rho, tau1, tau2, mvphi1, mvphi2);
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
                 std::vector<std::vector<double>> mvphis) {

  // Compute the asymptotic variance of the polychoric correlations

  int p = taus.size();
  arma::mat dacov(p, p);
  int k = 0;
  for(int i=0; i < (p-1); ++i) {
    for(int j=(i+1); j < p; ++j) {
      Rcpp::List x = COV2(poly(i, j), taus[i], taus[j], mvphis[i], mvphis[j]);
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
