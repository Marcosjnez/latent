/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 13/11/2025
 */

/*
 * Polychoric correlations (estimating the thresholds and correlations)
 */

double erf_inverse(double a) {
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

double dDnorm(double x) {

  // Derivative of the probability density of the normal distribution
  // const double SQRT2M_PI = std::sqrt(2 * M_PI);
  return -x*std::exp(-0.5*x*x) / SQRT2M_PI;

}

double Dnorm(double x) {

  // Probability density of the normal distribution
  // const double SQRT2M_PI = std::sqrt(2 * M_PI);
  return std::exp(-0.5*x*x) / SQRT2M_PI;

}

double Qnorm(double p) {

  // Quantile function of the normal distribution
  return M_SQRT2 * erf_inverse(2*p-1);

}

double Pnorm(double z) {

 // Cumulative probability of the normal distribution
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

double genz(const double sh, const double sk,
            const double pnorm_tau_h, const double pnorm_tau_k, const double r) {

  // Function to compute cummulative bivariate normal probabilities

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
    bvn = bvn * asr / (2 * TWOPI) + (1-pnorm_tau_h) * (1-pnorm_tau_k);
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
        bvn -= std::exp(-hk / 2) * std::sqrt(TWOPI) * Pnorm(-b / a) * b
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
        bvn += 1-pnorm_tau_h;
      } else {
        bvn += 1-pnorm_tau_k;
      }
    }
    if (r < 0) {
      bvn = -bvn + std::max(0.0, (1-pnorm_tau_h) - pnorm_tau_k);
    }
  }

  return bvn;

}

double ddbinorm(const double p, const double x, const double y) {

  /*
   * Function for the derivative of the bivariate normal density
   */

  if(!std::isfinite(x) || !std::isfinite(y)) {
    return 0;
  }

  const double z1 = x*x + y*y - 2*p*x*y;
  const double p2 = p*p;
  const double C = 1-p2;
  const double z2 = std::exp(-0.5*z1/C);
  const double dz1 = -2*x*y;
  const double dz2 = -z2 * 0.5*(dz1*C + 2*p*z1)/(C*C);
  const double denom = sqrt(C);
  const double ddenom = 1/(2*denom);
  const double dpd = (1/(2*M_PI)) * (dz2*denom + 2*ddenom*p*z2)/C;

  return dpd;
}

double dbinorm(double p, double x, double y) {

  /*
   * Function for the bivariate normal density
   */

  if(!std::isfinite(x) || !std::isfinite(y)) {
    return 0;
  }
  double z1 = x*x + y*y - 2*p*x*y;
  double p2 = p*p;
  double z2 = std::exp(-z1/2/(1-p2));
  double pd = 0.5*z2/M_PI/sqrt(1-p2);

  return pd;
}

const double neg_inf = -std::numeric_limits<double>::infinity();
const double pos_inf = std::numeric_limits<double>::infinity();
double pbinorm(const double rho,
               const double lower0, const double lower1,
               const double upper0, const double upper1,
               const double pnorm_lower0, const double pnorm_lower1,
               const double pnorm_upper0, const double pnorm_upper1) {

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
        return pnorm_upper1;
      }
      if(u2) {
        // return std::normal_distribution<>{0.0, 1.0}(upper0);
        return pnorm_upper0;
      }
      return genz(-upper0, -upper1, 1.00-pnorm_upper0, 1.00-pnorm_upper1, rho);
    }
    if(u1) {
      if(u2) {
        // return std::normal_distribution<>{0.0, 1.0}(-lower1);
        return 1.00 - pnorm_lower1;
      }
      // return std::normal_distribution<>{0.0, 1.0}(upper1) - std::normal_distribution<>{0.0, 1.0}(lower1);
      return pnorm_upper1 - pnorm_lower1;
    }
    if(u2) {
      return genz(-upper0, lower1, 1.00-pnorm_upper0, pnorm_lower1, -rho);
    }
    return genz(-upper0, -upper1, 1.00-pnorm_upper0, 1.00-pnorm_upper1, rho) -
      genz(-upper0, -lower1, 1.00-pnorm_upper0, 1.00-pnorm_lower1, rho);
  }

  if(u1) {
    if(u2) {
      if(l2) {
        // return std::normal_distribution<>{0.0, 1.0}(-lower0);
        return 1.00 - pnorm_lower0;
      }
      return genz(lower0, lower1, pnorm_lower0, pnorm_lower1, rho);
    }
    if(l2) {
      return genz(-upper1, lower0, 1.00-pnorm_upper1, pnorm_lower0, -rho);
    }
    return genz(lower0, lower1, pnorm_lower0, pnorm_lower1, rho) -
      genz(lower0, upper1, pnorm_lower0, pnorm_upper1, rho);
  }

  if(l2) {
    if(u2) {
      // return std::normal_distribution<>{0.0, 1.0}(upper0) - std::normal_distribution<>{0.0, 1.0}(lower0);
      return pnorm_upper0 - pnorm_lower0;
    }
    return genz(-upper0, -upper1, 1.00-pnorm_upper0, 1.00-pnorm_upper1, rho) -
      genz(-lower0, -upper1, 1.00-pnorm_lower0, 1.00-pnorm_upper1, rho);
  }

  if(u2) {
    return genz(lower0, lower1, pnorm_lower0, pnorm_lower1, rho) -
      genz(upper0, lower1, pnorm_upper0, pnorm_lower1, rho);
  }

  return genz(upper0, upper1, pnorm_upper0, pnorm_upper1, rho) -
    genz(lower0, upper1, pnorm_lower0, pnorm_upper1, rho) -
    genz(upper0, lower1, pnorm_upper0, pnorm_lower1, rho) +
    genz(lower0, lower1, pnorm_lower0, pnorm_lower1, rho);
}

void poly_derivs(double& df_dp, arma::vec& df_dtau1, arma::vec& df_dtau2,
                 double rho, arma::vec tau1, arma::vec tau2,
                 arma::vec pnorm_tau1, arma::vec pnorm_tau2,
                 std::vector<std::vector<int>> n) {

  // This function computes the polychoric loglik derivatives with respect to
  // the latent correlation and thresholds

  // df_dp = Derivative of the latent correlation
  // df_dtau1 = Derivatives of thresholds of the first variable
  // df_dtau2 = Derivatives of thresholds of the second variable
  // rho = latent correlation
  // tau1 = Vector of thresholds of the first variable (It must start at -Infinite and end at Infinite)
  // tau2 = Vector of thresholds of the second variable (It must start at -Infinite and end at Infinite)
  // pnorm_tau1 = pnorm of tau1
  // pnorm_tau2 = pnorm of tau2
  // n = contingency table

  int s = tau1.size()-1L;
  int r = tau2.size()-1L;
  arma::mat pbi(s, r);
  arma::mat dpbi_dp(s, r);
  double denominator = std::sqrt(1-rho*rho);
  for (size_t i = 0; i < s; ++i) { // Loop over the thresholds of first variable
    for (size_t j = 0; j < r; ++j) { // Loop over the thresholds of second variable
      int nij = n[i][j];
      if (nij == 0) continue;
      // CDF of the bivariate normal:
      pbi(i, j) = pbinorm(rho, tau1[i], tau2[j], tau1[i+1], tau2[j+1],
          pnorm_tau1[i], pnorm_tau2[j], pnorm_tau1[i+1], pnorm_tau2[j+1]);
      pbi(i,j) = std::max(pbi(i,j), 1e-16);
      // PDF of the Bivariate normal:
      dpbi_dp(i, j) = dbinorm(rho, tau1[i+1], tau2[j+1]) -
        dbinorm(rho, tau1[i], tau2[j+1]) -
        dbinorm(rho, tau1[i+1], tau2[j]) +
        dbinorm(rho, tau1[i], tau2[j]);
      // Derivative of the latent correlation:
      df_dp -= n[i][j] * dpbi_dp(i,j) / pbi(i,j);
    }
  }

  // Loop over the thresholds of first variable:
  for(int k=0; k < (s-1); ++k) {
    for(int j=0; j < r; ++j) {
      int nkj = n[k][j];
      if (nkj == 0) continue;
      double numerator1 = tau2[j+1]-rho*tau1[k+1];
      double numerator2 = tau2[j]-rho*tau1[k+1];
      double dpbi_dtau1 = Dnorm(tau1[k+1])*(Pnorm(numerator1/denominator) -
                               Pnorm(numerator2/denominator));
      // Derivatives of thresholds for the first variable:
      df_dtau1(k) -= dpbi_dtau1 * (n[k][j]/pbi(k,j)-n[k+1][j]/pbi(k+1,j));
    }
  }

  // Loop over the thresholds of second variable
  for(int m=0; m < (r-1); ++m) {
    for(int i=0; i < s; ++i) {
      int nmi = n[m][i];
      if (nmi == 0) continue;
      double numerator1 = tau1[i+1]-rho*tau2[m+1];
      double numerator2 = tau1[i]-rho*tau2[m+1];
      double dpbi_dtau2 = Dnorm(tau2[m+1])*(Pnorm(numerator1/denominator) -
                               Pnorm(numerator2/denominator));
      // Derivatives of thresholds for the second variable:
      df_dtau2(m) -= dpbi_dtau2 * (n[i][m]/pbi(i,m)-n[i][m+1]/pbi(i,m+1));
    }
  }

}

void poly_dderivs(double& ddf_dp, arma::vec& ddf_dtau1, arma::vec& ddf_dtau2,
                  double dp, arma::vec dtau1, arma::vec dtau2,
                  double rho, arma::vec tau1, arma::vec tau2,
                  arma::vec pnorm_tau1, arma::vec pnorm_tau2,
                  std::vector<std::vector<int>> n) {

  int s = tau1.size()-1L;
  int r = tau2.size()-1L;

  // ---- ADDED: guards + constants ----
  const double C = 1.0 - rho*rho;
  if (!(C > 0.0)) return;                 // avoid NaN when |rho|>=1
  const double eps = 1e-16;               // probability floor
  const double denominator = std::sqrt(C);
  const double ddenominator = (-rho/denominator) * dp; // directional derivative of sqrt(1-rho^2)

  // ---- ADDED: map directions dtau -> bounds (endpoints have zero direction) ----
  auto dt1 = [&](int idx) -> double {     // idx in [0..s], dtau1 corresponds to tau1[1..s-1]
    if (idx <= 0 || idx >= s) return 0.0;
    return dtau1(idx - 1);
  };
  auto dt2 = [&](int idx) -> double {     // idx in [0..r], dtau2 corresponds to tau2[1..r-1]
    if (idx <= 0 || idx >= r) return 0.0;
    return dtau2(idx - 1);
  };

  // ---- ADDED: safe partial derivatives of bivariate pdf wrt x and y (avoid 0*Inf at +/-Inf) ----
  auto dbinorm_dx = [&](double p, double x, double y) -> double {
    if (!std::isfinite(x) || !std::isfinite(y)) return 0.0;
    const double Cp = 1.0 - p*p;
    if (!(Cp > 0.0)) return 0.0;
    const double f = dbinorm(p, x, y);
    if (f == 0.0) return 0.0;
    return -((x - p*y)/Cp) * f;
  };
  auto dbinorm_dy = [&](double p, double x, double y) -> double {
    if (!std::isfinite(x) || !std::isfinite(y)) return 0.0;
    const double Cp = 1.0 - p*p;
    if (!(Cp > 0.0)) return 0.0;
    const double f = dbinorm(p, x, y);
    if (f == 0.0) return 0.0;
    return -((y - p*x)/Cp) * f;
  };

  // ---- ADDED: directional derivative of corner density φ2(rho,x,y) ----
  auto d_phi = [&](double x, double dx, double y, double dy) -> double {
    if (!std::isfinite(x) || !std::isfinite(y)) return 0.0;
    return ddbinorm(rho, x, y) * dp +
      dbinorm_dx(rho, x, y) * dx +
      dbinorm_dy(rho, x, y) * dy;
  };

  // ---- ADDED: directional derivative of u = (a - rho*t)/denominator ----
  auto du_rect = [&](double a, double da, double t, double dt) -> double {
    // u = (a - rho*t) / denominator
    const double num  = a - rho*t;
    const double dnum = da - dp*t - rho*dt;
    return dnum/denominator - (num/(denominator*denominator))*ddenominator;
  };

  arma::mat pbi(s, r);      // CDF of the bivariate normal
  arma::mat dpbi_dp(s, r);  // "PDF inclusion" = ∂pbi/∂rho
  arma::mat ddpbi_dp(s, r); // "d/d rho of ∂pbi/∂rho" inclusion at corners (kept as you had)

  // ---- ADDED: initialize (avoids use of uninitialized entries if you ever skip) ----
  pbi.zeros();
  dpbi_dp.zeros();
  ddpbi_dp.zeros();

  // ---- ADDED: dpbi_dir needs pbi/dpbi_dp computed; we’ll define it now and use after fill ----
  auto dpbi_dir = [&](int i, int j) -> double {

    // rho contribution
    double out = dpbi_dp(i, j) * dp;

    // bounds
    const double aL = tau1[i];
    const double aU = tau1[i+1];
    const double bL = tau2[j];
    const double bU = tau2[j+1];

    const double daL = dt1(i);
    const double daU = dt1(i+1);
    const double dbL = dt2(j);
    const double dbU = dt2(j+1);

    // aU
    if (daU != 0.0 && std::isfinite(aU)) {
      const double u1 = (bU - rho*aU)/denominator;
      const double u2 = (bL - rho*aU)/denominator;
      out += Dnorm(aU) * (Pnorm(u1) - Pnorm(u2)) * daU;
    }
    // aL
    if (daL != 0.0 && std::isfinite(aL)) {
      const double u1 = (bU - rho*aL)/denominator;
      const double u2 = (bL - rho*aL)/denominator;
      out -= Dnorm(aL) * (Pnorm(u1) - Pnorm(u2)) * daL;
    }
    // bU
    if (dbU != 0.0 && std::isfinite(bU)) {
      const double v1 = (aU - rho*bU)/denominator;
      const double v2 = (aL - rho*bU)/denominator;
      out += Dnorm(bU) * (Pnorm(v1) - Pnorm(v2)) * dbU;
    }
    // bL
    if (dbL != 0.0 && std::isfinite(bL)) {
      const double v1 = (aU - rho*bL)/denominator;
      const double v2 = (aL - rho*bL)/denominator;
      out -= Dnorm(bL) * (Pnorm(v1) - Pnorm(v2)) * dbL;
    }

    return out;
  };

  // ===== Your existing fill loop (kept), with ADDED probability floor =====
  for (size_t i = 0; i < (size_t)s; ++i) {
    for (size_t j = 0; j < (size_t)r; ++j) {
      // CDF of the bivariate normal:
      pbi(i, j) = pbinorm(rho, tau1[i], tau2[j], tau1[i+1], tau2[j+1],
          pnorm_tau1[i], pnorm_tau2[j], pnorm_tau1[i+1], pnorm_tau2[j+1]);
      pbi(i, j) = std::max((double)pbi(i, j), eps);

      // PDF inclusion:
      dpbi_dp(i, j) = dbinorm(rho, tau1[i+1], tau2[j+1]) -
        dbinorm(rho, tau1[i],   tau2[j+1]) -
        dbinorm(rho, tau1[i+1], tau2[j]) +
        dbinorm(rho, tau1[i],   tau2[j]);

      // Second derivative inclusion wrt rho (kept, but NOTE: not sufficient alone for full direction)
      ddpbi_dp(i, j) = ddbinorm(rho, tau1[i+1], tau2[j+1]) -
        ddbinorm(rho, tau1[i],   tau2[j+1]) -
        ddbinorm(rho, tau1[i+1], tau2[j]) +
        ddbinorm(rho, tau1[i],   tau2[j]);

      // df_dp -= n[i][j] * dpbi_dp(i,j) / pbi(i,j);

      // ---- FILLED: ddf_dp differential in direction (dp,dtau1,dtau2) ----
        {
          const int nij = n[i][j];
          if (nij != 0) {

            // d(pbi) full directional derivative (rho + moving bounds)
            const double dpbi = dpbi_dir((int)i, (int)j);

            // d(dpbi_dp): dp contribution + moving-corner contributions
            const double aL = tau1[i],     aU = tau1[i+1];
            const double bL = tau2[j],     bU = tau2[j+1];
            const double daL = dt1((int)i),     daU = dt1((int)i+1);
            const double dbL = dt2((int)j),     dbU = dt2((int)j+1);

            const double d_dpbi_dp =
              d_phi(aU, daU, bU, dbU) -
              d_phi(aL, daL, bU, dbU) -
              d_phi(aU, daU, bL, dbL) +
              d_phi(aL, daL, bL, dbL);

            const double pij = (double)pbi(i, j); // already floored
            const double gij = (double)dpbi_dp(i, j);

            // d( gij / pij ) = d_gij/pij - gij * d_pij / pij^2
            ddf_dp -= (double)nij * ( d_dpbi_dp / pij - gij * dpbi / (pij * pij) );
          }
        }
    }
  }

  // Loop over the thresholds of first variable:
  for(int k=0; k < (s-1); ++k) {
    for(int j=0; j < r; ++j) {

      double numerator1 = tau2[j+1]-rho*tau1[k+1];
      double numerator2 = tau2[j]-rho*tau1[k+1];
      double dpbi_dtau1 = Dnorm(tau1[k+1])*(Pnorm(numerator1/denominator) -
                                Pnorm(numerator2/denominator));

      // df_dtau1(k) -= dpbi_dtau1 * (n[k][j]/pbi(k,j)-n[k+1][j]/pbi(k+1,j));

      // ---- FILLED: ddf_dtau1(k) differential in direction (dp,dtau1,dtau2) ----
                                {
                                  const int n0 = n[k][j];
                                  const int n1 = n[k+1][j];
                                  if (n0 != 0 || n1 != 0) {

                                    const double t  = tau1[k+1];
                                    const double dt = dtau1(k);        // direction for tau1[k+1]

                                    // a = φ(t) * (Φ(uU)-Φ(uL))
                                    const double uU = (tau2[j+1] - rho*t)/denominator;
                                    const double uL = (tau2[j]   - rho*t)/denominator;
                                    const double PhiDiff = Pnorm(uU) - Pnorm(uL);
                                    const double a = dpbi_dtau1;

                                    // duU, duL
                                    const double duU = du_rect(tau2[j+1], dt2(j+1), t, dt);
                                    const double duL = du_rect(tau2[j],   dt2(j),   t, dt);

                                    // Avoid 0*Inf in φ(u)*du if u huge
                                    const double phiuU = Dnorm(uU);
                                    const double phiuL = Dnorm(uL);
                                    const double termU = (phiuU == 0.0 ? 0.0 : phiuU * duU);
                                    const double termL = (phiuL == 0.0 ? 0.0 : phiuL * duL);

                                    // da = dφ(t)*PhiDiff + φ(t)*( φ(uU)*duU - φ(uL)*duL )
                                    const double da =
                                      dDnorm(t) * dt * PhiDiff +
                                      Dnorm(t) * (termU - termL);

                                    // B = n0/p1 - n1/p2
                                    const double p1 = std::max((double)pbi(k,   j), eps);
                                    const double p2 = std::max((double)pbi(k+1, j), eps);
                                    const double B  = (double)n0 / p1 - (double)n1 / p2;

                                    // dB = -n0*dp1/p1^2 + n1*dp2/p2^2, where dp1 = d(pbi(k,j)), dp2 = d(pbi(k+1,j))
                                    const double dp1 = dpbi_dir(k,   j);
                                    const double dp2 = dpbi_dir(k+1, j);
                                    const double dB  = -(double)n0 * dp1 / (p1*p1) + (double)n1 * dp2 / (p2*p2);

                                    // d( -a * B ) = -(da*B + a*dB)
                                    ddf_dtau1(k) -= (da * B + a * dB);
                                  }
                                }
    }
  }

  // Loop over the thresholds of second variable
  for(int m=0; m < (r-1); ++m) {
    for(int i=0; i < s; ++i) {

      double numerator1 = tau1[i+1]-rho*tau2[m+1];
      double numerator2 = tau1[i]-rho*tau2[m+1];
      double dpbi_dtau2 = Dnorm(tau2[m+1])*(Pnorm(numerator1/denominator) -
                                Pnorm(numerator2/denominator));

      // df_dtau2(m) -= dpbi_dtau2 * (n[i][m]/pbi(i,m)-n[i][m+1]/pbi(i,m+1));

      // ---- FILLED: ddf_dtau2(m) differential in direction (dp,dtau1,dtau2) ----
                                {
                                  const int n0 = n[i][m];
                                  const int n1 = n[i][m+1];
                                  if (n0 != 0 || n1 != 0) {

                                    const double t  = tau2[m+1];
                                    const double dt = dtau2(m);        // direction for tau2[m+1]

                                    const double uU = (tau1[i+1] - rho*t)/denominator;
                                    const double uL = (tau1[i]   - rho*t)/denominator;
                                    const double PhiDiff = Pnorm(uU) - Pnorm(uL);
                                    const double a = dpbi_dtau2;

                                    const double duU = du_rect(tau1[i+1], dt1(i+1), t, dt);
                                    const double duL = du_rect(tau1[i],   dt1(i),   t, dt);

                                    const double phiuU = Dnorm(uU);
                                    const double phiuL = Dnorm(uL);
                                    const double termU = (phiuU == 0.0 ? 0.0 : phiuU * duU);
                                    const double termL = (phiuL == 0.0 ? 0.0 : phiuL * duL);

                                    const double da =
                                      dDnorm(t) * dt * PhiDiff +
                                      Dnorm(t) * (termU - termL);

                                    const double p1 = std::max((double)pbi(i, m),   eps);
                                    const double p2 = std::max((double)pbi(i, m+1), eps);
                                    const double B  = (double)n0 / p1 - (double)n1 / p2;

                                    const double dp1 = dpbi_dir(i, m);
                                    const double dp2 = dpbi_dir(i, m+1);
                                    const double dB  = -(double)n0 * dp1 / (p1*p1) + (double)n1 * dp2 / (p2*p2);

                                    ddf_dtau2(m) -= (da * B + a * dB);
                                  }
                                }
    }
  }
}

class polycor: public estimators {

public:

  arma::mat R;
  std::vector<arma::vec> taus;
  std::vector<arma::vec> pnorm_tau;
  std::vector<std::vector<std::vector<int>>> n;
  std::vector<arma::uvec> indices_taus;
  arma::uvec indices_R, lower_diag;
  int p;
  double loss;

  void param(arguments_optim& x) {

    for(int i=0; i < p; ++i) {
      taus[i] = x.transparameters(indices_taus[i]);
      taus[i] = arma::join_vert(taus[i], arma::vec({pos_inf}));
      taus[i] = arma::join_vert(arma::vec({neg_inf}), taus[i]);
      pnorm_tau[i] = 0.5 * arma::erfc(-taus[i] * M_SQRT1_2);
    }

    R.elem(lower_diag) = x.transparameters(indices_R);
    R = arma::symmatl(R);
    R.diag().ones(); // Ensure ones in the diagonal of the correlation matrix

  }

  void F(arguments_optim& x) {

    int m = 0L;
    loss = 0.00;
    for(int l=0; l < (p-1L); ++l) { // Loop across first variable
      const size_t s1 = taus[l].size()-1L;
      for(int k=(l+1L); k < p; ++k, ++m) { // Loop across second variable
        // In total, there are m variable pairs
        const size_t s2 = taus[k].size()-1L;
        for (size_t i = 0; i < s1; ++i) { // Loop across first variable thresholds
          for (size_t j = 0; j < s2; ++j) { // Loop across second variable thresholds
            int nmij = n[m][i][j];
            if (nmij == 0) continue;
            loss -= n[m][i][j] * arma::trunc_log(pbinorm(R(l,k),
                                                         taus[l](i), taus[k](j),
                                                         taus[l](i+1), taus[k](j+1),
                                                         pnorm_tau[l](i), pnorm_tau[k](j),
                                                         pnorm_tau[l](i+1), pnorm_tau[k](j+1)));
          }
        }
      }
    }

    x.f += loss;

  }

  void G(arguments_optim& x) {

    std::vector<arma::vec> dfdtaus(p);
    for (size_t i = 0; i < p; ++i) {
      dfdtaus[i].set_size(taus[i].n_elem-2L);
      dfdtaus[i].zeros();
    }
    arma::mat dfdp(p, p, arma::fill::zeros);

    int m = 0L;
    for(int l=0; l < (p-1L); ++l) { // Loop across first variable
      for(int k=(l+1L); k < p; ++k, ++m) { // Loop across second variable

        double dp = 0.0;
        arma::vec dtau1(taus[l].n_elem-2L, arma::fill::zeros);
        arma::vec dtau2(taus[k].n_elem-2L, arma::fill::zeros);

        poly_derivs(dp, dtau1, dtau2,
                    R(l,k), taus[l], taus[k],
                    pnorm_tau[l], pnorm_tau[k], n[m]);

        dfdp(l,k) += dp;
        dfdp(k,l) = dfdp(l,k);
        dfdtaus[l] += dtau1;
        dfdtaus[k] += dtau2;

      }
    }

    for(int i=0; i < p; ++i) {
      x.grad(indices_taus[i]) += dfdtaus[i];
    }

    x.grad(indices_R) += arma::vectorise(dfdp(lower_diag));

  }

  void dG(arguments_optim& x) {

    // Create the directions dtaus and dR:
    std::vector<arma::vec> dtaus(p);
    for(int i=0; i < p; ++i) {
      dtaus[i] = x.dtransparameters(indices_taus[i]);
    }

    arma::mat dR(p, p, arma::fill::zeros);
    dR.elem(lower_diag) = x.dtransparameters(indices_R);
    dR = arma::symmatl(dR);

    // compute the differentials:
    std::vector<arma::vec> ddfdtaus(p);
    for (size_t i = 0; i < p; ++i) {
      ddfdtaus[i].set_size(taus[i].n_elem-2L);
      ddfdtaus[i].zeros();
    }
    arma::mat ddfdp(p, p, arma::fill::zeros);

    int m = 0L;
    for(int l=0; l < (p-1L); ++l) { // Loop across first variable
      for(int k=(l+1L); k < p; ++k, ++m) { // Loop across second variable

        double ddp = 0.0;
        arma::vec ddtau1(taus[l].n_elem-2L, arma::fill::zeros);
        arma::vec ddtau2(taus[k].n_elem-2L, arma::fill::zeros);

        poly_dderivs(ddp, ddtau1, ddtau2,
                     dR(l,k), dtaus[l], dtaus[k],
                     R(l,k), taus[l], taus[k],
                     pnorm_tau[l], pnorm_tau[k], n[m]);

        ddfdp(l,k) += ddp;
        ddfdp(k,l) = ddfdp(l,k);
        ddfdtaus[l] += ddtau1;
        ddfdtaus[k] += ddtau2;

      }
    }

    // store the differentials in the x.dgrad object:
    for(int i=0; i < p; ++i) {
      x.dgrad(indices_taus[i]) += ddfdtaus[i];
    }
    x.dgrad(indices_R) += arma::vectorise(ddfdp(lower_diag));

  }

  void outcomes(arguments_optim& x) {

    doubles.resize(2);
    doubles[0] = loss;
    doubles[1] = -loss;

    matrices.resize(1);
    matrices[0] = R;
    // matrices[1] = loglik;

    list_vectors.resize(2);
    list_vectors[0] = taus;
    list_vectors[1] = pnorm_tau;

  };

};

polycor* choose_polycor(const Rcpp::List& estimator_setup) {

  polycor* myestimator = new polycor();

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  std::vector<arma::uvec> indices_taus = estimator_setup["indices_taus"];
  arma::uvec indices_R = estimator_setup["indices_R"];
  std::vector<std::vector<std::vector<int>>> n = estimator_setup["n"];
  int p = estimator_setup["p"];

  std::vector<arma::vec> taus(p);
  std::vector<arma::vec> pnorm_tau(p);
  arma::mat R(p, p, arma::fill::zeros);
  arma::uvec lower_diag = arma::trimatl_ind(arma::size(R));

  myestimator->indices = indices;
  myestimator->indices_taus = indices_taus;
  myestimator->indices_R = indices_R;
  myestimator->n = n;
  myestimator->p = p;
  myestimator->taus = taus;
  myestimator->pnorm_tau = pnorm_tau;
  myestimator->R = R;
  myestimator->lower_diag = lower_diag;

  return myestimator;

}

// void poly_dderivs(double& ddf_dp, arma::vec& ddf_dtau1, arma::vec& ddf_dtau2,
//                   double dp, arma::vec dtau1, arma::vec dtau2,
//                   double rho, arma::vec tau1, arma::vec tau2,
//                   arma::vec pnorm_tau1, arma::vec pnorm_tau2,
//                   std::vector<std::vector<int>> n) {
//
//   // This function computes the differential of the polychoric loglik derivatives
//
//   // ddf_dp = Directional derivative of the derivative of the latent correlation
//   // ddf_dtau1 = Directional derivative of the derivatives of thresholds of the first variable
//   // ddf_dtau2 = Directional derivative of the derivatives of thresholds of the second variable
//   // rho = latent correlation
//   // dp = Direction for the latent correlation
//   // dtau1 = Directions for the thresholds of the first variable
//   // dtau2 = Directions for the thresholds of the second variable
//   // tau1 = Vector of thresholds of the first variable (It must start at -Infinite and end at Infinite)
//   // tau2 = Vector of thresholds of the second variable (It must start at -Infinite and end at Infinite)
//   // pnorm_tau1 = pnorm of tau1
//   // pnorm_tau2 = pnorm of tau2
//   // n = contingency table
//
//   int s = tau1.size()-1L;
//   int r = tau2.size()-1L;
//   arma::mat pbi(s, r); // CDF of the bivariate normal
//   arma::mat dpbi_dp(s, r); // PDF of the Bivariate normal
//   arma::mat ddpbi_dp(s, r); // Derivative of the PDF of the Bivariate normal
//   double denominator = std::sqrt(1-rho*rho);
//   for (size_t i = 0; i < s; ++i) { // Loop over the thresholds of first variable
//     for (size_t j = 0; j < r; ++j) { // Loop over the thresholds of second variable
//       // CDF of the bivariate normal:
//       pbi(i, j) = pbinorm(rho, tau1[i], tau2[j], tau1[i+1], tau2[j+1],
//           pnorm_tau1[i], pnorm_tau2[j], pnorm_tau1[i+1], pnorm_tau2[j+1]);
//       // PDF of the Bivariate normal:
//       dpbi_dp(i, j) = dbinorm(rho, tau1[i+1], tau2[j+1]) -
//         dbinorm(rho, tau1[i], tau2[j+1]) -
//         dbinorm(rho, tau1[i+1], tau2[j]) +
//         dbinorm(rho, tau1[i], tau2[j]);
//       // Derivative of the PDF of the Bivariate normal:
//       ddpbi_dp(i, j) = ddbinorm(rho, tau1[i+1], tau2[j+1]) -
//         ddbinorm(rho, tau1[i], tau2[j+1]) -
//         ddbinorm(rho, tau1[i+1], tau2[j]) +
//         ddbinorm(rho, tau1[i], tau2[j]);
//       // df_dp -= n[i][j] * dpbi_dp(i,j) / pbi(i,j);
//       // ddf_dp -= help me compute the differential of df_dp in direction dp;
//     }
//   }
//
//   // Loop over the thresholds of first variable:
//   for(int k=0; k < (s-1); ++k) {
//     for(int j=0; j < r; ++j) {
//       double numerator1 = tau2[j+1]-rho*tau1[k+1];
//       double numerator2 = tau2[j]-rho*tau1[k+1];
//       double dpbi_dtau1 = Dnorm(tau1[k+1])*(Pnorm(numerator1/denominator) -
//                                 Pnorm(numerator2/denominator));
//       // df_dtau1(k) -= dpbi_dtau1 * (n[k][j]/pbi(k,j)-n[k+1][j]/pbi(k+1,j));
//       // ddf_dtau1(k) -= help me compute the differential of df_dtau1(k) in direction dtau1(k);
//     }
//   }
//
//   // Loop over the thresholds of second variable
//   for(int m=0; m < (r-1); ++m) {
//     for(int i=0; i < s; ++i) {
//       double numerator1 = tau1[i+1]-rho*tau2[m+1];
//       double numerator2 = tau1[i]-rho*tau2[m+1];
//       double dpbi_dtau2 = Dnorm(tau2[m+1])*(Pnorm(numerator1/denominator) -
//                                 Pnorm(numerator2/denominator));
//       // df_dtau2(m) -= dpbi_dtau2 * (n[i][m]/pbi(i,m)-n[i][m+1]/pbi(i,m+1));
//       // ddf_dtau2(m) -= help me compute the differential of df_dtau2(m) in direction dtau2(m);
//     }
//   }
//
// }
