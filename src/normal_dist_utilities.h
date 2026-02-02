/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 02/02/2026
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
