/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 13/11/2025
 */

/*
 * Polychoric correlations (estimating the thresholds and correlations)
 */

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
  arma::mat pbi(s, r); pbi.fill(1e-16);
  arma::mat dpbi_dp(s, r);
  double denominator = std::sqrt(1-rho*rho);
  for (size_t i = 0; i < s; ++i) { // Loop over the thresholds of first variable
    for (size_t j = 0; j < r; ++j) { // Loop over the thresholds of second variable
      double nij = n[i][j];
      if (nij < 1e-16) continue;
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
      df_dp -= nij * dpbi_dp(i,j) / pbi(i,j);
    }
  }

  // Loop over the thresholds of first variable:
  for(int k=0; k < (s-1); ++k) {
    for(int j=0; j < r; ++j) {
      double nkj = n[k][j];
      double nk1j = n[k+1][j];
      if (nkj < 1e-16 && nk1j < 1e-16) continue;
      double numerator1 = tau2[j+1]-rho*tau1[k+1];
      double numerator2 = tau2[j]-rho*tau1[k+1];
      double dpbi_dtau1 = Dnorm(tau1[k+1])*(Pnorm(numerator1/denominator) -
                               Pnorm(numerator2/denominator));
      // Derivatives of thresholds for the first variable:
      df_dtau1(k) -= dpbi_dtau1 * (nkj/pbi(k,j)-nk1j/pbi(k+1,j));
    }
  }

  // Loop over the thresholds of second variable
  for(int m=0; m < (r-1); ++m) {
    for(int i=0; i < s; ++i) {
      double nim = n[i][m];
      double nim1 = n[i][m+1];
      if (nim < 1e-16 && nim1 < 1e-16) continue;
      double numerator1 = tau1[i+1]-rho*tau2[m+1];
      double numerator2 = tau1[i]-rho*tau2[m+1];
      double dpbi_dtau2 = Dnorm(tau2[m+1])*(Pnorm(numerator1/denominator) -
                               Pnorm(numerator2/denominator));
      // Derivatives of thresholds for the second variable:
      df_dtau2(m) -= dpbi_dtau2 * (nim/pbi(i,m)-nim1/pbi(i,m+1));
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
  double loss, N;

  void param(arguments_optim& x) {

    for(int i=0; i < p; ++i) {
      taus[i] = x.transparameters(indices_taus[i]);
      taus[i] = arma::join_vert(taus[i], arma::vec({pos_inf}));
      taus[i] = arma::join_vert(arma::vec({neg_inf}), taus[i]);
      pnorm_tau[i] = 0.5 * arma::erfc(-taus[i] * M_SQRT1_2);
    }

    R.elem(lower_diag) = x.transparameters(indices_R);
    R = arma::symmatl(R);
    // R.diag().ones(); // Ensure ones in the diagonal of the correlation matrix

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
            if (nmij < 1e-16) continue;
            loss -= n[m][i][j] * arma::trunc_log(pbinorm(R(l,k),
                                                         taus[l](i), taus[k](j),
                                                         taus[l](i+1), taus[k](j+1),
                                                         pnorm_tau[l](i), pnorm_tau[k](j),
                                                         pnorm_tau[l](i+1), pnorm_tau[k](j+1)));
          }
        }
      }
    }

    x.f += loss/N; // loglik divided by N

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
        // Rprintf("dfdtaus[l] = %.10f\n", dtau1);
        // Rprintf("dfdtaus[k] = %.10f\n", dtau2);

      }
    }

    for(int i=0; i < p; ++i) {
      x.grad(indices_taus[i]) += dfdtaus[i]/N;
    }

    x.grad(indices_R) += arma::vectorise(dfdp(lower_diag))/N;

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
      x.dgrad(indices_taus[i]) += ddfdtaus[i]/N;
    }
    x.dgrad(indices_R) += arma::vectorise(ddfdp(lower_diag))/N;

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
  double N = estimator_setup["N"];

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
  myestimator->N = N;

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
