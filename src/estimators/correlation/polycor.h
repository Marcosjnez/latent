/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 13/11/2025
 */

/*
 * Polychoric correlations (estimating the taus and correlations)
 */

class polycor: public estimators {

public:

  arma::mat R;
  std::vector<arma::vec> taus; // CAMBIAR A STD::VECTOR<DOUBLE>
  std::vector<arma::vec> mvphi;
  std::vector<std::vector<std::vector<int>>> n;
  std::vector<arma::uvec> indices_taus;
  arma::uvec indices_R, lower_diag;
  int p;
  double loss = 0.00, loglik = 0.00;

  void param(arguments_optim& x) {

    for(int i=0; i < p; ++i) {
      taus[i] = x.transparameters(indices_taus[i]);
      taus[i] = arma::join_vert(taus[i], arma::vec({pos_inf}));
      taus[i] = arma::join_vert(arma::vec({neg_inf}), taus[i]);
      mvphi[i] = 0.5 * arma::erfc(-taus[i] * M_SQRT1_2);
    }

    R.elem(lower_diag) = x.transparameters(indices_R);
    R = arma::symmatl(R);

  }

  void F(arguments_optim& x) {

    int m = 0L;
    loss = 0.00;
    for(size_t l=0; l < (p-1L); ++l) {
      const size_t s1 = taus[l].size()-1L;
      for(int k=(l+1L); k < p; ++k) {
        const size_t s2 = taus[k].size()-1L;
        for (size_t i = 0; i < s1; ++i) {
          for (size_t j = 0; j < s2; ++j) {
            loss -= n[m][i][j] * arma::trunc_log(pbinorm(taus[l](i), taus[k](j),
                                                         taus[l](i+1), taus[k](j+1),
                                                         R(l,k),
                                                         mvphi[l](i), mvphi[k](j),
                                                         mvphi[l](i+1), mvphi[k](j+1)));
          }
        }
        ++m;
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

    int h = 0L;
    for(size_t l=0; l < (p-1L); ++l) {
      for(int k=(l+1L); k < p; ++k) {

        double dp = 0.0;
        arma::vec dtau1(taus[l].n_elem-2L, arma::fill::zeros);
        arma::vec dtau2(taus[k].n_elem-2L, arma::fill::zeros);

        poly_derivs(dp, dtau1, dtau2,
                    R(l,k), taus[l], taus[k],
                    mvphi[l], mvphi[k], n[h]);

        dfdp(l,k) += 0.5*dp;
        dfdp(k,l) = dfdp(l,k);
        dfdtaus[l] += dtau1;
        dfdtaus[k] += dtau2;
        ++h;

      }
    }

    for(int i=0; i < p; ++i) {
      x.grad(indices_taus[i]) += dfdtaus[i];
    }

    dfdp *= 2;
    dfdp.diag() *= 0.5;
    x.grad(indices_R) += arma::vectorise(dfdp(lower_diag));

  }

  void dG(arguments_optim& x) {

    x.dgrad.elem(indices_R) += x.dtransparameters(indices[0]);

  }

  void outcomes(arguments_optim& x) {

    doubles.resize(2);
    loglik = -loss;
    doubles[0] = loss;
    doubles[1] = loglik;

    matrices.resize(1);
    matrices[0] = R;
    // matrices[1] = loglik;

    list_vectors.resize(2);
    list_vectors[0] = taus;
    list_vectors[1] = mvphi;

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
  std::vector<arma::vec> mvphi(p);
  arma::mat R(p, p, arma::fill::zeros);
  arma::uvec lower_diag = arma::trimatl_ind(arma::size(R));

  myestimator->indices = indices;
  myestimator->indices_taus = indices_taus;
  myestimator->indices_R = indices_R;
  myestimator->n = n;
  myestimator->p = p;
  myestimator->taus = taus;
  myestimator->mvphi = mvphi;
  myestimator->R = R;
  myestimator->lower_diag = lower_diag;

  return myestimator;

}
