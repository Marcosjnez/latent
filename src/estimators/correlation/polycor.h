/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 17/07/2025
 */

/*
 * Polychoric correlation logarithm likelihood
 */

class polycor: public estimators {

public:

  arma::mat R;
  std::vector<arma::vec> taus; // CAMBIAR A STD::VECTOR<DOUBLE>
  std::vector<arma::vec> mvphi;
  std::vector<std::vector<std::vector<int>>> n;
  int p;
  double f0;

  void param() {

    for(int i=0; i < p; ++i) {
      arma::uvec taus_indices = indices[i+1L];
      taus[i] = transparameters(taus_indices);
      taus[i]  = arma::join_vert(taus[i], arma::vec({pos_inf}));
      taus[i]  = arma::join_vert(arma::vec({neg_inf}), taus[i]);
      mvphi[i] = 0.5 * arma::erfc(-taus[i] * M_SQRT1_2);
    }
    R = arma::reshape(transparameters(indices[p+1L]), p, p);

    int m = 0L;

    f0 = 0.0;
    for(size_t l=0; l < (p-1L); ++l) {
      const size_t s1 = taus[l].size()-1L;
      for(int k=(l+1L); k < p; ++k) {
        const size_t s2 = taus[k].size()-1L;
        for (size_t i = 0; i < s1; ++i) {
          for (size_t j = 0; j < s2; ++j) {
            f0 -= n[m][i][j] * arma::trunc_log(pbinorm(taus[l][i], taus[k][j],
                                               taus[l][i+1], taus[k][j+1],
                                               R(l,k),
                                               mvphi[l][i], mvphi[k][j],
                                               mvphi[l][i+1], mvphi[k][j+1]));
          }
        }
        ++m;
      }
    }

  }

  void F() {

    f = f0;

  }

  void G() {

    grad.set_size(transparameters.n_elem);
    grad.zeros();

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
                    R(l,k), taus[l], taus[k], mvphi[l], mvphi[k], n[h]);

        dfdp(l,k) += 0.5*dp;
        dfdp(k,l) = dfdp(l,k);
        dfdtaus[l] += dtau1;
        dfdtaus[k] += dtau2;
        ++h;
      }
    }

    for(int i=0; i < p; ++i) {
      arma::uvec taus_indices = indices[i+1L];
      grad(taus_indices) += dfdtaus[i];
    }

    grad(indices[p+1L]) += arma::vectorise(dfdp);

    // arma::vec v1;
    // for (size_t i = 0; i < p; ++i) {
    //   v1 = arma::join_vert(v1, dfdtaus[i]);
    // }
    // // Then, vectorize the matrix:
    // arma::vec v2 = arma::vectorise(dfdp);
    // // Concatenate both:
    // grad += arma::join_vert(v1, v2);

  }

  void dG() {}

  void E() {}

  void M() {}

  void H() {}

  void outcomes() {

    doubles.resize(1);

    vectors.resize(1);

    matrices.resize(1);
    matrices[0] = R;
    // matrices[1] = loglik;

    cubes.resize(1);

    list_vectors.resize(2);
    list_vectors[0] = taus;
    list_vectors[1] = mvphi;

    list_matrices.resize(1);

  };

};

polycor* choose_polycor(const Rcpp::List& estimator_setup) {

  polycor* myestimator = new polycor();

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  // std::vector<arma::vec> mvphi = estimator_setup["mvphi"];
  std::vector<std::vector<std::vector<int>>> n = estimator_setup["n"];
  int p = estimator_setup["p"];

  std::vector<arma::vec> taus(p);
  std::vector<arma::vec> mvphi(p);

  myestimator->indices = indices;
  myestimator->n = n;
  myestimator->p = p;
  myestimator->taus = taus;
  myestimator->mvphi = mvphi;

  return myestimator;

}
