/*
 * Author: Marcos Jiménez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 05/03/2023
 *
 */

// Criteria for polychoric correlations

class cor_criterion {

public:

  virtual void F(arguments_cor& x) = 0;

  virtual void gcor(arguments_cor& x) = 0;

  virtual void dgcor(arguments_cor& x) = 0;

};

/*
 * polychoric family
 */

class poly_fam: public cor_criterion {

public:

  void F(arguments_cor& x) {

    x.f = 0;
    int K = 0L;
    for(size_t l=0; l < (x.q-1L); ++l) {
      for(int k=(l+1L); k < x.q; ++k) {
        for (size_t i = 0; i < x.s[l]; ++i) {
          for (size_t j = 0; j < x.s[k]; ++j) {
            // CDF of the bivariate normal:
            double ppi = pbinorm(x.taus[l][i], x.taus[k][j], x.taus[l][i + 1], x.taus[k][j + 1], x.correlation(l, k),
                                 x.mvphi[l][i], x.mvphi[k][j], x.mvphi[l][i+1], x.mvphi[k][j+1]);
            x.f -= x.n[K][i][j] * std::log(ppi) / x.nobs;
          }
        }
        ++K;
      }
    }
    x.f /= x.n_pairs;
  }

  void gcor(arguments_cor& x) {

    x.gcor.zeros();
    int K = 0L;
    for(size_t l=0; l < (x.q-1L); ++l) {
      for(int k=(l+1L); k < x.q; ++k) {
        for (size_t i = 0; i < x.s[l]; ++i) {
          for (size_t j = 0; j < x.s[k]; ++j) {
            // CDF of the bivariate normal:
            double ppi = pbinorm(x.taus[l][i], x.taus[k][j], x.taus[l][i + 1], x.taus[k][j + 1], x.correlation(l, k),
                  x.mvphi[l][i], x.mvphi[k][j], x.mvphi[l][i+1], x.mvphi[k][j+1]);
            // PDF of the Bivariate normal:
            double gij = dbinorm(x.correlation(l, k), x.taus[l][i+1], x.taus[k][j+1]) -
              dbinorm(x.correlation(l, k), x.taus[l][i], x.taus[k][j+1]) -
              dbinorm(x.correlation(l, k), x.taus[l][i+1], x.taus[k][j]) +
              dbinorm(x.correlation(l, k), x.taus[l][i], x.taus[k][j]);
            if(ppi < 1e-09) ppi = 1e-09; // Avoid division by zero
            x.gcor(l, k) -= x.n[K][i][j] / ppi * gij / x.nobs; // Update Gradient
            x.gcor(k, l) = x.gcor(l, k);
          }
        }
        ++K;
      }
    }
    x.gcor /= x.n_pairs;
  }

  void dgcor(arguments_cor& x) {

    x.dgcor.zeros();
    int K = 0L;
    for(size_t l=0; l < (x.q-1L); ++l) {
      for(int k=(l+1L); k < x.q; ++k) {
        for (size_t i = 0; i < x.s[l]; ++i) {
          for (size_t j = 0; j < x.s[k]; ++j) {
            // CDF of the bivariate normal:
            double ppi = pbinorm(x.taus[l][i], x.taus[k][j], x.taus[l][i + 1], x.taus[k][j + 1], x.correlation(l, k),
                                 x.mvphi[l][i], x.mvphi[k][j], x.mvphi[l][i+1], x.mvphi[k][j+1]);
            // PDF of the Bivariate normal:
            double gij = dbinorm(x.correlation(l, k), x.taus[l][i+1], x.taus[k][j+1]) -
              dbinorm(x.correlation(l, k), x.taus[l][i], x.taus[k][j+1]) -
              dbinorm(x.correlation(l, k), x.taus[l][i+1], x.taus[k][j]) +
              dbinorm(x.correlation(l, k), x.taus[l][i], x.taus[k][j]);
            // Derivative of the PDF of the Bivariate normal:
            double hij = ddbinorm(x.correlation(l, k), x.taus[l][i+1], x.taus[k][j+1]) -
              ddbinorm(x.correlation(l, k), x.taus[l][i], x.taus[k][j+1]) -
              ddbinorm(x.correlation(l, k), x.taus[l][i+1], x.taus[k][j]) +
              ddbinorm(x.correlation(l, k), x.taus[l][i], x.taus[k][j]);
            if(ppi < 1e-09) ppi = 1e-09; // Avoid division by zero
            // double term = hij - gij*x.correlation(l, k);
            x.dgcor(l, k) += x.n[K][i][j]*(gij*gij - ppi*hij)/(ppi*ppi) / x.nobs * x.dcor(l, k);
            x.dgcor(k, l) = x.dgcor(l, k);
          }
        }
        ++K;
      }
    }
    x.dgcor /= x.n_pairs;
  }

};

// Choose the cor criteria:

cor_criterion* choose_cor_criterion(std::string cor_fam) {

  cor_criterion *criterion;

  if (cor_fam == "poly") {

    criterion = new poly_fam();

  } else if(cor_fam == "none") {

  } else {

    Rcpp::stop("Available correlations: pearson, poly");

  }

  return criterion;

}
