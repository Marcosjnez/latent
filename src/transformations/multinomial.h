/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 06/10/2025
 */

// Logarithm multinomial probability transformation:

class multinomial:public transformations {

public:

  arma::vec trans, logtrans, df_dloglik;
  std::vector<arma::mat> peta, eta, df_dpeta, indices;
  std::vector<arma::mat> dpeta, deta;
  arma::uvec K;
  arma::mat y;
  int S, J, I;
  int n_in, n_out;

  void transform(arguments_optim& x) {

    trans = x.transparameters(indices_in[0]);
    logtrans = arma::trunc_log(trans);

    int l=0;
    for(int j=0; j < J; ++j) {
      for(int i=0; i < I; ++i) {
        for(int k=0; k < K[j]; ++k, ++l) {
          peta[j](k,i) = trans(l);
          eta[j](k,i) = logtrans(l);
        }
      }
    }

    arma::cube loglik(S, J, I, arma::fill::zeros);
    for(int i=0; i < I; ++i) {
      for(int j=0; j < J; ++j) {
        for(int s=0; s < S; ++s) {
          if (std::isnan(y(s,j))) continue;
          int value = y(s,j);
          loglik(s,j,i) = eta[j](value,i);
        }
      }
    }

    x.transparameters.elem(indices_out[0]) = arma::vectorise(loglik);

  }

  void update_grad(arguments_optim& x) {

    df_dloglik = x.grad(indices_out[0]);

    for(int j=0; j < J; ++j) {
      df_dpeta[j].zeros();
    }

    int l=0L;
    for(int i=0; i < I; ++i) {
      for(int j=0; j < J; ++j) {
        for(int s=0; s < S; ++s, ++l) {
          if (std::isnan(y(s,j))) continue;
          int value = y(s,j);
          df_dpeta[j](value,i) += df_dloglik(l)/peta[j](value,i);
        }
      }
    }

    arma::vec v = arma::vec();  // initialize empty vector
    for (int j = 0; j < J; ++j) {
      v = arma::join_cols(v, arma::vectorise(df_dpeta[j]));
    }

    x.grad(indices_in[0]) += v;

  }

  void dtransform(arguments_optim& x) {

    arma::vec dtrans_in = x.dtransparameters(indices_in[0]);

    int l=0;
    for(int j=0; j < J; ++j) {
      for(int i=0; i < I; ++i) {
        for(int k=0; k < K[j]; ++k, ++l) {
          dpeta[j](k,i) = dtrans_in(l);
          deta[j](k,i) = dtrans_in(l)/trans(l);
        }
      }
    }

    arma::cube dloglik(S, J, I, arma::fill::zeros);
    for(int i=0; i < I; ++i) {
      for(int j=0; j < J; ++j) {
        for(int s=0; s < S; ++s) {
          if (std::isnan(y(s,j))) continue;
          int value = y(s,j);
          dloglik(s,j,i) = deta[j](value,i);
        }
      }
    }

    x.dtransparameters.elem(indices_out[0]) = arma::vectorise(dloglik);

  }

  void update_dgrad(arguments_optim& x) {

    arma::vec ddf_dloglik = x.dgrad(indices_out[0]);

    std::vector<arma::mat> dgrad_in = df_dpeta;
    for(int j=0; j < J; ++j) {
      dgrad_in[j].zeros();
    }

    int l=0L;
    for(int i=0; i < I; ++i) {
      for(int j=0; j < J; ++j) {
        for(int s=0; s < S; ++s, ++l) {
          if (std::isnan(y(s,j))) continue;
          int value = y(s,j);
          // dgrad_in[j](value,i) += df_dloglik(l)/peta[j](value,i);
          dgrad_in[j](value,i) += (peta[j](value,i)*ddf_dloglik(l) -
            df_dloglik(l)*dpeta[j](value,i))/(peta[j](value,i)*peta[j](value,i));
        }
      }
    }

    arma::vec v = arma::vec();  // initialize empty vector
    for (int j = 0; j < J; ++j) {
      v = arma::join_cols(v, arma::vectorise(dgrad_in[j]));
    }

    x.dgrad(indices_in[0]) += v;

  }

  void jacobian(arguments_optim& x) {

    jacob.set_size(n_out, n_in);
    jacob.zeros();

    int l = 0L;
    for(int i=0; i < I; ++i) {
      for(int j=0; j < J; ++j) {
        for(int s=0; s < S; ++s, ++l) {
          if (std::isnan(y(s,j))) continue;
          int value = y(s,j);
          int index = indices[j](value,i);
          jacob(l, index) = 1/peta[j](value,i);
        }
      }
    }

  }

  void update_vcov(arguments_optim& x) {

  }

  void dconstraints(arguments_optim& x) {

    constraints = false;

  }

  void outcomes(arguments_optim& x) {

    int p = indices_out[0].n_elem;
    arma::vec chisq_p(p, arma::fill::value(1.00));

    vectors.resize(2);
    vectors[0] = dconstr;
    vectors[1] = chisq_p;

    matrices.resize(2);
    matrices[0] = jacob;
    matrices[1] = sum_djacob;

  }

};

multinomial* choose_multinomial(const Rcpp::List& trans_setup) {

  multinomial* mytrans = new multinomial();

  std::vector<arma::uvec> indices_in = trans_setup["indices_in"];
  std::vector<arma::uvec> indices_out = trans_setup["indices_out"];
  arma::mat y = trans_setup["y"];
  arma::uvec K = trans_setup["K"];
  int S = trans_setup["S"];
  int J = trans_setup["J"];
  int I = trans_setup["I"];

  int n_in  = indices_in[0].n_elem;
  int n_out = indices_out[0].n_elem;
  // arma::uvec peta_indices = indices_in[1]; // SxJxIxK indices
  // arma::uvec removeNAs = indices_in[2];
  std::vector<arma::mat> peta(J);
  for(int j=0; j < J; ++j) {
    peta[j].resize(K(j), I);
    peta[j].zeros();
  }
  std::vector<arma::mat> indices = peta;
  int l=0;
  for(int j=0; j < J; ++j) {
    for(int i=0; i < I; ++i) {
      for(int k=0; k < K[j]; ++k, ++l) {
        indices[j](k,i) = l;
      }
    }
  }

  mytrans->indices_in = indices_in;
  mytrans->indices_out = indices_out;
  mytrans->n_in = n_in;
  mytrans->n_out = n_out;
  // mytrans->peta_indices = peta_indices;
  // mytrans->removeNAs = removeNAs;
  mytrans->y = y;
  mytrans->S = S;
  mytrans->J = J;
  mytrans->I = I;
  mytrans->K = K;
  mytrans->peta = peta;
  mytrans->eta = peta;
  mytrans->dpeta = peta;
  mytrans->deta = peta;
  mytrans->df_dpeta = peta;
  mytrans->indices = indices;

  return mytrans;

}
