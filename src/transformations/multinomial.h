/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 14/07/2026
 */

// Logarithm multinomial probability transformation:

class multinomial: public transformations {

public:

  int S, J, I, n_in, n_out;
  arma::uvec indices_peta, indices_loglik, indices_in, indices_out, K;
  arma::vec trans, logtrans;
  std::vector<arma::mat> peta, eta, df_dpeta, indices;
  std::vector<arma::mat> dpeta, deta;
  arma::mat y, loglik, jacob;

  void transform(arguments_optim& x) {

    trans = x.transparameters.elem(indices_peta);
    logtrans = arma::trunc_log(trans);
    loglik = arma::reshape(x.transparameters.elem(indices_loglik), S, I);

    int l = 0;
    for (int j = 0; j < J; ++j) {
      for (int i = 0; i < I; ++i) {
        for (int k = 0; k < K[j]; ++k, ++l) {
          peta[j](k,i) = trans(l);
          eta[j](k,i) = logtrans(l);
        }
      }
    }

    for (int i = 0; i < I; ++i) {
      for (int s = 0; s < S; ++s) {
        double ll = 0.0;
        for (int j = 0; j < J; ++j) {
          if (std::isnan(y(s,j))) continue;
          int value = y(s,j);
          ll += eta[j](value,i);
        }
        loglik(s,i) += ll;
      }
    }

    x.transparameters.elem(indices_out) = arma::vectorise(loglik);

  }

  void update_grad(arguments_optim& x) {

    // arma::vec grad_out = x.grad_init(indices_out);
    arma::vec grad_out = x.grad(indices_out);

    for (int j = 0; j < J; ++j) {
      df_dpeta[j].zeros();
    }

    int l = 0L;
    for (int i = 0; i < I; ++i) {
      for (int s = 0; s < S; ++s, ++l) {
        double gk = grad_out(l);
        for (int j = 0; j < J; ++j) {
          if (std::isnan(y(s,j))) continue;
          int value = y(s,j);
          df_dpeta[j](value,i) += gk / peta[j](value,i);
        }
      }
    }

    // arma::mat df_dloglik(S, I, arma::fill::zeros);
    arma::vec v;
    for (int j = 0; j < J; ++j) {
      v = arma::join_cols(v, arma::vectorise(df_dpeta[j]));
      // arma::vec(df_dloglik.memptr(), df_dloglik.n_elem, false, true) += grad_out;
      // df_dloglik += grad_out;
    }

    x.grad.elem(indices_peta) += v;
    // x.grad.elem(indices_loglik) += grad_out;
    // x.grad.elem(indices_loglik) += arma::vectorise(df_dloglik);

  }

  void dtransform(arguments_optim& x) {

    arma::vec dtrans_in = x.dtransparameters.elem(indices_peta);

    int l = 0;
    for (int j = 0; j < J; ++j) {
      for (int i = 0; i < I; ++i) {
        for (int k = 0; k < K[j]; ++k, ++l) {
          dpeta[j](k,i) = dtrans_in(l);
          deta[j](k,i) = dtrans_in(l) / trans(l);
        }
      }
    }

    arma::mat dloglik = arma::reshape(x.dtransparameters.elem(indices_loglik), S, I);

    for (int i = 0; i < I; ++i) {
      for (int s = 0; s < S; ++s) {
        double dll = 0.0;
        for (int j = 0; j < J; ++j) {
          if (std::isnan(y(s,j))) continue;
          int value = y(s,j);
          dll += deta[j](value,i);
        }
        dloglik(s,i) += dll;
      }
    }

    x.dtransparameters.elem(indices_out) = arma::vectorise(dloglik);

  }

  void update_dgrad(arguments_optim& x) {

    arma::vec grad_out = x.grad.elem(indices_out);
    arma::vec dgrad_out = x.dgrad.elem(indices_out);

    std::vector<arma::mat> dgrad_in = df_dpeta;
    for (int j = 0; j < J; ++j) {
      dgrad_in[j].zeros();
    }

    int l = 0L;
    for (int i = 0; i < I; ++i) {
      for (int s = 0; s < S; ++s, ++l) {
        double gk = grad_out(l);
        double dgk = dgrad_out(l);

        for (int j = 0; j < J; ++j) {
          if (std::isnan(y(s,j))) continue;
          int value = y(s,j);
          double p = peta[j](value,i);
          double dp = dpeta[j](value,i);

          dgrad_in[j](value,i) += (p * dgk - gk * dp) / (p * p);
        }
      }
    }

    arma::vec v;
    for (int j = 0; j < J; ++j) {
      v = arma::join_cols(v, arma::vectorise(dgrad_in[j]));
    }

    x.dgrad.elem(indices_peta) += v;
    // x.dgrad.elem(indices_loglik) += dgrad_out;

  }

  void jacobian(arguments_optim& x) {

    jacob.set_size(n_out, n_in);
    jacob.zeros();

    const int offset_peta = 0;
    const int offset_loglik = indices_peta.n_elem;

    int l = 0L;
    for (int i = 0; i < I; ++i) {
      for (int s = 0; s < S; ++s, ++l) {

        jacob(l, offset_loglik + l) = 1.0;

        for (int j = 0; j < J; ++j) {
          if (std::isnan(y(s,j))) continue;
          int value = y(s,j);
          int index = indices[j](value,i);
          jacob(l, offset_peta + index) = 1.0 / peta[j](value,i);
        }
      }
    }

  }

  void update_vcov(arguments_optim& x) {

    indices_in = arma::join_cols(indices_peta, indices_loglik);
    x.vcov(indices_out, indices_out) =
      jacob * x.vcov(indices_in, indices_in) * jacob.t();

  }

  void dconstraints(arguments_optim& x) {

    constraints = false;

  }

  void outcomes(arguments_optim& x) {

    matrices.resize(1);
    matrices[0] = jacob;
    names_matrices.resize(1);
    names_matrices[0] = "jacobian";

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

  arma::uvec indices_peta = indices_in[0];
  arma::uvec indices_loglik = indices_in[1];

  int n_in = indices_peta.n_elem + indices_loglik.n_elem;
  int n_out = indices_out[0].n_elem;

  std::vector<arma::mat> peta(J);
  for (int j = 0; j < J; ++j) {
    peta[j].resize(K(j), I);
    peta[j].zeros();
  }

  std::vector<arma::mat> indices = peta;

  int l = 0;
  for (int j = 0; j < J; ++j) {
    for (int i = 0; i < I; ++i) {
      for (int k = 0; k < K[j]; ++k, ++l) {
        indices[j](k,i) = l;
      }
    }
  }

  mytrans->indices_peta = indices_peta;
  mytrans->indices_loglik = indices_loglik;
  mytrans->indices_out = indices_out[0];
  mytrans->n_in = n_in;
  mytrans->n_out = n_out;
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
