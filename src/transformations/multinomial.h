/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 06/10/2025
 */

// Logarithm multinomial probability transformation:

class multinomial:public transformations {

public:

  arma::vec trans, logtrans;
  std::vector<arma::mat> peta, eta, dpeta, indices;
  arma::uvec K;
  arma::vec dloglik;
  arma::mat y;
  int S, J, I;
  int n_in, n_out;

  void transform(arguments_optim& x) {

    trans = x.transparameters(indices_in[0]);
    logtrans = arma::trunc_log(trans);

    arma::cube loglik(S, J, I, arma::fill::zeros);
    int l=0;
    for(int j=0; j < J; ++j) {
      for(int i=0; i < I; ++i) {
        for(int k=0; k < K[j]; ++k, ++l) {
          peta[j](k,i) = trans(l);
          eta[j](k,i) = logtrans(l);
        }
      }
    }

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

    // jacob.set_size(transparameters.n_elem, parameters.n_elem);
    // jacob.zeros();

    dloglik = x.grad(indices_out[0]);

    for(int j=0; j < J; ++j) {
      dpeta[j].zeros();
    }

    int l=0L;
    for(int i=0; i < I; ++i) {
      for(int j=0; j < J; ++j) {
        for(int s=0; s < S; ++s, ++l) {
          if (std::isnan(y(s,j))) continue;
          int value = y(s,j);
          dpeta[j](value,i) += dloglik(l)/peta[j](value,i);
        }
      }
    }

    arma::vec v = arma::vec();  // initialize empty vector
    for (int j = 0; j < J; ++j) {
      v = arma::join_cols(v, arma::vectorise(dpeta[j]));
    }

    x.grad(indices_in[0]) += v;

  }

  void update_dparam(arguments_optim& x) {

  }

  void update_dgrad(arguments_optim& x) {

  }

  void update_hess(arguments_optim& x) {

    jacob.set_size(n_out, n_in);
    jacob.zeros();
    sum_djacob.set_size(n_in, n_in);
    sum_djacob.zeros();

    int l = 0L;
    for(int i=0; i < I; ++i) {
      for(int j=0; j < J; ++j) {
        for(int s=0; s < S; ++s, ++l) {
          if (std::isnan(y(s,j))) continue;
          int value = y(s,j);
          int index = indices[j](value,i);
          jacob(l, index) = 1/peta[j](value,i);
          sum_djacob(index, index) -= dloglik(l)/(peta[j](value,i)*peta[j](value,i));
        }
      }
    }

    // hess_in = jacob.t() * hess_out * jacob + sum_djacob;

  }

  void update_vcov(arguments_optim& x) {

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

  void dconstraints(arguments_optim& x) {

    constraints = false;

  }

  void M(arguments_optim& x) {

    // New number of subjects in each class:
    arma::vec freqs_i = arma::sum(x.freqs, 0).t();
    // Standardize the frequencies in each class:
    arma::mat w = x.freqs;
    w.each_row() /= freqs_i.t(); // Columns sum up to one

    for(int j=0; j < J; ++j) {
      peta[j].zeros();
      for(int i=0; i < I; ++i) {
        for(int s=0; s < S; ++s) {
          if (std::isnan(y(s,j))) continue;
          int value = y(s,j);
          peta[j](value,i) += w(s,i);
        }
        eta[j].col(i) = arma::trunc_log(peta[j].col(i));
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

    for(int j=0; j < J; ++j) {
      for(int i=0; i < I; ++i) {
        eta[j].col(i) -= eta[j](0,i);
      }
    }

    arma::vec allpeta = arma::vec();
    arma::vec alleta = arma::vec();
    for (int j = 0; j < J; ++j) {
      allpeta = arma::join_cols(allpeta, arma::vectorise(peta[j]));
      alleta = arma::join_cols(alleta, arma::vectorise(eta[j]));
    }

    x.transparameters.elem(indices_in[1]) = alleta;
    x.transparameters.elem(indices_in[0]) = allpeta;
    x.transparameters.elem(indices_out[0]) = arma::vectorise(loglik);

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
  std::vector<arma::mat> eta = peta;
  std::vector<arma::mat> dpeta = peta;
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
  mytrans->eta = eta;
  mytrans->dpeta = dpeta;
  mytrans->indices = indices;

  return mytrans;

}
