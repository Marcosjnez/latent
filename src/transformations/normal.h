/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 14/07/2026
 */

const double LOG2M_PI05 = 0.5*std::log(2 * M_PI);

// Logarithm gaussian density transformation:

class normal: public transformations {

public:

  int S, J, I, n_in, n_out;
  arma::uvec indices_mu, indices_sigma, indices_loglik, indices_in, indices_out;
  arma::mat y, mu, sigma2, loglik, sigma4, sigma6, dmu, dsigma2, jacob;

  void transform(arguments_optim& x) {

    mu = arma::reshape(x.transparameters(indices_mu), J, I);
    sigma2 = arma::reshape(x.transparameters(indices_sigma), J, I);
    loglik = arma::reshape(x.transparameters(indices_loglik), S, I);

    sigma2.elem(arma::find(sigma2 < arma::datum::eps)).fill(arma::datum::eps);

    sigma4 = sigma2 % sigma2;
    arma::mat log_sigma2 = arma::trunc_log(sigma2);

    for (int i = 0; i < I; ++i) {
      for (int s = 0; s < S; ++s) {
        double ll = 0.0;
        for (int j = 0; j < J; ++j) {
          if (std::isnan(y(s,j))) continue;
          double x_ = y(s,j) - mu(j,i);
          ll += -0.5 * x_ * x_ / sigma2(j,i) - 0.5 * log_sigma2(j,i) - LOG2M_PI05;
        }
        loglik(s,i) += ll;
      }
    }

    x.transparameters.elem(indices_out) = arma::vectorise(loglik);

  }

  void update_grad(arguments_optim& x) {

    arma::mat df_dmu(J, I, arma::fill::zeros);
    arma::mat df_dsigma2(J, I, arma::fill::zeros);

    arma::vec grad_out = x.grad(indices_out);

    int k = 0L;
    for (int i = 0; i < I; ++i) {
      for (int s = 0; s < S; ++s, ++k) {
        double gk = grad_out(k);
        for (int j = 0; j < J; ++j) {
          if (std::isnan(y(s,j))) continue;
          double x_ = y(s,j) - mu(j,i);
          df_dmu(j,i) += gk * x_ / sigma2(j,i);
          df_dsigma2(j,i) += gk * (x_ * x_ - sigma2(j,i)) / (2.0 * sigma4(j,i));
        }
      }
    }

    x.grad(indices_mu) += arma::vectorise(df_dmu);
    x.grad(indices_sigma) += arma::vectorise(df_dsigma2);

  }

  void dtransform(arguments_optim& x) {

    dmu = arma::reshape(x.dtransparameters(indices_mu), J, I);
    dsigma2 = arma::reshape(x.dtransparameters(indices_sigma), J, I);
    arma::mat dloglik = arma::reshape(x.dtransparameters(indices_loglik), S, I);

    for (int i = 0; i < I; ++i) {
      for (int s = 0; s < S; ++s) {
        double dll = 0.0;
        for (int j = 0; j < J; ++j) {
          if (std::isnan(y(s,j))) continue;
          double x_ = y(s,j) - mu(j,i);
          dll +=
            (x_ / sigma2(j,i)) * dmu(j,i) +
            ((x_ * x_ - sigma2(j,i)) / (2.0 * sigma4(j,i))) * dsigma2(j,i);
        }
        dloglik(s,i) += dll;
      }
    }

    x.dtransparameters.elem(indices_out) = arma::vectorise(dloglik);

  }

  void update_dgrad(arguments_optim& x) {

    sigma6 = sigma4 % sigma2;

    arma::mat ddf_dmu(J, I, arma::fill::zeros);
    arma::mat ddf_dsigma2(J, I, arma::fill::zeros);

    arma::vec grad_out = x.grad(indices_out);
    arma::vec dgrad_out = x.dgrad(indices_out);

    int k = 0L;
    for (int i = 0; i < I; ++i) {
      for (int s = 0; s < S; ++s, ++k) {

        double gk = grad_out(k);
        double dgk = dgrad_out(k);

        for (int j = 0; j < J; ++j) {
          if (std::isnan(y(s,j))) continue;

          double x_ = y(s,j) - mu(j,i);
          double v = sigma2(j,i);
          double v2 = sigma4(j,i);
          double v3 = sigma6(j,i);

          double f_mu = x_ / v;
          double f_sigma2 = (x_ * x_ - v) / (2.0 * v2);

          double df_mu = -dmu(j,i) / v - (x_ / v2) * dsigma2(j,i);

          double df_sigma2 =
            -(x_ / v2) * dmu(j,i) +
            ((0.5 * v - x_ * x_) / v3) * dsigma2(j,i);

          ddf_dmu(j,i) += dgk * f_mu + gk * df_mu;
          ddf_dsigma2(j,i) += dgk * f_sigma2 + gk * df_sigma2;
        }
      }
    }

    x.dgrad(indices_mu) += arma::vectorise(ddf_dmu);
    x.dgrad(indices_sigma) += arma::vectorise(ddf_dsigma2);

  }

  void jacobian(arguments_optim& x) {

    jacob.set_size(n_out, n_in);
    jacob.zeros();

    const int offset_mu = 0;
    const int offset_sigma = J * I;
    const int offset_loglik = 2 * J * I;

    int k = 0L;
    for (int i = 0; i < I; ++i) {
      for (int s = 0; s < S; ++s, ++k) {

        jacob(k, offset_loglik + k) = 1.0;

        for (int j = 0; j < J; ++j) {
          if (std::isnan(y(s,j))) continue;

          const double x_ = y(s,j) - mu(j,i);
          const int ij = j + i * J;

          jacob(k, offset_mu + ij) = x_ / sigma2(j,i);
          jacob(k, offset_sigma + ij) = (x_ * x_ - sigma2(j,i)) / (2.0 * sigma4(j,i));
        }
      }
    }

  }

  void update_vcov(arguments_optim& x) {

    indices_in = arma::join_cols(
      indices_mu,
      arma::join_cols(indices_sigma, indices_loglik)
    );

    // x.vcov(indices_out, indices_out) =
    //   jacob * x.vcov(indices_in, indices_in) * jacob.t();

    arma::mat vcov_in(indices_in.n_elem, indices_in.n_elem);

    for(arma::uword j = 0L; j < indices_in.n_elem; ++j) {
      for(arma::uword i = 0L; i < indices_in.n_elem; ++i) {
        vcov_in(i, j) = x.vcov(indices_in[i], indices_in[j]);
      }
    }

    arma::mat vcov_out = jacob * vcov_in * jacob.t();

    for(arma::uword j = 0L; j < indices_out.n_elem; ++j) {
      for(arma::uword i = 0L; i < indices_out.n_elem; ++i) {
        x.vcov(indices_out[i], indices_out[j]) = vcov_out(i, j);
      }
    }

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

normal* choose_normal(const Rcpp::List& trans_setup) {

  normal* mytrans = new normal();

  std::vector<arma::uvec> indices_in = trans_setup["indices_in"];
  std::vector<arma::uvec> indices_out = trans_setup["indices_out"];
  arma::mat y = trans_setup["y"];
  int S = trans_setup["S"];
  int J = trans_setup["J"];
  int I = trans_setup["I"];

  arma::uvec indices_mu = indices_in[0];
  arma::uvec indices_sigma = indices_in[1];
  arma::uvec indices_loglik = indices_in[2];

  int n_in = indices_mu.n_elem + indices_sigma.n_elem + indices_loglik.n_elem;
  int n_out = indices_out[0].n_elem;

  mytrans->indices_mu = indices_mu;
  mytrans->indices_sigma = indices_sigma;
  mytrans->indices_loglik = indices_loglik;
  mytrans->indices_out = indices_out[0];
  mytrans->n_in = n_in;
  mytrans->n_out = n_out;
  mytrans->y = y;
  mytrans->S = S;
  mytrans->J = J;
  mytrans->I = I;

  return mytrans;

}
