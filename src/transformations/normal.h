/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 06/10/2025
 */

const double LOG2M_PI05 = 0.5*std::log(2 * M_PI);

// Logarithm gaussian density transformation:

class normal:public transformations {

public:

  arma::uvec mu_indices;
  arma::uvec sigma_indices;
  arma::mat y, mu, sigma, sigma2, sigma3, sigma4;
  arma::mat dmu, dsigma;
  int S, J, I;
  int n_out;
  int n_in;

  void transform(arguments_optim& x) {

    arma::cube loglik(S, J, I, arma::fill::zeros);
    mu = arma::reshape(x.transparameters(mu_indices), J, I);
    sigma = arma::reshape(x.transparameters(sigma_indices), J, I);
    sigma2 = sigma % sigma;
    sigma3 = sigma2 % sigma;
    arma::mat log_sigma = arma::trunc_log(sigma);

    for(int i=0; i < I; ++i) {
      for(int j=0; j < J; ++j) {
        for(int s=0; s < S; ++s) {
          if (std::isnan(y(s,j))) continue;
          double x = y(s,j)-mu(j,i);
          loglik(s,j,i) = -0.5 * x * x / sigma2(j,i) - log_sigma(j,i) - LOG2M_PI05;
        }
      }
    }

    x.transparameters.elem(indices_out[0]) = arma::vectorise(loglik);

  }

  void update_grad(arguments_optim& x) {

    arma::mat df_dmu(J, I, arma::fill::zeros);
    arma::mat df_dsigma(J, I, arma::fill::zeros);

    arma::vec grad = x.grad(indices_out[0]);

    int k=0L;
    for(int i=0; i < I; ++i) {
      for(int j=0; j < J; ++j) {
        for(int s=0; s < S; ++s, ++k) {
          if (std::isnan(y(s,j))) continue;
          double x = (y(s,j)-mu(j,i));
          df_dmu(j,i) += grad(k) * x/sigma2(j,i);
          df_dsigma(j,i) += grad(k) * (x*x - sigma2(j,i))/sigma3(j,i);
        }
      }
    }

    x.grad(mu_indices) += arma::vectorise(df_dmu);
    x.grad(sigma_indices) += arma::vectorise(df_dsigma);

  }

  void dtransform(arguments_optim& x) {

    dmu    = arma::reshape(x.dtransparameters(mu_indices),    J, I);
    dsigma = arma::reshape(x.dtransparameters(sigma_indices), J, I);

    arma::cube dloglik(S, J, I, arma::fill::zeros);
    for (int i = 0; i < I; ++i) {
      for (int j = 0; j < J; ++j) {
        const double dmu_ji    = dmu(j,i);
        const double dsigma_ji = dsigma(j,i);
        for (int s = 0; s < S; ++s) {
          if (std::isnan(y(s,j))) continue;
          double x_ = y(s,j) - mu(j,i);
          dloglik(s,j,i) =
            (x_ / sigma2(j,i)) * dmu_ji +
            ((x_*x_ - sigma2(j,i)) / sigma3(j,i)) * dsigma_ji;
        }
      }
    }

    x.dtransparameters.elem(indices_out[0]) = arma::vectorise(dloglik);

  }

  void update_dgrad(arguments_optim& x) {

    // arma::mat ddf_dmu(J, I, arma::fill::zeros);
    // arma::mat ddf_dsigma(J, I, arma::fill::zeros);
    //
    // arma::vec dgrad_out = x.dgrad(indices_out[0]);
    //
    // int k=0L;
    // for(int i=0; i < I; ++i) {
    //   for(int j=0; j < J; ++j) {
    //     for(int s=0; s < S; ++s, ++k) {
    //       if (std::isnan(y(s,j))) continue;
    //       // ddf_dmu(j,i) +=
    //       // ddf_dsigma(j,i) +=
    //     }
    //   }
    // }
    //
    // x.dgrad(mu_indices) += arma::vectorise(ddf_dmu);
    // x.dgrad(sigma_indices) += arma::vectorise(ddf_dsigma);

    sigma4 = sigma2 % sigma2;

    arma::mat ddf_dmu(J, I, arma::fill::zeros);
    arma::mat ddf_dsigma(J, I, arma::fill::zeros);

    arma::vec grad_out    = x.grad(indices_out[0]);
    arma::vec dgrad_out   = x.dgrad(indices_out[0]);

    int k = 0L;
    for (int i = 0; i < I; ++i) {
      for (int j = 0; j < J; ++j) {

        double dmu_ji    = dmu(j,i);
        double dsigma_ji = dsigma(j,i);

        for (int s = 0; s < S; ++s, ++k) {
          if (std::isnan(y(s,j))) continue;

          double x_   = y(s,j) - mu(j,i);
          double s2   = sigma2(j,i);
          double s3   = sigma3(j,i);
          double s4   = sigma4(j,i);

          double f_mu     = x_ / s2;
          double f_sigma  = (x_*x_ - s2) / s3;

          double df_mu    = -dmu_ji / s2 - (2.0 * x_ / s3) * dsigma_ji;
          double df_sigma = -(2.0 * x_ / s3) * dmu_ji +
            ((s2 - 3.0 * x_*x_) / s4) * dsigma_ji;

          double gk   = grad_out(k);
          double dgk  = dgrad_out(k);

          ddf_dmu(j,i)    += dgk * f_mu    + gk * df_mu;
          ddf_dsigma(j,i) += dgk * f_sigma + gk * df_sigma;
        }
      }
    }

    x.dgrad(mu_indices)    += arma::vectorise(ddf_dmu);
    x.dgrad(sigma_indices) += arma::vectorise(ddf_dsigma);

  }

  void jacobian(arguments_optim& x) {

    // Initialize the jacobian:
    jacob.set_size(n_out, n_in);
    jacob.zeros();

    int ij_m = 0L;
    int ij_s = J*I;
    int k = 0L;
    for (int i = 0; i < I; ++i) {
      for (int j = 0; j < J; ++j, ++ij_m, ++ij_s) {

        for (int s = 0; s < S; ++s, ++k) {
          if (std::isnan(y(s,j))) continue;
          const double x  = y(s,j) - mu(j,i);
          jacob(k, ij_m) = x/sigma2(j,i);
          jacob(k, ij_s) = (x*x - sigma2(j,i))/sigma3(j,i);
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
    arma::vec chisq_p(p, arma::fill::value(2.00));

    vectors.resize(2);
    vectors[0] = dconstr;
    vectors[1] = chisq_p;

    matrices.resize(2);
    matrices[0] = jacob;
    matrices[1] = sum_djacob;

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

  arma::uvec mu_indices = indices_in[1];
  arma::uvec sigma_indices = indices_in[2];
  // arma::uvec item_indices = indices_in[3];
  int n_in  = indices_in[0].n_elem;
  int n_out = indices_out[0].n_elem;

  mytrans->indices_in = indices_in;
  mytrans->indices_out = indices_out;
  mytrans->n_in = n_in;
  mytrans->n_out = n_out;

  mytrans->mu_indices = mu_indices;
  mytrans->sigma_indices = sigma_indices;
  // mytrans->item_indices = item_indices;
  mytrans->y = y;
  mytrans->S = S;
  mytrans->J = J;
  mytrans->I = I;

  return mytrans;

}
