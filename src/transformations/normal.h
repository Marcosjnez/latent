/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 14/02/2026
 */

const double LOG2M_PI05 = 0.5*std::log(2 * M_PI);

// Logarithm gaussian density transformation:

class normal:public transformations {

public:

  arma::uvec indices_mu, indices_sigma, indices_in, indices_out;
  arma::mat y, mu, sigma, sigma2, sigma3, sigma4, dmu, dsigma, jacob;
  int S, J, I, n_in, n_out;

  void transform(arguments_optim& x) {

    arma::cube loglik(S, J, I, arma::fill::zeros);
    mu = arma::reshape(x.transparameters(indices_mu), J, I);
    sigma = arma::reshape(x.transparameters(indices_sigma), J, I);
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

    x.transparameters.elem(indices_out) = arma::vectorise(loglik);

  }

  void update_grad(arguments_optim& x) {

    arma::mat df_dmu(J, I, arma::fill::zeros);
    arma::mat df_dsigma(J, I, arma::fill::zeros);

    arma::vec grad = x.grad(indices_out);

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

    x.grad(indices_mu) += arma::vectorise(df_dmu);
    x.grad(indices_sigma) += arma::vectorise(df_dsigma);

  }

  void dtransform(arguments_optim& x) {

    dmu    = arma::reshape(x.dtransparameters(indices_mu),    J, I);
    dsigma = arma::reshape(x.dtransparameters(indices_sigma), J, I);

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

    x.dtransparameters.elem(indices_out) = arma::vectorise(dloglik);

  }

  void update_dgrad(arguments_optim& x) {

    sigma4 = sigma2 % sigma2;

    arma::mat ddf_dmu(J, I, arma::fill::zeros);
    arma::mat ddf_dsigma(J, I, arma::fill::zeros);

    arma::vec grad_out    = x.grad(indices_out);
    arma::vec dgrad_out   = x.dgrad(indices_out);

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

    x.dgrad(indices_mu)    += arma::vectorise(ddf_dmu);
    x.dgrad(indices_sigma) += arma::vectorise(ddf_dsigma);

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

    indices_in = arma::join_cols(indices_mu, indices_sigma);
    x.vcov(indices_out, indices_out) = jacob * x.vcov(indices_in, indices_in) * jacob.t();

  }

  void dconstraints(arguments_optim& x) {

    constraints = false;

  }

  void outcomes(arguments_optim& x) {

    int p = indices_out.n_elem;
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
  // arma::uvec indices_mu = trans_setup["indices_mu"];
  // arma::uvec indices_sigma = trans_setup["indices_sigma"];
  // arma::uvec indices_out = trans_setup["indices_out"];
  arma::mat y = trans_setup["y"];
  int S = trans_setup["S"];
  int J = trans_setup["J"];
  int I = trans_setup["I"];

  arma::uvec indices_mu = indices_in[0];
  arma::uvec indices_sigma = indices_in[1];
  int n_in  = indices_mu.n_elem + indices_sigma.n_elem;
  int n_out = indices_out[0].n_elem;

  mytrans->indices_mu = indices_mu;
  mytrans->indices_sigma = indices_sigma;
  mytrans->indices_out = indices_out[0];
  mytrans->n_in = n_in;
  mytrans->n_out = n_out;
  mytrans->y = y;
  mytrans->S = S;
  mytrans->J = J;
  mytrans->I = I;

  return mytrans;

}
