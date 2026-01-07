/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 27/10/2025
 */

// Logarithm beta density transformation:

class beta:public transformations {

public:

  arma::uvec mu_indices;
  arma::uvec sigma_indices;
  arma::mat y, mu, sigma, sigma2, sigma3, sigma4;
  int S, J, I;
  int n_out;
  int n_in;

  void transform(arguments_optim& x) {

    arma::cube loglik(S, J, I, arma::fill::zeros);
    mu = arma::reshape(x.transparameters(mu_indices), J, I);
    sigma = arma::reshape(x.transparameters(sigma_indices), J, I);
    sigma2 = sigma % sigma;
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

    sigma3 = sigma2 % sigma;

    arma::mat dmu(J, I, arma::fill::zeros);
    arma::mat dsigma(J, I, arma::fill::zeros);

    arma::vec grad = x.grad(indices_out[0]);

    int k=0L;
    for(int i=0; i < I; ++i) {
      for(int j=0; j < J; ++j) {
        for(int s=0; s < S; ++s, ++k) {
          if (std::isnan(y(s,j))) continue;
          double x = (y(s,j)-mu(j,i));
          // arma::uword kk = static_cast<arma::uword>(k);
          // arma::uword idx = indices_out[0](kk);
          dmu(j,i) += grad(k) * x/sigma2(j,i);
          dsigma(j,i) += grad(k) * (x*x - sigma2(j,i))/sigma3(j,i);
        }
      }
    }

    // Fill the gradient:
    x.grad(mu_indices) += arma::vectorise(dmu);
    x.grad(sigma_indices) += arma::vectorise(dsigma);
    // grad_out = jacob.t() * grad_in;

  }

  void dtransform(arguments_optim& x) {

  }

  void update_dgrad(arguments_optim& x) {

  }

  void jacobian(arguments_optim& x) {

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

beta* choose_beta(const Rcpp::List& trans_setup) {

  beta* mytrans = new beta();

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
