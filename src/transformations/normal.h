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

  void update_dgrad(arguments_optim& x) {

  }

  void update_hess(arguments_optim& x) {

    // Initialize the jacobian:
    jacob.set_size(n_out, n_in);
    jacob.zeros();

    // Initialize sum_djacob:
    sum_djacob.set_size(n_in, n_in);
    sum_djacob.zeros();

    arma::mat sigma4 = sigma3 % sigma;
    arma::vec grad_out = x.grad(indices_out[0]);
    arma::cube Grad_in(grad_out.memptr(), S, J, I, false);

    int ij_m = 0L;
    int ij_s = J*I;
    int k = 0L;
    for (int i = 0; i < I; ++i) {
      for (int j = 0; j < J; ++j, ++ij_m, ++ij_s) {

        // Accumulate G0, G1, G2 over s (skip NaNs in y)
        double G0 = 0.0, G1 = 0.0, G2 = 0.0;
        for (int s = 0; s < S; ++s, ++k) {
          if (std::isnan(y(s,j))) continue;
          const double x  = y(s,j) - mu(j,i);
          jacob(k, ij_m) = x/sigma2(j,i);
          jacob(k, ij_s) = (x*x - sigma2(j,i))/sigma3(j,i);
          G0 += Grad_in(s,j,i);
          G1 += Grad_in(s,j,i) * x;
          G2 += Grad_in(s,j,i) * x * x;
        }

        // Diagonal blocks
        sum_djacob(ij_m , ij_m) += - G0 / sigma2(j,i);
        sum_djacob(ij_s, ij_s) +=   G0 / sigma2(j,i) - 3.0 * G2 / sigma4(j,i);

        // Off-diagonal (symmetric) block
        const double mixed = - 2.0 * G1 / sigma3(j,i);
        sum_djacob(ij_m, ij_s) += mixed;
        sum_djacob(ij_s, ij_m) += mixed;

      }
    }

    // hess_in = jacob.t() * hess_out * jacob + sum_djacob;

  }

  void update_vcov(arguments_optim& x) {

    // Initialize the jacobian:
    jacob.set_size(n_out, n_in);
    jacob.zeros();

    int ij_m = 0L;
    int ij_s = J*I;
    int k = 0L;
    for (int i = 0; i < I; ++i) {
      for (int j = 0; j < J; ++j, ++ij_m, ++ij_s) {

        double G0 = 0.0, G1 = 0.0, G2 = 0.0;
        for (int s = 0; s < S; ++s, ++k) {
          if (std::isnan(y(s,j))) continue;
          const double x  = y(s,j) - mu(j,i);
          jacob(k, ij_m) = x/sigma2(j,i);
          jacob(k, ij_s) = (x*x - sigma2(j,i))/sigma3(j,i);
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

    mu.resize(J, I); mu.zeros();
    sigma.resize(J, I); sigma.zeros();
    sigma2.resize(J, I); sigma.zeros();

    // for(int i=0; i < I; ++i) {
    //   for(int j=0; j < J; ++j) {
    //     // for(int s=0; s < S; ++s) {
    //     //   if (std::isnan(y(s,j))) continue;
    //     //   mu(j,i) += y(s,j) * w(s,i); // y.t() * w
    //     //   double x = y(s,j) - mu(j,i);
    //     //   sigma2(j,i) += x * x * w(s,i);
    //     // }
    //     mu(j,i) = arma::accu(y.col(j) % w.col(i));
    //     arma::vec xdiff = y.col(j) - mu(j,i);
    //     sigma2(j,i) = arma::accu(xdiff % xdiff % w.col(i));
    //   }
    // }

    for (int i = 0; i < I; ++i) {
      for (int j = 0; j < J; ++j) {
        arma::vec yj = y.col(j);
        arma::vec wi = w.col(i);
        arma::uvec idx = arma::find_finite(yj);
        mu(j,i)     = arma::dot(yj.elem(idx), wi.elem(idx));
        sigma2(j,i) = arma::dot(arma::square(yj.elem(idx) - mu(j,i)), wi.elem(idx));
      }
    }

    sigma = arma::sqrt(sigma2);

    // Block of transformations:

    arma::cube loglik(S, J, I, arma::fill::zeros);
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

    x.transparameters(mu_indices) = arma::vectorise(mu);
    x.transparameters(sigma_indices) = arma::vectorise(sigma);
    x.transparameters(indices_out[0]) = arma::vectorise(loglik);

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
