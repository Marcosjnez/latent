/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 27/08/2025
 */

const double LOG2M_PI05 = 0.5*std::log(2 * M_PI);

// Logarithm gaussian density transformation:

arma::vec ldnorm(arma::vec x, arma::vec mean, arma::vec sd) {
  x -= mean;
  arma::vec sd2 = sd % sd;
  arma::vec lpd = -0.5 * x % x/sd2 - arma::trunc_log(SQRT2M_PI*sd);
  return lpd;
}

class normal:public transformations {

public:

  arma::uvec mu_indices;
  arma::uvec sigma_indices;
  arma::uvec removeNAs;
  arma::vec y, mu, x, x2, sigma, sigma2, sigma3, sigma4;

  void transform() {

    mu_indices = indices_in[1];
    sigma_indices = indices_in[2];
    mu = transparameters.elem(mu_indices);
    sigma = transparameters.elem(sigma_indices);

    x = y - mu;
    sigma2 = sigma % sigma;
    x2 = arma::square(x);
    transparameters = -0.5 * x2/sigma2 -
      arma::trunc_log(sigma) - LOG2M_PI05;

    // Set to zero the missing data:
    removeNAs = indices_in[3];
    transparameters(removeNAs).zeros();

  }

  void update_grad() {

    // const arma::uword n_out = transparameters.n_elem;
    // const arma::uword n_in  = indices_in[0].n_elem;
    // jacob.set_size(n_out, n_in);
    // jacob.zeros();

    sigma3 = sigma2 % sigma;
    arma::vec dmu = grad_in % x/sigma2;
    arma::vec dsigma = grad_in % (x2 - sigma2) / sigma3;

    // Set to zero the missing data:
    dmu.elem(removeNAs).zeros();
    dsigma.elem(removeNAs).zeros();

    grad_out.resize(indices_in[0].n_elem); grad_out.zeros();
    grad_out(mu_indices) += dmu;
    grad_out(sigma_indices) += dsigma;

    // The same but much slower:
    // jacob.set_size(transparameters.n_elem, indices_in[0].n_elem);
    // jacob.zeros();
    // jacob.cols(mu_indices) += arma::diagmat(x/sigma2);
    // jacob.cols(sigma_indices) += arma::diagmat((x2 - sigma2) / sigma3);
    // grad_out = jacob.t() * grad_in;

  }

  void update_hess() {

    const arma::uword n_out = transparameters.n_elem;
    const arma::uword n_in  = indices_in[0].n_elem;

    jacob.set_size(n_out, n_in);
    jacob.zeros();
    sum_djacob.set_size(n_in, n_in);
    sum_djacob.zeros();

    // jacob.cols(mu_indices) += arma::diagmat(x/sigma2);
    // jacob.cols(sigma_indices) += arma::diagmat((x2 - sigma2) / sigma3);
    // Vectorised version:
    arma::vec dmu = x / sigma2;
    arma::vec dsigma = (x2 - sigma2) / sigma3;

    // Set to zero the missing data:
    dmu.elem(removeNAs).zeros();
    dsigma.elem(removeNAs).zeros();

    arma::uvec rows = arma::regspace<arma::uvec>(0, n_out - 1);
    arma::uvec lin_mu    = rows + mu_indices    * n_out;
    arma::uvec lin_sigma = rows + sigma_indices * n_out;
    jacob.elem(lin_mu)    += dmu;
    jacob.elem(lin_sigma) += dsigma;

    // sigma4 = sigma3 % sigma;
    // sum_djacob(mu_indices, mu_indices) += arma::diagmat(-1.0/sigma2 % grad_in);
    // sum_djacob(sigma_indices, sigma_indices) += arma::diagmat((1/sigma2 - 3.0*x2/sigma4) % grad_in);
    // sum_djacob(mu_indices, sigma_indices) += arma::diagmat(-2.0 % x /sigma3 % grad_in);
    // sum_djacob(sigma_indices, mu_indices) = sum_djacob(mu_indices, sigma_indices);
    // Vectorised version:
    sigma4 = sigma3 % sigma;
    // diagonal (mu,mu)
    arma::vec hmu = -1.0 / sigma2;
    arma::vec hsigma = 1.0 / sigma2 - 3.0 * x2 / sigma4;
    arma::vec hmusigma = -2.0 * x / sigma3;
    // Set to zero the missing data:
    hmu.elem(removeNAs).zeros();
    hsigma.elem(removeNAs).zeros();
    hmusigma.elem(removeNAs).zeros();

    lin_mu = mu_indices + mu_indices * n_in;
    sum_djacob.elem(lin_mu) += hmu % grad_in;
    // diagonal (sigma,sigma)
    lin_sigma = sigma_indices + sigma_indices * n_in;
    sum_djacob.elem(lin_sigma) += hsigma % grad_in;
    // off-diagonal (mu,sigma)
    arma::uvec lin_mu_sigma = mu_indices + sigma_indices * n_in;
    arma::vec vals_mu_sigma = hmusigma % grad_in;
    sum_djacob.elem(lin_mu_sigma) += vals_mu_sigma;
    // symmetric counterpart (sigma,mu)
    arma::uvec lin_sigma_mu = sigma_indices + mu_indices * n_in;
    sum_djacob.elem(lin_sigma_mu) = sum_djacob.elem(lin_mu_sigma);

    // hess_in = jacob.t() * hess_out * jacob + sum_djacob;

  }

  void update_vcov() {

  }

  void dconstraints() {

    constraints = false;

  }

  void outcomes() {

    int p = transparameters.n_elem;
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
  arma::vec y = trans_setup["y"];

  mytrans->indices_in = indices_in;
  mytrans->indices_out = indices_out;
  mytrans->y = y;

  return mytrans;

}
