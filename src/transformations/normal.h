/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 17/08/2025
 */

const double LOG2M_PI05 = 0.5*std::log(2 * M_PI);

// Logarithm gaussian density transformation:

arma::vec ldnorm(arma::vec x, arma::vec mean, arma::vec sd) {
  x -= mean;
  arma::vec sd2 = sd % sd;
  arma::vec lpd = -0.5 * x % x/sd2 - arma::trunc_log(SQRT2M_PI*sd);
  // lpd.elem( arma::find_nonfinite(lpd) ).zeros();
  return lpd;
}

// [[Rcpp::export]]
arma::mat reshape_block(const arma::vec& v, std::size_t k) {
  // Reshape a vector as a block matrix with k columns
  std::size_t n = v.n_elem;
  std::size_t bs = n / k;
  arma::mat M(n, k, arma::fill::zeros);
  for (std::size_t t = 0; t < n; ++t) {
    M(t, t / bs) = v(t);
  }
  return M;
}

class normal:public transformations {

public:

  arma::uvec mu_indices;
  arma::uvec sigma_indices;
  arma::vec y, mu, x, x2, sigma, sigma2, sigma3, sigma4;
  // arma::vec newgrad;

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

  }

  void jacobian() {

    // jacob.set_size(transparameters.n_elem, parameters.n_elem);
    // jacob.zeros();

    sigma3 = sigma2 % sigma;
    arma::vec dmu = grad % x/sigma2;
    arma::vec dsigma = grad % (x2 - sigma2) / sigma3;

    grad.resize(indices_in[0].n_elem); grad.zeros();
    grad(mu_indices) += dmu;
    grad(sigma_indices) += dsigma;

  }

  void d2jacobian() {

    // jacob2.set_size(y.n_elem, parameters.n_elem, parameters.n_elem);

    arma::mat h(indices_in[0].n_elem, indices_in[0].n_elem, arma::fill::zeros);

    sigma4 = sigma3 % sigma;
    h(mu_indices, mu_indices) += arma::diagmat(-1/sigma2);
    h(sigma_indices, sigma_indices) += arma::diagmat(1/sigma2 - 3*x2/sigma4);
    h(mu_indices, sigma_indices) = arma::diagmat(-2*sigma % x /sigma4);
    h(sigma_indices, mu_indices) = h(mu_indices, sigma_indices);

  }

  void dconstraints() {

  }

  void outcomes() {

    matrices.resize(1);
    matrices[0] = jacob;

    cubes.resize(1);
    cubes[0] = jacob2;

  }

};

normal* choose_normal(const Rcpp::List& trans_setup) {

  normal* mytrans = new normal();

  std::vector<arma::uvec> indices_in = trans_setup["indices_in"];
  std::vector<arma::uvec> indices_out = trans_setup["indices_out"];
  arma::vec y = trans_setup["y"];

  // arma::vec newgrad(indices_in[0].n_elem);

  mytrans->indices_in = indices_in;
  mytrans->indices_out = indices_out;
  mytrans->y = y;
  // mytrans->newgrad = newgrad;

  return mytrans;

}
