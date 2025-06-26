/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 09/06/2025
 */

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

class lgaussian:public transformations {

public:

  arma::uvec means_indices;
  arma::uvec sds_indices;
  arma::vec y;

  arma::vec x, means, sds, vars;

  void transform() {

    // Rprintf("43 \n");
    means = parameters.elem(means_indices);
    sds = parameters.elem(sds_indices);

    // Rprintf("47 \n");
    x = y - means;
    vars = sds % sds;
    // transparameters = ldnorm(y, means, sds);
    transparameters = -0.5 * x % x/vars - arma::trunc_log(SQRT2M_PI*sds);
    // Rprintf("50 \n");

  }

  void jacobian() {

    // jacob.set_size(transparameters.n_elem, parameters.n_elem);
    // jacob.zeros();
    //
    arma::vec exp_transparameters = arma::trunc_exp(transparameters);
    // arma::mat dmeans = arma::diagmat(x/vars % exp_transparameters);
    // arma::mat dsds = arma::diagmat(((x % x)-vars)/(vars % sds) %
    //   exp_transparameters);
    //
    // // Rprintf("68 \n");
    // jacob.cols(means_indices) = dmeans;
    // jacob.cols(sds_indices) = dsds;

    arma::vec dparameters(parameters.n_elem);
    dparameters(means_indices) = grad % x/vars % exp_transparameters;
    dparameters(sds_indices) = grad % ((x % x)-vars)/(vars % sds) %
      exp_transparameters;
    g = dparameters;
    // Rprintf("74 \n");
    // Rf_error("75");
  }

  void d2jacobian() {

    // jacob2.set_size(y.n_elem, parameters.n_elem, parameters.n_elem);
    // arma::vec d2means = transparameters %
    //   (y % y - 2*means*y + means % means - vars) / (vars % vars);
    // arma::vec d2sds = transparameters % ((1/vars - 3*x % x / (vars % vars)) +
    //   (x % x - vars) % (x % x - vars) / (vars % vars % vars));
    // arma::vec d2meanssds = transparameters % x % (x % x - 3*vars) /
    //   (vars % vars % sds);
    //
    // for (std::size_t i = 0; i < parameters.n_elem; ++i) {
    //   for (std::size_t j = i; j < parameters.n_elem; ++j) {
    //
    //   }
    //
    // }
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

lgaussian* choose_lgaussian(const Rcpp::List& trans_setup) {

  lgaussian* mytrans = new lgaussian();

  arma::uvec indices = trans_setup["indices"];
  arma::uvec target_indices = trans_setup["target_indices"];
  arma::vec y = trans_setup["y"];
  arma::uvec means_indices = trans_setup["means_indices"];
  arma::uvec sds_indices = trans_setup["sds_indices"];

  mytrans->indices = indices;
  mytrans->target_indices = target_indices;
  mytrans->y = y;
  mytrans->means_indices = means_indices;
  mytrans->sds_indices = sds_indices;

  return mytrans;

}
