/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 31/08/2025
 */

/*
 * Constant prior 3 (for the variance of the items that belongs to
 * the same class membership)
 */

class bayesconst3: public estimators {

public:

  int K;
  double alpha, constant, prod_vars;
  arma::vec vars, varshat, sds, logvars;

  void param(arguments_optim& x) {

    sds = x.transparameters(indices[0]);
    vars = sds % sds;
    logvars = arma::trunc_log(vars);
    // prod_vars = arma::prod(vars);
    constant = alpha/(K + 0.00);

  }

  void F(arguments_optim& x) {

    f = -0.5*constant * (arma::accu(logvars) + arma::accu(varshat/vars));
    x.f -= f;

  }

  void G(arguments_optim& x) {

    x.grad.elem(indices[0]) -= constant * (varshat/(vars % sds) - 1/sds);

  }

  void dG(arguments_optim& x) {

    // dg.set_size(transparameters.n_elem); dg.zeros();

  }

  void H(arguments_optim& x) {

    arma::vec d2 = -constant * (1/vars - 3*varshat/(vars % vars));
    x.hess(indices[0], indices[0]) += diagmat(d2);

  }

  void outcomes(arguments_optim& x) {

    doubles.resize(1);
    doubles[0] = f;

    // vectors.resize(1);
    //
    // matrices.resize(1);
    //
    // cubes.resize(1);

  }

};

bayesconst3* choose_bayesconst3(const Rcpp::List& estimator_setup) {

  bayesconst3* myestimator = new bayesconst3();

  double alpha = estimator_setup["alpha"];
  int K = estimator_setup["K"];
  arma::vec varshat = estimator_setup["varshat"];
  std::vector<arma::uvec> indices = estimator_setup["indices"];

  myestimator->alpha = alpha;
  myestimator->K = K;
  myestimator->varshat = varshat;
  myestimator->indices = indices;

  return myestimator;

}
