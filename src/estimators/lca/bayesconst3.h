/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 24/08/2025
 */

/*
 * Constant prior 3 (for the variance of the items that belongs to
 * the same class membership)
 */

class bayesconst3: public estimators {

public:

  int K;
  double alpha, constant, prod_vars;
  arma::vec vars, varshat, sds;

  void param() {

    sds = transparameters;
    vars = sds % sds;
    prod_vars = arma::prod(vars);
    constant = alpha/K;

  }

  void F() {

    f = 0.5*constant * (arma::trunc_log(prod_vars) + arma::accu(varshat/vars));

  }

  void G() {

    grad = -constant * (varshat/(vars % sds) - 1/sds);

  }

  void dG() {

    dg.set_size(transparameters.n_elem); dg.zeros();

  }

  void H() {

    arma::vec d2 = -constant * (1/vars - 3*varshat/(vars % vars));
    hess = diagmat(d2);

  }

  void E() { // Update the parameter estimates

  }

  void M() { // Update the posterior probabilities

  }

  void outcomes() {

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
