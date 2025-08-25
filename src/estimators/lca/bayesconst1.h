/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 25/08/2025
 */

/*
 * Constant prior 1 (for the class membership probabilities)
 */

class bayesconst1: public estimators {

public:

  int K;
  double alpha, constant;
  arma::vec logtrans;
  arma::vec constant_logtrans;

  void param() {

    logtrans = arma::trunc_log(transparameters);
    constant = alpha/(K + 0.00);
    constant_logtrans = constant*logtrans;

  }

  void F() {

    f = -arma::accu(constant_logtrans);

  }

  void G() {

    grad.set_size(transparameters.n_elem); grad.zeros();
    grad = -constant/transparameters;

  }

  void dG() {

    dg.set_size(transparameters.n_elem); dg.zeros();

  }

  void H() {

    hess.set_size(transparameters.n_elem, transparameters.n_elem);
    hess.zeros();

    arma::vec d2constant_logtrans = constant/(transparameters % transparameters);
    hess = diagmat(d2constant_logtrans);

  }

  void E() { // Update the parameter estimates

  }

  void M() { // Update the posterior probabilities

  }

  void outcomes() {

    doubles.resize(1);
    doubles[0] = f;

    // vectors.resize(1);
    //
    // matrices.resize(1);
    //
    // cubes.resize(1);

  }

};

bayesconst1* choose_bayesconst1(const Rcpp::List& estimator_setup) {

  bayesconst1* myestimator = new bayesconst1();

  int K = estimator_setup["K"];
  double alpha = estimator_setup["alpha"];
  std::vector<arma::uvec> indices = estimator_setup["indices"];

  myestimator->K = K;
  myestimator->alpha = alpha;
  myestimator->indices = indices;

  return myestimator;

}
