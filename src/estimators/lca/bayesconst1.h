/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 23/08/2025
 */

/*
 * Constant prior 1 (for the class membership probabilities)
 */

class bayesconst1: public estimators {

public:

  double constant;
  arma::vec logtrans;
  arma::vec constant_logtrans;

  void param() {

    logtrans = arma::trunc_log(transparameters);
    double K = logtrans.n_elem + 0.00;
    constant = 1.00/K;
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

    // vectors.resize(1);
    //
    // matrices.resize(1);
    //
    // cubes.resize(1);

  }

};

bayesconst1* choose_bayesconst1(const Rcpp::List& estimator_setup) {

  bayesconst1* myestimator = new bayesconst1();

  std::vector<arma::uvec> indices = estimator_setup["indices"];

  myestimator->indices = indices;

  return myestimator;

}
