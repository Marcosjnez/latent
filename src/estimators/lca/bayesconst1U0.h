/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 24/09/2025
 */

/*
 * Constant prior 1 (for the class membership probabilities)
 */

class bayesconst1U0: public estimators {

public:

  int K;
  int U;
  double alpha, constant;
  arma::vec trans, logtrans;
  arma::vec constant_logtrans;

  void param(arguments_optim& x) {

    trans = x.transparameters(indices[0]);
    logtrans = arma::trunc_log(trans);
    constant = alpha/(K*U + 0.00);
    constant_logtrans = constant*logtrans;

  }

  void F(arguments_optim& x) {

    f = arma::accu(constant_logtrans);
    x.f -= f;

  }

  void G(arguments_optim& x) {

    x.grad.elem(indices[0]) -= constant/trans;

  }

  void dG(arguments_optim& x) {

    dg.set_size(transparameters.n_elem); dg.zeros();

  }

  void H(arguments_optim& x) {

    arma::vec d2constant_logtrans = constant/(trans % trans);
    x.hess(indices[0], indices[0]) += diagmat(d2constant_logtrans);

  }

  void E(arguments_optim& x) { // Update the parameter estimates

    x.loglik = f;

  }

  void M(arguments_optim& x) { // Update the posterior probabilities

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

bayesconst1U0* choose_bayesconst1U0(const Rcpp::List& estimator_setup) {

  bayesconst1U0* myestimator = new bayesconst1U0();

  int K = estimator_setup["K"];
  int U = estimator_setup["U"];
  double alpha = estimator_setup["alpha"];
  std::vector<arma::uvec> indices = estimator_setup["indices"];

  myestimator->K = K;
  myestimator->U = U;
  myestimator->alpha = alpha;
  myestimator->indices = indices;

  return myestimator;

}
