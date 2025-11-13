/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 31/08/2025
 */

/*
 * Constant prior 2 (for the conditional item probabilities of multinomial items)
 */

class bayesconst2: public estimators {

public:

  double alpha;
  int K;
  arma::vec constant, pihat;
  arma::vec trans, logtrans;
  arma::vec constant_logtrans;

  void param(arguments_optim& x) {

    trans = x.transparameters(indices[0]);
    logtrans = arma::trunc_log(trans);
    // double K = logtrans.n_elem + 0.00;
    constant = pihat * (alpha/(K + 0.00));
    constant_logtrans = constant % logtrans;

  }

  void F(arguments_optim& x) {

    f = arma::accu(constant_logtrans);
    x.f -= f;

  }

  void G(arguments_optim& x) {

    x.grad.elem(indices[0]) -= constant/trans;

  }

  void dG(arguments_optim& x) {

    // dg.set_size(transparameters.n_elem); dg.zeros();

  }

  void H(arguments_optim& x) {

    arma::vec d2constant_logtrans = constant/(trans % trans);
    x.hess(indices[0], indices[0]) += diagmat(d2constant_logtrans);

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

bayesconst2* choose_bayesconst2(const Rcpp::List& estimator_setup) {

  bayesconst2* myestimator = new bayesconst2();

  double alpha = estimator_setup["alpha"];
  int K = estimator_setup["K"];
  arma::vec pihat = estimator_setup["pihat"];
  std::vector<arma::uvec> indices = estimator_setup["indices"];

  myestimator->K = K;
  myestimator->alpha = alpha;
  myestimator->pihat = pihat;
  myestimator->indices = indices;

  return myestimator;

}
