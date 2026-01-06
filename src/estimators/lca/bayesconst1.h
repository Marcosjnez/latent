/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 27/10/2025
 */

/*
 * Constant prior 1 (for the class membership probabilities)
 */

class bayesconst1: public estimators {

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

    arma::vec dtrans = x.dtransparameters(indices[0]);
    x.dgrad.elem(indices[0]) += constant*dtrans/(trans % trans);

  }

  void outcomes(arguments_optim& x) {

    doubles.resize(2);
    doubles[0] = f;
    doubles[1] = -f;

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
  int U = estimator_setup["U"];
  double alpha = estimator_setup["alpha"];
  std::vector<arma::uvec> indices = estimator_setup["indices"];

  myestimator->K = K;
  myestimator->U = U;
  myestimator->alpha = alpha;
  myestimator->indices = indices;

  return myestimator;

}
