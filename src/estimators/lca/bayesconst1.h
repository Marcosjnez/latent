/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 15/02/2026
 */

/*
 * Constant prior 1 (for the class membership probabilities)
 */

class bayesconst1: public estimators {

public:

  int K;
  int U;
  double alpha, constant, N;
  arma::uvec indices;
  arma::vec trans, logtrans;
  arma::vec constant_logtrans;

  void param(arguments_optim& x) {

    trans = x.transparameters(indices);
    logtrans = arma::trunc_log(trans);
    constant = alpha/(K*U + 0.00);
    constant_logtrans = constant*logtrans;

  }

  void F(arguments_optim& x) {

    f = arma::accu(constant_logtrans);
    x.f -= f/N;

  }

  void G(arguments_optim& x) {

    x.grad.elem(indices) -= constant/trans/N;

  }

  void dG(arguments_optim& x) {

    arma::vec dtrans = x.dtransparameters(indices);
    x.dgrad.elem(indices) += constant*dtrans/(trans % trans)/N;

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

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  int K = estimator_setup["K"];
  int U = estimator_setup["U"];
  double alpha = estimator_setup["alpha"];
  double N = estimator_setup["N"];

  myestimator->indices = indices[0];
  myestimator->K = K;
  myestimator->U = U;
  myestimator->alpha = alpha;
  myestimator->N = N;

  return myestimator;

}
