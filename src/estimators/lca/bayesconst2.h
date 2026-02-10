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

  double alpha, N;
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
    x.f -= f/N;

  }

  void G(arguments_optim& x) {

    x.grad.elem(indices[0]) -= constant/trans/N;

  }

  void dG(arguments_optim& x) {

    arma::vec dtrans = x.dtransparameters(indices[0]);
    x.dgrad.elem(indices[0]) += constant % dtrans/(trans % trans)/N;

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

bayesconst2* choose_bayesconst2(const Rcpp::List& estimator_setup) {

  bayesconst2* myestimator = new bayesconst2();

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  double alpha = estimator_setup["alpha"];
  int K = estimator_setup["K"];
  arma::vec pihat = estimator_setup["pihat"];
  double N = estimator_setup["N"];

  myestimator->indices = indices;
  myestimator->K = K;
  myestimator->alpha = alpha;
  myestimator->pihat = pihat;
  myestimator->N = N;

  return myestimator;

}
