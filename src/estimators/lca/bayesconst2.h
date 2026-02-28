/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 15/02/2026
 */

/*
 * Constant prior 2 (for the conditional item probabilities of multinomial items)
 */

class bayesconst2: public estimators {

public:

  int K;
  double alpha, N;
  arma::uvec indices;
  arma::vec constant, pihat;
  arma::vec trans, logtrans;
  arma::vec constant_logtrans;

  void param(arguments_optim& x) {

    trans = x.transparameters(indices);
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

    x.grad.elem(indices) -= constant/trans/N;

  }

  void dG(arguments_optim& x) {

    arma::vec dtrans = x.dtransparameters(indices);
    x.dgrad.elem(indices) += constant % dtrans/(trans % trans)/N;

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
  int K = estimator_setup["K"];
  double alpha = estimator_setup["alpha"];
  double N = estimator_setup["N"];
  arma::vec pihat = estimator_setup["pihat"];

  myestimator->indices = indices[0];
  myestimator->K = K;
  myestimator->alpha = alpha;
  myestimator->pihat = pihat;
  myestimator->N = N;

  return myestimator;

}
