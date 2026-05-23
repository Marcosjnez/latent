/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 15/02/2026
 */

/*
 * Constant prior 3 (for the variance of the items that belongs to
 * the same class membership)
 */

class bayesconst3: public estimators {

public:

  int K;
  double alpha, constant, prod_vars, N, loss;
  arma::uvec indices;
  arma::vec vars, varshat, sds, logvars;

  void param(arguments_optim& x) {

    sds = x.transparameters(indices);
    vars = sds % sds;
    logvars = arma::trunc_log(vars);
    constant = alpha/(K + 0.00);

  }

  void F(arguments_optim& x) {

    loss = 0.5*constant * (arma::accu(logvars) + arma::accu(varshat/vars));
    x.f += loss;

  }

  void G(arguments_optim& x) {

    x.grad.elem(indices) -= constant * (varshat/(vars % sds) - 1/sds);

  }

  void dG(arguments_optim& x) {

    arma::vec dsds = x.dtransparameters(indices);
    x.dgrad.elem(indices) += constant * (3 * varshat % dsds / (vars % vars) - dsds / vars);

  }

  void outcomes(arguments_optim& x) {

    doubles.resize(7);
    doubles[0] =  0.00;        // loss actual model
    doubles[1] =  0.00;        // loss independence model
    doubles[2] =  0.00;        // loss saturated model
    doubles[3] =  0.00;        // loglik actual model
    doubles[4] =  0.00;        // loglik independence model
    doubles[5] =  0.00;        // loglik saturated model
    doubles[6] =  loss;        // penalty

    // vectors.resize(1);
    //
    // matrices.resize(1);
    //
    // cubes.resize(1);

  }

};

bayesconst3* choose_bayesconst3(const Rcpp::List& estimator_setup) {

  bayesconst3* myestimator = new bayesconst3();

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  int K = estimator_setup["K"];
  double alpha = estimator_setup["alpha"];
  double N = estimator_setup["N"];
  arma::vec varshat = estimator_setup["varshat"];

  myestimator->indices = indices[0];
  myestimator->alpha = alpha;
  myestimator->K = K;
  myestimator->varshat = varshat;
  myestimator->N = N;

  return myestimator;

}
