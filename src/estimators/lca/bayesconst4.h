/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 13/07/2026
 */

/*
 * Constant prior 4 (for the error variance-covariance matrix of the
 * multivariate gaussian items)
 */

class bayesconst4: public estimators {

public:

  int K, J;
  double alpha, constant, prod_vars, N, logdetSigma, loss;
  arma::uvec indices_Sigma;
  arma::mat Sigma, invSigma, D;

  void param(arguments_optim& x) {

    Sigma = arma::reshape(x.transparameters.elem(indices_Sigma), J, J);

    if(!Sigma.is_sympd()) {
      arma::vec eigval;
      arma::mat eigvec;
      eig_sym(eigval, eigvec, Sigma);
      arma::vec d = arma::clamp(eigval, 0.00001, eigval.max());
      Sigma = eigvec * arma::diagmat(d) * eigvec.t();
    }

    constant = alpha / (K + 0.00);
    logdetSigma = arma::log_det_sympd(Sigma);
    invSigma = arma::inv_sympd(Sigma);

  }

  void F(arguments_optim& x) {

    loss = 0.5 * constant * (logdetSigma + arma::trace(D * invSigma));
    x.f += loss;

  }

  void G(arguments_optim& x) {

    arma::mat gSigma =
      0.5 * constant * (invSigma - invSigma * D * invSigma);

    x.grad.elem(indices_Sigma) += arma::vectorise(gSigma);

  }

  void dG(arguments_optim& x) {

    arma::mat dSigma = arma::reshape(x.dtransparameters.elem(indices_Sigma), J, J);
    arma::mat dInvSigma = -invSigma * dSigma * invSigma;

    arma::mat dgSigma =
      0.5 * constant *
      (dInvSigma - dInvSigma * D * invSigma - invSigma * D * dInvSigma);

    x.dgrad.elem(indices_Sigma) += arma::vectorise(dgSigma);

  }

  void outcomes(arguments_optim& x) {

    doubles.resize(1);
    doubles[0] = loss;

    names_doubles.resize(1);
    names_doubles[0] = "penalty";

  }

};

bayesconst4* choose_bayesconst4(const Rcpp::List& estimator_setup) {

  bayesconst4* myestimator = new bayesconst4();

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  double alpha = estimator_setup["alpha"];
  int K = estimator_setup["K"];
  int J = estimator_setup["J"];
  arma::mat D = estimator_setup["D"];

  myestimator->indices_Sigma = indices[0];
  myestimator->alpha = alpha;
  myestimator->K = K;
  myestimator->J = J;
  myestimator->D = D;

  return myestimator;

}
