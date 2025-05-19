/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 16/05/2025
 */

/*
 * Combine the likelihoods of Latent class models
 */

class latentloglik_combination: public estimators {

public:

  int nclasses; // Number of latent classes
  arma::vec lclasses, classes;
  arma::uvec indices_classes, indices_target_classes;
  arma::mat jointlogp;
  arma::vec loglik_case, logliks;
  int S;

  void param() {

    latentpars(indices_target_classes) = transparameters(indices_classes);
    loglatentpars = trunc_log(latentpars);

    // for (arma::uword j = 0; j < lclasses.n_elem; ++j) {
    //     Rprintf("%g ", lclasses[j]);
    //   }

    for(int s=0; s < S; ++s) {
      jointlogp.row(s) = latentloglik.row(s) + loglatentpars.t();
      double max_vector = jointlogp.row(s).max();
      loglik_case[s] = max_vector + arma::trunc_log(arma::accu(arma::trunc_exp(jointlogp.row(s) - max_vector)));
      logliks[s] = n[s] * loglik_case[s];
    }

    latentloglik.zeros();

    // arma::vec v = logliks;
    // for (arma::uword j = 0; j < v.n_elem; ++j) {
    //   Rprintf("%g \n", v[j]);
    // }
  }

  void F() {

    f = -arma::accu(logliks);

    // Rprintf("%g \n", f);
  }

  void G() {

    // Rf_error("58");
    arma::vec gclasses(nclasses, arma::fill::zeros);
    for(int s=0; s < S; ++s) {
      gclasses += -n[s]*arma::trunc_exp(-loglik_case[s] + latentloglik.row(s).t());
    }

    // for (arma::uword j = 0; j < gclasses.n_elem; ++j) {
    //   Rprintf("%g \n", gclasses[j]);
    // }

    grad = gclasses(indices_target_classes);

  }

  void dG() {

    dg.set_size(parameters.n_elem); dg.zeros();

  }

  void H() {

    hessian.set_size(parameters.n_elem, parameters.n_elem); hessian.zeros();

  }

  void E() { // Update the parameter estimates and loglik of each observation

    // Rf_error("EM algorithm not available");
    // latentloglik.zeros();

  }

  void M() { // Update the posterior probabilities

    // Rf_error("EM algorithm not available");

  }

  void outcomes() {

    vectors.resize(4);
    vectors[0] = latentpars;
    vectors[1] = loglatentpars;
    vectors[2] = loglik_case;
    vectors[3] = n;

    matrices.resize(2);
    matrices[0] = posterior;
    matrices[1] = latentloglik;

  }

};

latentloglik_combination* choose_latentloglik_combination(Rcpp::List estimator_setup) {

  latentloglik_combination* myestimator = new latentloglik_combination();

  arma::vec classes = estimator_setup["classes"];
  arma::uvec indices_classes = estimator_setup["indices_classes"];
  arma::uvec indices_target_classes = estimator_setup["indices_target_classes"];
  int S = estimator_setup["S"];
  int nclasses = estimator_setup["nclasses"];
  arma::vec n = estimator_setup["n"];
  arma::uvec indices = estimator_setup["indices"];

  arma::vec loglik_case(S, arma::fill::zeros);
  arma::mat latentloglik(S, nclasses, arma::fill::zeros);
  arma::vec latentpars(nclasses, arma::fill::zeros);
  arma::vec loglatentpars(nclasses, arma::fill::zeros);
  arma::mat jointlogp(S, nclasses, arma::fill::zeros);
  arma::vec logliks(S, arma::fill::zeros);

  myestimator->nclasses = nclasses;
  myestimator->classes = classes;
  myestimator->lclasses = classes;
  myestimator->indices_classes = indices_classes;
  myestimator->indices_target_classes = indices_target_classes;
  myestimator->S = S;
  myestimator->n = n;
  myestimator->indices = indices;
  myestimator->loglik_case = loglik_case;
  myestimator->logliks = logliks;
  // myestimator->posterior = posterior;
  myestimator->latentloglik = latentloglik;
  myestimator->latentpars = latentpars;
  myestimator->loglatentpars = loglatentpars;
  myestimator->jointlogp = jointlogp;
  myestimator->evalEM = false;

  return myestimator;

}
