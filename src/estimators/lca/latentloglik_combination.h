/*
 * Author: Marcos Jimenez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 03/02/2025
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

    // NEEDS latentloglik
    // PROVIDES latentpars and jointlogp

    classes(indices_target_classes) = parameters(indices_classes);
    latentpars = trunc_log(classes);

    // for (arma::uword j = 0; j < lclasses.n_elem; ++j) {
    //     Rprintf("%g ", lclasses[j]);
    //   }

    for(int s=0; s < S; ++s) {
      jointlogp.row(s) = latentloglik.row(s) + latentpars.t();
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

    g.set_size(parameters.n_elem); g.zeros();
    arma::vec gclasses(nclasses, arma::fill::zeros);
    for(int s=0; s < S; ++s) {
      gclasses += -n[s]*arma::trunc_exp(-loglik_case[s] + latentloglik.row(s).t());
    }

    // arma::vec v = gclasses;
    // for (arma::uword j = 0; j < v.n_elem; ++j) {
    //   Rprintf("%g \n", v[j]);
    // }

    g = gclasses(indices_target_classes);

    // arma::vec v = dale;
    // for (arma::uword j = 0; j < v.n_elem; ++j) {
    //   Rprintf("%g \n", v[j]);
    // }
  }

  void dG() {

    dg.set_size(parameters.n_elem); dg.zeros();

  }

  void H() {

    hessian.set_size(parameters.n_elem, parameters.n_elem); hessian.zeros();

  }

  void E() { // Update the parameter estimates and loglik of each observation

    // Rf_error("EM algorithm not available");

  }

  void M() { // Update the posterior probabilities

    // Rf_error("EM algorithm not available");

  }

  void outcomes() {

  }

};
