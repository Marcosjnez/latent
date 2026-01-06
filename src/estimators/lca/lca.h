/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 27/10/2025
 */

/*
 * Latent class analysis with expanded class probabilities
 */

class lca: public estimators {

public:

  int S; // rows of Y (number of response patterns)
  int J; // cols of Y (number of items)
  int I; // Number of latent classes
  arma::uvec indices_classes;
  arma::uvec indices_itemloglik;
  arma::uvec indices_theta;
  arma::vec weights; // Number of repetitions of each response pattern
  arma::vec logweights;
  arma::mat classes;
  arma::mat logclasses;
  arma::vec logliks;
  arma::vec loglik_case;
  arma::mat jointlogp;
  arma::cube itemloglik, gitemloglik, ditemloglik;
  arma::mat latentloglik;
  arma::mat posterior, logposterior;
  std::vector<arma::uvec> hess_indices;

  void param(arguments_optim& x) {

    arma::vec cl = x.transparameters.elem(indices_classes);
    classes = arma::reshape(cl, S, I);
    logclasses = arma::trunc_log(classes);

    arma::vec values = x.transparameters.elem(indices_itemloglik);
    std::memcpy(itemloglik.memptr(), values.memptr(), sizeof(double) * values.n_elem);

    latentloglik.zeros();

    for(int s=0; s < S; ++s) {
      for(int i=0; i < I; ++i) {
        for(int j=0; j < J; ++j) {
          latentloglik(s, i) += itemloglik(s, j, i);
        }
      }
      jointlogp.row(s) = latentloglik.row(s) + logclasses.row(s);
      double max_vector = jointlogp.row(s).max();
      loglik_case(s) = max_vector +
        arma::trunc_log(arma::accu(arma::trunc_exp(jointlogp.row(s) - max_vector)));
      logposterior.row(s) = jointlogp.row(s) - loglik_case(s);
      logliks(s) = weights(s) * loglik_case(s);
    }

  }

  void F(arguments_optim& x) {

    f = arma::accu(logliks);
    x.f -= f;

  }

  void G(arguments_optim& x) {

    arma::mat gclasses(S, I, arma::fill::zeros);
    arma::mat glatentloglik(S, I, arma::fill::zeros);
    for(int s=0; s < S; ++s) {
      for(int i=0; i < I; ++i) {
        gclasses(s, i) -= arma::trunc_exp(logweights[s] + latentloglik(s, i) - loglik_case(s));
        glatentloglik(s, i) -= arma::trunc_exp(logweights[s] + logposterior(s, i));
      }
    }

    // Replicate glatentloglik across slices:
    for (arma::uword k = 0; k < I; ++k) {
      gitemloglik.slice(k).each_col() = glatentloglik.col(k);
    }

    x.grad.elem(indices_classes) += arma::vectorise(gclasses);
    x.grad.elem(indices_itemloglik) += arma::vectorise(gitemloglik);

  }

  void dG(arguments_optim& x) {

    arma::mat dclasses = arma::reshape(x.dtransparameters(indices_classes), S, I);
    arma::vec dvalues = x.dtransparameters(indices_itemloglik);
    std::memcpy(ditemloglik.memptr(), dvalues.memptr(), sizeof(double) * dvalues.n_elem);

    arma::mat dgclasses;
    arma::cube dgitems(S, J, I, arma::fill::zeros);
    x.dgrad.elem(indices_classes) += arma::vectorise(dgclasses);
    x.dgrad.elem(indices_itemloglik) += arma::vectorise(dgitems);

  }

  void outcomes(arguments_optim& x) {

    doubles.resize(2);
    doubles[0] = f;
    doubles[1] = -f;

    vectors.resize(3);
    vectors[0] = loglik_case;
    vectors[1] = logliks;
    vectors[2] = weights;

    matrices.resize(5);
    matrices[0] = latentloglik;
    matrices[1] = logposterior;
    matrices[2] = jointlogp;
    matrices[3] = classes;
    matrices[4] = logclasses;

    cubes.resize(2);
    cubes[0] = itemloglik;
    cubes[1] = gitemloglik;

  }

};

lca* choose_lca(const Rcpp::List& estimator_setup) {

  lca* myestimator = new lca();

  int S = estimator_setup["S"];
  int J = estimator_setup["J"];
  int I = estimator_setup["I"];
  std::vector<arma::uvec> indices = estimator_setup["indices"];
  arma::vec weights = estimator_setup["weights"];
  std::vector<arma::uvec> hess_indices = estimator_setup["hess_indices"];

  arma::mat classes(S, I);
  arma::vec logliks(S, arma::fill::zeros);
  arma::vec loglik_case(S, arma::fill::zeros);
  arma::cube itemloglik(S, J, I, arma::fill::zeros);
  arma::cube gitemloglik(S, J, I, arma::fill::zeros);
  arma::cube ditemloglik(S, J, I, arma::fill::zeros);
  arma::mat latentloglik(S, I, arma::fill::zeros);
  arma::mat jointlogp(S, I, arma::fill::zeros);
  arma::mat posterior(S, I, arma::fill::zeros);
  arma::mat logposterior(S, I, arma::fill::zeros);
  arma::vec logweights = arma::trunc_log(weights);
  arma::uvec indices_classes = indices[1];
  arma::uvec indices_itemloglik = indices[2];
  arma::uvec indices_theta = indices[3];

  myestimator->S = S;
  myestimator->J = J;
  myestimator->I = I;
  myestimator->weights = weights;
  myestimator->logweights = logweights;
  myestimator->indices = indices;
  myestimator->indices_classes = indices_classes;
  myestimator->indices_itemloglik = indices_itemloglik;
  myestimator->indices_theta = indices_theta;
  myestimator->hess_indices = hess_indices;

  myestimator->classes = classes;
  myestimator->logliks = logliks;
  myestimator->loglik_case = loglik_case;
  myestimator->itemloglik = itemloglik;
  myestimator->gitemloglik = gitemloglik;
  myestimator->ditemloglik = ditemloglik;
  myestimator->latentloglik = latentloglik;
  myestimator->jointlogp = jointlogp;
  myestimator->posterior = posterior;
  myestimator->logposterior = logposterior;

  return myestimator;

}
