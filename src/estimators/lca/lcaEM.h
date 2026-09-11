/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 11/09/2026
 */

// Latent class analysis for Expectation-Maximization

class lcaEM: public estimators {

public:

  int S;
  int I;

  arma::vec weights;

  arma::uvec indices_classes, indices_classloglik;

  arma::mat classes;
  arma::mat logclasses;
  arma::mat classloglik;
  arma::mat joint_classloglik;

  arma::vec loglik_case;

  arma::mat posterior;
  arma::mat logposterior;

  // Weighted posterior probabilities frozen during the E step:
  arma::mat weighted_posterior_em;

  double loss;
  double q_loss;


  void param(arguments_optim& x) override {

    classes = arma::reshape(x.transparameters.elem(indices_classes), S, I);
    classloglik = arma::reshape(x.transparameters.elem(indices_classloglik), S, I);

    logclasses = arma::trunc_log(classes);
    joint_classloglik = classloglik + logclasses;

  }


  // E step: freeze posterior probabilities
  void E(arguments_optim& x) override {

    weighted_posterior_em = posterior;
    weighted_posterior_em.each_col() %= weights;
    x.posterior = posterior;

  }


  // Observed-data negative log-likelihood
  void observed_F(arguments_optim& x) override {

    loss = 0.0;

    for(int s=0; s < S; ++s) {

      double max_vector = joint_classloglik.row(s).max();

      loglik_case(s) =
        max_vector +
        arma::trunc_log(
          arma::accu(
            arma::trunc_exp(
              joint_classloglik.row(s) - max_vector
            )
          )
        );

      logposterior.row(s) =
        joint_classloglik.row(s) - loglik_case(s);

      posterior.row(s) =
        arma::trunc_exp(logposterior.row(s));

      loss -= weights(s)*loglik_case(s);

    }

    x.f += loss;

  }


  // Gradient of the observed-data negative log-likelihood
  void observed_G(arguments_optim& x) override {

    arma::mat weighted_posterior = posterior;
    weighted_posterior.each_col() %= weights;

    x.grad.elem(indices_classes) -=
      arma::vectorise(weighted_posterior / classes);

    x.grad.elem(indices_classloglik) -=
      arma::vectorise(weighted_posterior);

  }


  // Expected complete-data negative log-likelihood:
  //
  // Q = - sum_s w_s sum_i tau_si *
  //       [log(class_si) + classloglik_si]
  //
  // weighted_posterior_em is kept fixed during the complete M step.
  void F(arguments_optim& x) override {

    q_loss =
      -arma::accu(weighted_posterior_em % joint_classloglik);

    x.f += q_loss;

  }


  // Gradient of the Q-function
  void G(arguments_optim& x) override {

    x.grad.elem(indices_classes) -=
      arma::vectorise(weighted_posterior_em / classes);

    x.grad.elem(indices_classloglik) -=
      arma::vectorise(weighted_posterior_em);

  }


  // Differential of the gradient of the Q-function
  void dG(arguments_optim& x) override {

    arma::mat dclasses =
      arma::reshape(
        x.dtransparameters.elem(indices_classes), S, I
      );

    x.dgrad.elem(indices_classes) +=
      arma::vectorise(
        weighted_posterior_em % dclasses /
          (classes % classes)
      );

    // The Q-function is linear in classloglik, so its contribution to
    // dgrad is zero. Second derivatives with respect to the conditional
    // parameters are subsequently obtained by the likelihood transformations
    // in update_dgrad().

  }


  void outcomes(arguments_optim& x) override {

    doubles.resize(2);
    doubles[0] = loss;
    doubles[1] = -loss;

    names_doubles.resize(2);
    names_doubles[0] = "loss";
    names_doubles[1] = "loglik";

    vectors.resize(1);
    vectors[0] = loglik_case;

    names_vectors.resize(1);
    names_vectors[0] = "loglik_case";

    matrices.resize(1);
    matrices[0] = logposterior;

    names_matrices.resize(1);
    names_matrices[0] = "logposterior";

  }

};

lcaEM* choose_lcaEM(const Rcpp::List& estimator_setup) {

  lcaEM* myestimator = new lcaEM();

  std::vector<arma::uvec> indices = estimator_setup["indices"];

  int S = estimator_setup["S"];
  int I = estimator_setup["I"];

  arma::vec weights = estimator_setup["weights"];

  arma::mat classes(S, I, arma::fill::zeros);
  arma::mat classloglik(S, I, arma::fill::zeros);
  arma::mat logclasses(S, I, arma::fill::zeros);
  arma::mat joint_classloglik(S, I, arma::fill::zeros);

  arma::mat posterior(S, I, arma::fill::zeros);
  arma::mat logposterior(S, I, arma::fill::zeros);
  arma::mat weighted_posterior_em(S, I, arma::fill::zeros);

  arma::vec loglik_case(S, arma::fill::zeros);

  myestimator->S = S;
  myestimator->I = I;

  myestimator->weights = weights;

  myestimator->indices_classes = indices[0];
  myestimator->indices_classloglik = indices[1];

  myestimator->classes = classes;
  myestimator->classloglik = classloglik;
  myestimator->logclasses = logclasses;
  myestimator->joint_classloglik = joint_classloglik;

  myestimator->posterior = posterior;
  myestimator->logposterior = logposterior;
  myestimator->weighted_posterior_em = weighted_posterior_em;

  myestimator->loglik_case = loglik_case;

  myestimator->loss = 0.0;
  myestimator->q_loss = 0.0;

  return myestimator;

}
