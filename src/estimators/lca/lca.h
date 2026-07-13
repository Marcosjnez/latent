/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 13/07/2026
 */

/*
 * Latent class analysis with expanded class probabilities
 */

class lca: public estimators {

public:

  int S;
  int I;
  arma::vec weights;
  arma::vec logweights;
  arma::uvec indices_classes, indices_classloglik;
  arma::mat classes;
  arma::mat logclasses;
  arma::mat classloglik;
  arma::mat joint_classloglik;
  arma::vec loglik_case;
  arma::vec logliks;
  arma::mat posterior, logposterior;
  double loss, N;

  void param(arguments_optim& x) {

    classes = arma::reshape(x.transparameters.elem(indices_classes), S, I);
    classloglik = arma::reshape(x.transparameters.elem(indices_classloglik), S, I);

    logclasses = arma::trunc_log(classes);

    for (int s = 0; s < S; ++s) {

      joint_classloglik.row(s) = classloglik.row(s) + logclasses.row(s);

      double max_vector = joint_classloglik.row(s).max();

      loglik_case(s) = max_vector +
        arma::trunc_log(
          arma::accu(
            arma::trunc_exp(joint_classloglik.row(s) - max_vector)
          )
        );

      logposterior.row(s) = joint_classloglik.row(s) - loglik_case(s);
      posterior.row(s) = arma::trunc_exp(logposterior.row(s));
      logliks(s) = weights(s) * loglik_case(s);

    }

  }

  void F(arguments_optim& x) {

    loss = -arma::accu(logliks);
    x.f += loss;

  }

  void G(arguments_optim& x) {

    arma::mat df_dclasses(S, I, arma::fill::zeros);
    arma::mat df_dclassloglik(S, I, arma::fill::zeros);

    for (int s = 0; s < S; ++s) {
      for (int i = 0; i < I; ++i) {

        df_dclasses(s, i) -=
          arma::trunc_exp(logweights(s) + classloglik(s, i) - loglik_case(s));

        df_dclassloglik(s, i) -=
          arma::trunc_exp(logweights(s) + logposterior(s, i));

      }
    }

    x.grad.elem(indices_classes) += arma::vectorise(df_dclasses);
    x.grad.elem(indices_classloglik) += arma::vectorise(df_dclassloglik);

  }

  void dG(arguments_optim& x) {

    arma::mat dclasses = arma::reshape(x.dtransparameters.elem(indices_classes), S, I);
    arma::mat dclassloglik = arma::reshape(x.dtransparameters.elem(indices_classloglik), S, I);

    arma::mat ddf_dclasses(S, I, arma::fill::zeros);
    arma::mat ddf_dclassloglik(S, I, arma::fill::zeros);

    arma::mat dlogclasses(S, I, arma::fill::zeros);
    arma::mat djoint_classloglik(S, I, arma::fill::zeros);
    arma::vec dloglik_case(S, arma::fill::zeros);
    arma::mat dlogposterior(S, I, arma::fill::zeros);

    for (int s = 0; s < S; ++s) {
      for (int i = 0; i < I; ++i) {
        dlogclasses(s, i) = dclasses(s, i) / classes(s, i);
        djoint_classloglik(s, i) = dclassloglik(s, i) + dlogclasses(s, i);
      }
    }

    for (int s = 0; s < S; ++s) {
      double acc = 0.0;
      for (int i = 0; i < I; ++i) {
        acc += posterior(s, i) * djoint_classloglik(s, i);
      }
      dloglik_case(s) = acc;
    }

    for (int s = 0; s < S; ++s) {
      for (int i = 0; i < I; ++i) {

        dlogposterior(s, i) = djoint_classloglik(s, i) - dloglik_case(s);

        double term_classes =
          arma::trunc_exp(logweights(s) + classloglik(s, i) - loglik_case(s));

        double term_classloglik =
          arma::trunc_exp(logweights(s) + logposterior(s, i));

        ddf_dclasses(s, i) -=
          term_classes * (dclassloglik(s, i) - dloglik_case(s));

        ddf_dclassloglik(s, i) -=
          term_classloglik * dlogposterior(s, i);

      }
    }

    x.dgrad.elem(indices_classes) += arma::vectorise(ddf_dclasses);
    x.dgrad.elem(indices_classloglik) += arma::vectorise(ddf_dclassloglik);

  }

  void outcomes(arguments_optim& x) {

    doubles.resize(7);
    double loglik = -loss;
    doubles[0] =  loss;
    doubles[1] =  0.00;
    doubles[2] =  0.00;
    doubles[3] =  loglik;
    doubles[4] =  0.00;
    doubles[5] =  0.00;
    doubles[6] =  0.00;

    vectors.resize(1);
    vectors[0] = loglik_case;

    matrices.resize(1);
    matrices[0] = logposterior;

  }

};

lca* choose_lca(const Rcpp::List& estimator_setup) {

  lca* myestimator = new lca();

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  int S = estimator_setup["S"];
  int I = estimator_setup["I"];
  arma::vec weights = estimator_setup["weights"];

  double N = arma::accu(weights);

  arma::mat classes(S, I, arma::fill::zeros);
  arma::mat classloglik(S, I, arma::fill::zeros);
  arma::mat logclasses(S, I, arma::fill::zeros);
  arma::mat joint_classloglik(S, I, arma::fill::zeros);
  arma::mat posterior(S, I, arma::fill::zeros);
  arma::mat logposterior(S, I, arma::fill::zeros);
  arma::vec logliks(S, arma::fill::zeros);
  arma::vec loglik_case(S, arma::fill::zeros);
  arma::vec logweights = arma::trunc_log(weights);

  myestimator->S = S;
  myestimator->I = I;
  myestimator->weights = weights;
  myestimator->logweights = logweights;
  myestimator->indices_classes = indices[0];
  myestimator->indices_classloglik = indices[1];
  myestimator->classes = classes;
  myestimator->classloglik = classloglik;
  myestimator->logclasses = logclasses;
  myestimator->joint_classloglik = joint_classloglik;
  myestimator->posterior = posterior;
  myestimator->logposterior = logposterior;
  myestimator->logliks = logliks;
  myestimator->loglik_case = loglik_case;
  myestimator->N = N;
  myestimator->loss = 0.0;

  return myestimator;

}
