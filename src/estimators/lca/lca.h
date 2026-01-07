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
  arma::vec weights; // Number of repetitions of each response pattern
  arma::vec logweights;
  arma::uvec indices_classes;
  arma::uvec indices_cubeloglik;
  arma::cube cubeloglik; // items loglik conditional on classes
  arma::mat classes;     // probabilities of classes
  arma::mat logclasses;
  arma::mat classloglik;
  arma::mat joint_classloglik;
  arma::vec loglik_case; // loglik of each response pattern
  arma::vec logliks; // accumulated loglik contribution of each response pattern
  arma::mat posterior, logposterior;

  void param(arguments_optim& x) {

    // Load the parameters:
    classes = arma::reshape(x.transparameters.elem(indices_classes), S, I);
    arma::vec fill = x.transparameters.elem(indices_cubeloglik);
    std::memcpy(cubeloglik.memptr(), fill.memptr(), sizeof(double)*fill.n_elem);

    // Find the logarithm likelihood by case:
    logclasses = arma::trunc_log(classes);
    classloglik.zeros();
    for(int s=0; s < S; ++s) { // response patterns
      for(int i=0; i < I; ++i) { // classes
        for(int j=0; j < J; ++j) { // items
          classloglik(s, i) += cubeloglik(s, j, i);
        }
      }
      joint_classloglik.row(s) = classloglik.row(s) + logclasses.row(s);
      double max_vector = joint_classloglik.row(s).max();
      loglik_case(s) = max_vector +
        arma::trunc_log(arma::accu(arma::trunc_exp(joint_classloglik.row(s) - max_vector)));
      logposterior.row(s) = joint_classloglik.row(s) - loglik_case(s);
      logliks(s) = weights(s) * loglik_case(s);
    }

  }

  void F(arguments_optim& x) {

    f = arma::accu(logliks);
    x.f -= f;

  }

  void G(arguments_optim& x) {

    // Initialize derivative structures:
    arma::mat df_dclasses(S, I, arma::fill::zeros);
    arma::mat df_dclassloglik(S, I, arma::fill::zeros);
    arma::cube df_dcubeloglik(S, J, I, arma::fill::zeros);

    for(int s=0; s < S; ++s) { // response patterns
      for(int i=0; i < I; ++i) { // classes
        df_dclasses(s, i) -= arma::trunc_exp(logweights[s] + classloglik(s, i) - loglik_case(s));
        df_dclassloglik(s, i) -= arma::trunc_exp(logweights[s] + logposterior(s, i));
      }
    }

    // Replicate df_dclassloglik across items belonging to the same class:
    for (arma::uword k = 0; k < I; ++k) {
      df_dcubeloglik.slice(k).each_col() = df_dclassloglik.col(k);
    }

    x.grad.elem(indices_classes) += arma::vectorise(df_dclasses);
    x.grad.elem(indices_cubeloglik) += arma::vectorise(df_dcubeloglik);

  }

  void dG(arguments_optim& x) {

    // Load the directions:
    arma::mat dclasses = arma::reshape(x.dtransparameters(indices_classes), S, I);
    arma::vec fill = x.dtransparameters(indices_cubeloglik);
    arma::cube dcubeloglik(S, J, I, arma::fill::zeros);
    std::memcpy(dcubeloglik.memptr(), fill.memptr(), sizeof(double)*fill.n_elem);

    // Initialize the differential structures for the gradient:
    arma::mat ddf_dclasses(S, I, arma::fill::zeros);
    arma::mat ddf_dclassloglik(S, I, arma::fill::zeros);
    arma::cube ddf_dcubeloglik(S, J, I, arma::fill::zeros);

    // Differentials of intermediate quantities:
    arma::mat dclassloglik(S, I, arma::fill::zeros);
    arma::mat dlogclasses(S, I, arma::fill::zeros);
    arma::mat djoint_classloglik(S, I, arma::fill::zeros);
    arma::vec dloglik_case(S, arma::fill::zeros);
    arma::mat dlogposterior(S, I, arma::fill::zeros);

    // Differential of dclassloglik(s,i) = sum_j dcubeloglik(s,j,i)
    for (int s = 0; s < S; ++s) {
      for (int i = 0; i < I; ++i) {
        for (int j = 0; j < J; ++j) {
          dclassloglik(s, i) += dcubeloglik(s, j, i);
        }
      }
    }

    // Differentials of dlogclasses and djoint_classloglik:
    for (int s = 0; s < S; ++s) {
      for (int i = 0; i < I; ++i) {
        double C_si  = classes(s, i);
        double dC_si = dclasses(s, i);
        dlogclasses(s, i) = dC_si / C_si;
        djoint_classloglik(s, i) = dclassloglik(s, i) + dlogclasses(s, i);
      }
    }

    // Differential of dloglik_case(s) = sum_i posterior(s,i) * djoint_classloglik(s,i)
    for (int s = 0; s < S; ++s) {
      double acc = 0.0;
      for (int i = 0; i < I; ++i) {
        double post_si = arma::trunc_exp(logposterior(s, i));
        acc += post_si * djoint_classloglik(s, i);
      }
      dloglik_case(s) = acc;
    }

    // Differential of logposterior:
    for (int s = 0; s < S; ++s) {
      for (int i = 0; i < I; ++i) {
        dlogposterior(s, i) = djoint_classloglik(s, i) - dloglik_case(s);
      }
    }

    // Final differentials:
    for (int s = 0; s < S; ++s) {
      for (int i = 0; i < I; ++i) {

        double term1 = arma::trunc_exp(logweights[s] + classloglik(s, i) - loglik_case(s));
        ddf_dclasses(s, i) -= term1 * (dclassloglik(s, i) - dloglik_case(s));

        double term2 = arma::trunc_exp(logweights[s] + logposterior(s, i));
        ddf_dclassloglik(s, i) -= term2 * dlogposterior(s, i);
      }
    }

    // Replicate ddf_dclassloglik across items belonging to the same class:
    for (arma::uword k = 0; k < (arma::uword) I; ++k) {
      ddf_dcubeloglik.slice(k).each_col() = ddf_dclassloglik.col(k);
    }

    x.dgrad.elem(indices_classes)     += arma::vectorise(ddf_dclasses);
    x.dgrad.elem(indices_cubeloglik)  += arma::vectorise(ddf_dcubeloglik);

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
    matrices[0] = classloglik;
    matrices[1] = logposterior;
    matrices[2] = joint_classloglik;
    matrices[3] = classes;
    matrices[4] = logclasses;

    cubes.resize(1);
    cubes[0] = cubeloglik;

  }

};

lca* choose_lca(const Rcpp::List& estimator_setup) {

  lca* myestimator = new lca();

  int S = estimator_setup["S"];
  int J = estimator_setup["J"];
  int I = estimator_setup["I"];
  std::vector<arma::uvec> indices = estimator_setup["indices"];
  arma::vec weights = estimator_setup["weights"];

  arma::mat classes(S, I);
  arma::vec logliks(S, arma::fill::zeros);
  arma::vec loglik_case(S, arma::fill::zeros);
  arma::cube cubeloglik(S, J, I, arma::fill::zeros);
  arma::mat classloglik(S, I, arma::fill::zeros);
  arma::mat joint_classloglik(S, I, arma::fill::zeros);
  arma::mat posterior(S, I, arma::fill::zeros);
  arma::mat logposterior(S, I, arma::fill::zeros);
  arma::vec logweights = arma::trunc_log(weights);
  arma::uvec indices_classes = indices[1];
  arma::uvec indices_cubeloglik = indices[2];

  myestimator->S = S;
  myestimator->J = J;
  myestimator->I = I;
  myestimator->weights = weights;
  myestimator->logweights = logweights;
  myestimator->indices = indices;
  myestimator->indices_classes = indices_classes;
  myestimator->indices_cubeloglik = indices_cubeloglik;

  myestimator->classes = classes;
  myestimator->logliks = logliks;
  myestimator->loglik_case = loglik_case;
  myestimator->cubeloglik = cubeloglik;
  myestimator->classloglik = classloglik;
  myestimator->joint_classloglik = joint_classloglik;
  myestimator->posterior = posterior;
  myestimator->logposterior = logposterior;

  return myestimator;

}
