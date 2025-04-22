/*
 * Author: Marcos Jimenez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 11/04/2025
 */

class estimators {

public:

  double f, loglik;
  arma::vec parameters, dparameters;
  arma::mat g, dg;
  arma::mat hessian, dparam_dS, modhessian;
  arma::mat posterior, latentloglik;
  arma::vec uniquenesses;
  // arma::mat lambda, phi, psi, model, residuals;
  int nhessian, nS;
  arma::vec n, latentpars; // P(X = c) // classes_hat
  std::vector<arma::mat> conditionals; // conditionals_hat

  std::vector<double> doubles;
  std::vector<arma::vec> vectors;
  std::vector<arma::mat> matrices;
  std::vector<std::vector<arma::mat>> list_matrices;

  arma::uvec indices;
  int q;
  bool evalE = false;

  virtual void param() = 0;

  virtual void F() = 0;

  virtual void G() = 0;

  virtual void dG() = 0;

  virtual void E() = 0;

  virtual void M() = 0;

  virtual void H() = 0;

  virtual void outcomes() = 0;

};

#include "estimators/efa/ml_efa.h"
#include "estimators/efa/uls_efa.h"
#include "estimators/efa/dwls_lt.h"

#include "estimators/rotation/cf.h"
#include "estimators/rotation/oblimin.h"
#include "estimators/rotation/geomin.h"
#include "estimators/rotation/varimax.h"
#include "estimators/rotation/varimin.h"
#include "estimators/rotation/target.h"
#include "estimators/rotation/xtarget.h"
#include "estimators/rotation/lclf.h"

#include "estimators/lca/lca_multinomial.h"
#include "estimators/lca/lca_multinomial_softmax.h"
#include "estimators/lca/lca_gaussian.h"
#include "estimators/lca/lca_gaussian_softmax.h"
#include "estimators/lca/latentloglik_combination.h"
#include "estimators/lca/latentloglik_combination_softmax.h"

#include "estimators/lreg/lreg.h"

#include "estimators/cfa/cfa_dwls.h"
#include "estimators/cfa/cfa_ml.h"

// Choose the estimator:

estimators* choose_estimator(Rcpp::List estimator_setup, estimators* xestimator) {

  estimators* criterion;
  std::string estimator = estimator_setup["estimator"];

  if (estimator == "ml_efa") {

    ml_efa* myestimator = new ml_efa();

    int p = estimator_setup["p"];
    int q = estimator_setup["q"];
    int nhessian = estimator_setup["nhessian"];
    int nS = estimator_setup["nS"];
    arma::uvec vp = consecutive(0, p-1);
    arma::uvec indices_diag_p = p*vp + vp;
    arma::uvec vq = consecutive(0, q-1);
    arma::uvec indices_diag_q = q*vq + vq;
    arma::mat Ip(p, p, arma::fill::eye); // Define identity
    arma::mat Iq(q, q, arma::fill::eye); // Define identity

    // Provide these:
    arma::uvec indices = estimator_setup["indices"];
    arma::uvec lambda_indices = estimator_setup["lambda_indices"];
    arma::uvec targetlambda_indices = estimator_setup["targetlambda_indices"];
    arma::uvec phi_indices = estimator_setup["phi_indices"];
    arma::uvec targetphi_indices = estimator_setup["targetphi_indices"];
    arma::uvec psi_indices = estimator_setup["psi_indices"];
    arma::uvec targetpsi_indices = estimator_setup["targetpsi_indices"];
    arma::mat R = estimator_setup["R"];

    myestimator->R = R;
    myestimator->p = p;
    myestimator->q = q;
    myestimator->nhessian = nhessian;
    myestimator->nS = nS;
    myestimator->indices_diag_p = indices_diag_p;
    myestimator->indices_diag_q = indices_diag_q;
    myestimator->Ip = Ip;
    myestimator->Iq = Iq;

    myestimator->indices = indices;
    myestimator->lambda_indices = lambda_indices;
    myestimator->phi_indices = phi_indices;
    myestimator->psi_indices = psi_indices;
    myestimator->targetlambda_indices = targetlambda_indices;
    myestimator->targetphi_indices = targetphi_indices;
    myestimator->targetpsi_indices = targetpsi_indices;

    criterion = myestimator;

  } else if (estimator == "uls_efa") {

    uls_efa* myestimator = new uls_efa();

    int p = estimator_setup["p"];
    int q = estimator_setup["q"];
    int nhessian = estimator_setup["nhessian"];
    int nS = estimator_setup["nS"];
    arma::uvec vp = consecutive(0, p-1);
    arma::uvec indices_diag_p = p*vp + vp;
    arma::uvec vq = consecutive(0, q-1);
    arma::uvec indices_diag_q = q*vq + vq;
    arma::mat Ip(p, p, arma::fill::eye); // Define identity
    arma::mat Iq(q, q, arma::fill::eye); // Define identity

    // Provide these:
    arma::uvec indices = estimator_setup["indices"];
    arma::uvec lambda_indices = estimator_setup["lambda_indices"];
    arma::uvec targetlambda_indices = estimator_setup["targetlambda_indices"];
    arma::uvec phi_indices = estimator_setup["phi_indices"];
    arma::uvec targetphi_indices = estimator_setup["targetphi_indices"];
    arma::uvec psi_indices = estimator_setup["psi_indices"];
    arma::uvec targetpsi_indices = estimator_setup["targetpsi_indices"];
    arma::mat R = estimator_setup["R"];

    myestimator->R = R;
    myestimator->p = p;
    myestimator->q = q;
    myestimator->nhessian = nhessian;
    myestimator->nS = nS;
    myestimator->indices_diag_p = indices_diag_p;
    myestimator->indices_diag_q = indices_diag_q;
    myestimator->Ip = Ip;
    myestimator->Iq = Iq;

    myestimator->indices = indices;
    myestimator->lambda_indices = lambda_indices;
    myestimator->phi_indices = phi_indices;
    myestimator->psi_indices = psi_indices;
    myestimator->targetlambda_indices = targetlambda_indices;
    myestimator->targetphi_indices = targetphi_indices;
    myestimator->targetpsi_indices = targetpsi_indices;

    criterion = myestimator;

  } else if (estimator == "dwls_lt") {

    dwls_lt* myestimator = new dwls_lt(); // Creating dwls_lt object dynamically

    int p = estimator_setup["p"];
    int q = estimator_setup["q"];
    int nhessian = estimator_setup["nhessian"];
    int nS = estimator_setup["nS"];
    arma::uvec vp = consecutive(0, p-1);
    arma::uvec indices_diag_p = p*vp + vp;
    arma::uvec vq = consecutive(0, q-1);
    arma::uvec indices_diag_q = q*vq + vq;
    arma::mat Ip(p, p, arma::fill::eye); // Define identity pxp
    arma::mat Iq(q, q, arma::fill::eye); // Define identity qxq

    // Provide these:
    arma::uvec indices = estimator_setup["indices"];
    arma::uvec lambda_indices = estimator_setup["lambda_indices"];
    arma::uvec targetlambda_indices = estimator_setup["targetlambda_indices"];
    arma::uvec phi_indices = estimator_setup["phi_indices"];
    arma::uvec targetphi_indices = estimator_setup["targetphi_indices"];
    arma::uvec psi_indices = estimator_setup["psi_indices"];
    arma::uvec targetpsi_indices = estimator_setup["targetpsi_indices"];
    arma::mat R = estimator_setup["R"];
    arma::mat W = estimator_setup["W"];
    arma::mat lambda(p, q, arma::fill::zeros);

    myestimator->R = R;
    myestimator->W = W;
    myestimator->lambda = lambda;
    myestimator->p = p;
    myestimator->q = q;
    myestimator->nhessian = nhessian;
    myestimator->nS = nS;
    myestimator->indices_diag_p = indices_diag_p;
    myestimator->indices_diag_q = indices_diag_q;
    myestimator->Ip = Ip;
    myestimator->Iq = Iq;

    myestimator->indices = indices;
    myestimator->lambda_indices = lambda_indices;
    myestimator->phi_indices = phi_indices;
    myestimator->psi_indices = psi_indices;
    myestimator->targetlambda_indices = targetlambda_indices;
    myestimator->targetphi_indices = targetphi_indices;
    myestimator->targetpsi_indices = targetpsi_indices;

    criterion = myestimator;

  } else if(estimator == "cfa_dwls") {

    cfa_dwls* myestimator = new cfa_dwls();

    arma::mat R = estimator_setup["R"];
    arma::mat lambda = estimator_setup["lambda"];
    arma::mat phi = estimator_setup["phi"];
    arma::mat psi = estimator_setup["psi"];
    arma::mat W = estimator_setup["W"];
    int nparam = estimator_setup["nparam"];
    arma::uvec lambda_indexes = estimator_setup["lambda_indexes"];
    arma::uvec phi_indexes = estimator_setup["phi_indexes"];
    arma::uvec psi_indexes = estimator_setup["psi_indexes"];
    arma::uvec target_indexes = estimator_setup["target_indexes"];
    arma::uvec targetphi_indexes = estimator_setup["targetphi_indexes"];
    arma::uvec targetpsi_indexes = estimator_setup["targetpsi_indexes"];
    arma::uvec indices = estimator_setup["indices"];

    myestimator->R = R;
    myestimator->lambda = lambda;
    myestimator->phi = phi;
    myestimator->psi = psi;
    myestimator->W = W;
    myestimator->nparam = nparam;
    myestimator->lambda_indexes = lambda_indexes;
    myestimator->phi_indexes = phi_indexes;
    myestimator->psi_indexes = psi_indexes;
    myestimator->target_indexes = target_indexes;
    myestimator->targetphi_indexes = targetphi_indexes;
    myestimator->targetpsi_indexes = targetpsi_indexes;
    myestimator->indices = indices;

    criterion = myestimator;

  } else if(estimator == "cfa_ml") {

    cfa_ml* myestimator = new cfa_ml();

    arma::mat R = estimator_setup["R"];
    arma::mat lambda = estimator_setup["lambda"];
    arma::mat phi = estimator_setup["phi"];
    arma::mat psi = estimator_setup["psi"];
    double logdetR = estimator_setup["logdetR"];
    int nparam = estimator_setup["nparam"];
    int p = estimator_setup["p"];
    arma::uvec lambda_indexes = estimator_setup["lambda_indexes"];
    arma::uvec phi_indexes = estimator_setup["phi_indexes"];
    arma::uvec psi_indexes = estimator_setup["psi_indexes"];
    arma::uvec target_indexes = estimator_setup["target_indexes"];
    arma::uvec targetphi_indexes = estimator_setup["targetphi_indexes"];
    arma::uvec targetpsi_indexes = estimator_setup["targetpsi_indexes"];
    arma::uvec indices = estimator_setup["indices"];

    myestimator->R = R;
    myestimator->logdetR = logdetR;
    myestimator->lambda = lambda;
    myestimator->phi = phi;
    myestimator->psi = psi;
    myestimator->nparam = nparam;
    myestimator->p = p;
    myestimator->lambda_indexes = lambda_indexes;
    myestimator->phi_indexes = phi_indexes;
    myestimator->psi_indexes = psi_indexes;
    myestimator->target_indexes = target_indexes;
    myestimator->targetphi_indexes = targetphi_indexes;
    myestimator->targetpsi_indexes = targetpsi_indexes;
    myestimator->indices = indices;

    criterion = myestimator;

  } else if (estimator == "cf") {

    cf* myestimator = new cf();

    arma::mat lambda = estimator_setup["lambda"];
    bool orth = estimator_setup["orth"];
    double k = estimator_setup["k"];
    arma::uvec indices = estimator_setup["indices"];

    int p = lambda.n_rows;
    int q = lambda.n_cols;
    arma::mat M(p, p, arma::fill::ones);
    M.diag(0).zeros();
    arma::mat N(q, q, arma::fill::ones);
    N.diag(0).zeros();

    myestimator->lambda = lambda;
    myestimator->orth = orth;
    myestimator->indices = indices;
    myestimator->q = q;
    myestimator->Mm = M;
    myestimator->N = N;
    myestimator->k = k;

    criterion = myestimator;

  } else if (estimator == "oblimin") {

    oblimin* myestimator = new oblimin();

    arma::mat lambda = estimator_setup["lambda"];
    bool orth = estimator_setup["orth"];
    double gamma = estimator_setup["gamma"];
    arma::uvec indices = estimator_setup["indices"];

    int p = lambda.n_rows;
    int q = lambda.n_cols;
    arma::mat N(q, q, arma::fill::ones);
    N.diag(0).zeros();
    arma::mat I(p, p, arma::fill::eye), gamma_C(p, p, arma::fill::ones);
    gamma_C *= (gamma/p);
    arma::mat I_gamma_C = (I - gamma_C);

    myestimator->lambda = lambda;
    myestimator->orth = orth;
    myestimator->indices = indices;
    myestimator->q = q;
    myestimator->N = N;
    myestimator->I_gamma_C = I_gamma_C;

    criterion = myestimator;

  } else if (estimator == "geomin") {

    geomin* myestimator = new geomin();

    arma::mat lambda = estimator_setup["lambda"];
    bool orth = estimator_setup["orth"];
    double epsilon = estimator_setup["epsilon"];
    arma::uvec indices = estimator_setup["indices"];

    int p = lambda.n_rows;
    int q = lambda.n_cols;

    myestimator->lambda = lambda;
    myestimator->orth = orth;
    myestimator->indices = indices;
    myestimator->q = q;
    myestimator->epsilon = epsilon;

    criterion = myestimator;

  } else if (estimator == "varimax") {

    varimax* myestimator = new varimax();

    arma::mat lambda = estimator_setup["lambda"];
    bool orth = estimator_setup["orth"];
    arma::uvec indices = estimator_setup["indices"];

    int p = lambda.n_rows;
    int q = lambda.n_cols;

    arma::vec v(p, arma::fill::ones);
    arma::mat I(p, p, arma::fill::eye);
    arma::mat H = I - v * v.t() / (p + 0.0);

    myestimator->lambda = lambda;
    myestimator->orth = orth;
    myestimator->indices = indices;
    myestimator->Hh = H;
    myestimator->q = q;

    criterion = myestimator;

  } else if (estimator == "varimin") {

    varimin* myestimator = new varimin();

    arma::mat lambda = estimator_setup["lambda"];
    bool orth = estimator_setup["orth"];
    arma::uvec indices = estimator_setup["indices"];

    int p = lambda.n_rows;
    int q = lambda.n_cols;

    arma::vec v(p, arma::fill::ones);
    arma::mat I(p, p, arma::fill::eye);
    arma::mat H = I - v * v.t() / (p + 0.0);

    myestimator->lambda = lambda;
    myestimator->orth = orth;
    myestimator->indices = indices;
    myestimator->Hh = H;
    myestimator->q = q;

    criterion = myestimator;

  } else if (estimator == "target") {

    target* myestimator = new target();

    arma::mat lambda = estimator_setup["lambda"];
    bool orth = estimator_setup["orth"];
    arma::uvec indices = estimator_setup["indices"];
    arma::mat Target = estimator_setup["Target"];
    arma::mat Weight = estimator_setup["Weight"];

    arma::mat Weight2 = Weight % Weight;
    int p = lambda.n_rows;
    int q = lambda.n_cols;

    myestimator->lambda = lambda;
    myestimator->orth = orth;
    myestimator->indices = indices;
    myestimator->Target = Target;
    myestimator->Weight = Weight;
    myestimator->Weight2 = Weight2;
    myestimator->q = q;

    criterion = myestimator;

  } else if (estimator == "xtarget") {

    xtarget* myestimator = new xtarget();

    arma::mat lambda = estimator_setup["lambda"];
    bool orth = estimator_setup["orth"];
    arma::uvec indices = estimator_setup["indices"];
    arma::mat Target = estimator_setup["Target"];
    arma::mat Weight = estimator_setup["Weight"];
    arma::mat PhiTarget = estimator_setup["PhiTarget"];
    arma::mat PhiWeight = estimator_setup["PhiWeight"];
    double w = estimator_setup["w"];

    arma::mat Weight2 = Weight % Weight;
    arma::mat PhiWeight2 = PhiWeight % PhiWeight;
    int p = lambda.n_rows;
    int q = lambda.n_cols;

    myestimator->lambda = lambda;
    myestimator->orth = orth;
    myestimator->indices = indices;
    myestimator->Target = Target;
    myestimator->Weight = Weight;
    myestimator->Weight2 = Weight2;
    myestimator->PhiTarget = PhiTarget;
    myestimator->PhiWeight = PhiWeight;
    myestimator->PhiWeight2 = PhiWeight2;
    myestimator->w = w;
    myestimator->q = q;

    criterion = myestimator;

  } else if (estimator == "lclf") {

    lclf* myestimator = new lclf();

    arma::mat lambda = estimator_setup["lambda"];
    bool orth = estimator_setup["orth"];
    double epsilon = estimator_setup["epsilon"];
    arma::uvec indices = estimator_setup["indices"];

    int p = lambda.n_rows;
    int q = lambda.n_cols;

    myestimator->lambda = lambda;
    myestimator->orth = orth;
    myestimator->indices = indices;
    myestimator->q = q;
    myestimator->epsilon = epsilon;

    criterion = myestimator;

  } else if (estimator == "lca_multinomial_softmax") {

    lca_multinomial_softmax* myestimator = new lca_multinomial_softmax();

    arma::mat Y = estimator_setup["Y"];
    int S = estimator_setup["S"];
    int J = estimator_setup["J"];
    arma::uvec K = estimator_setup["K"];
    int nclasses = estimator_setup["nclasses"];
    arma::vec n = estimator_setup["n"];
    arma::vec classes = estimator_setup["classes"];
    // arma::vec logclasses = estimator_setup["logclasses"];
    double constant_class = estimator_setup["constant_class"];
    bool fix_logit = estimator_setup["fix_logit"];
    std::vector<arma::mat> conditionals = estimator_setup["conditionals"];
    // std::vector<arma::mat> logconditionals = estimator_setup["logconditionals"];
    std::vector<arma::vec> constant_cond = estimator_setup["constant_cond"];;
    // std::vector<arma::mat> gconditionals = estimator_setup["gconditionals"];
    arma::uvec indices_classes = estimator_setup["indices_classes"];
    arma::uvec indices_target_classes = estimator_setup["indices_target_classes"];
    std::vector<std::vector<arma::uvec>> indices_conditionals = estimator_setup["indices_conditionals"];
    std::vector<std::vector<arma::uvec>> indices_target_conditionals = estimator_setup["indices_target_conditionals"];
    arma::uvec indices = estimator_setup["indices"];
    // arma::uvec targets = estimator_setup["targets"];
    arma::uvec indices_conditionals2 = estimator_setup["indices_conditionals2"];
    arma::uvec indices_target_conditionals2 = estimator_setup["indices_target_conditionals2"];

    arma::cube loglik(J, nclasses, S, arma::fill::zeros);
    arma::vec logliks(S, arma::fill::zeros);
    arma::vec mid(S, arma::fill::zeros);
    arma::vec loglik_case(S, arma::fill::zeros);
    arma::mat pYXX(S, nclasses, arma::fill::zeros);
    arma::mat latentloglik(S, nclasses, arma::fill::zeros);
    arma::mat jointlogp(S, nclasses, arma::fill::zeros);

    myestimator->Y = Y;
    myestimator->S = S;
    myestimator->J = J;
    myestimator->K = K;
    myestimator->nclasses = nclasses;
    myestimator->n = n;
    myestimator->lclasses = classes;
    myestimator->constant_class = constant_class;
    myestimator->fix_logit = fix_logit;
    myestimator->conditionals = conditionals;
    myestimator->constant_cond = constant_cond;
    myestimator->indices_classes = indices_classes;
    myestimator->indices_target_classes = indices_target_classes;
    myestimator->indices_conditionals = indices_conditionals;
    myestimator->indices_target_conditionals = indices_target_conditionals;
    myestimator->indices = indices;
    myestimator->indices_conditionals2 = indices_conditionals2;
    myestimator->indices_target_conditionals2 = indices_target_conditionals2;

    myestimator->logconditionals = conditionals;
    myestimator->gconditionals = conditionals;
    myestimator->loglik = loglik;
    myestimator->logliks = logliks;
    myestimator->mid = mid;
    myestimator->loglik_case = loglik_case;
    myestimator->pYXX = pYXX;
    myestimator->latentloglik = latentloglik;
    myestimator->jointlogp = jointlogp;
    myestimator->evalE = false;

    criterion = myestimator;

  } else if (estimator == "lca_multinomial") {

    lca_multinomial* myestimator = new lca_multinomial();

    arma::mat Y = estimator_setup["Y"];
    int S = estimator_setup["S"];
    int J = estimator_setup["J"];
    arma::uvec K = estimator_setup["K"];
    int nclasses = estimator_setup["nclasses"];
    arma::vec n = estimator_setup["n"];
    arma::vec classes = estimator_setup["classes"];
    // arma::vec logclasses = estimator_setup["logclasses"];
    double constant_class = estimator_setup["constant_class"];
    bool fix_logit = estimator_setup["fix_logit"];
    std::vector<arma::mat> conditionals = estimator_setup["conditionals"];
    // std::vector<arma::mat> logconditionals = estimator_setup["logconditionals"];
    std::vector<arma::vec> constant_cond = estimator_setup["constant_cond"];;
    // std::vector<arma::mat> gconditionals = estimator_setup["gconditionals"];
    arma::uvec indices_classes = estimator_setup["indices_classes"];
    arma::uvec indices_target_classes = estimator_setup["indices_target_classes"];
    std::vector<std::vector<arma::uvec>> indices_conditionals = estimator_setup["indices_conditionals"];
    std::vector<std::vector<arma::uvec>> indices_target_conditionals = estimator_setup["indices_target_conditionals"];
    arma::uvec indices = estimator_setup["indices"];
    // arma::uvec targets = estimator_setup["targets"];
    arma::uvec indices_conditionals2 = estimator_setup["indices_conditionals2"];
    arma::uvec indices_target_conditionals2 = estimator_setup["indices_target_conditionals2"];

    arma::cube loglik(J, nclasses, S, arma::fill::zeros);
    arma::vec logliks(S, arma::fill::zeros);
    arma::vec mid(S, arma::fill::zeros);
    arma::vec loglik_case(S, arma::fill::zeros);
    arma::mat pYXX(S, nclasses, arma::fill::zeros);
    // arma::mat posterior(S, nclasses, arma::fill::randu);
    // posterior = arma::normalise(posterior, 1);
    arma::mat latentloglik(S, nclasses, arma::fill::zeros);
    arma::mat jointlogp(S, nclasses, arma::fill::zeros);

    myestimator->Y = Y;
    myestimator->S = S;
    myestimator->J = J;
    myestimator->K = K;
    myestimator->nclasses = nclasses;
    myestimator->n = n;
    myestimator->classes = classes;
    myestimator->constant_class = constant_class;
    myestimator->fix_logit = fix_logit;
    myestimator->conditionals = conditionals;
    myestimator->constant_cond = constant_cond;
    myestimator->indices_classes = indices_classes;
    myestimator->indices_target_classes = indices_target_classes;
    myestimator->indices_conditionals = indices_conditionals;
    myestimator->indices_target_conditionals = indices_target_conditionals;
    myestimator->indices = indices;
    myestimator->indices_conditionals2 = indices_conditionals2;
    myestimator->indices_target_conditionals2 = indices_target_conditionals2;

    myestimator->logconditionals = conditionals;
    myestimator->gconditionals = conditionals;
    myestimator->loglik = loglik;
    myestimator->logliks = logliks;
    myestimator->mid = mid;
    myestimator->loglik_case = loglik_case;
    myestimator->pYXX = pYXX;
    // myestimator->posterior = posterior;
    myestimator->latentloglik = latentloglik;
    myestimator->jointlogp = jointlogp;
    myestimator->evalE = true;

    criterion = myestimator;

  } else if (estimator == "lca_gaussian_softmax") {

    lca_gaussian_softmax* myestimator = new lca_gaussian_softmax();

    arma::mat Y = estimator_setup["Y"];
    int S = estimator_setup["S"];
    int J = estimator_setup["J"];
    int nclasses = estimator_setup["nclasses"];
    arma::vec n = estimator_setup["n"];
    arma::vec classes = estimator_setup["classes"];
    double constant_class = estimator_setup["constant_class"];
    bool fix_logit = estimator_setup["fix_logit"];
    std::vector<arma::mat> conditionals = estimator_setup["conditionals"];
    arma::uvec indices_classes = estimator_setup["indices_classes"];
    arma::uvec indices_target_classes = estimator_setup["indices_target_classes"];
    std::vector<std::vector<arma::uvec>> indices_conditionals = estimator_setup["indices_conditionals"];
    std::vector<std::vector<arma::uvec>> indices_target_conditionals = estimator_setup["indices_target_conditionals"];
    arma::uvec indices = estimator_setup["indices"];
    arma::uvec indices_conditionals2 = estimator_setup["indices_conditionals2"];
    arma::uvec indices_target_conditionals2 = estimator_setup["indices_target_conditionals2"];

    arma::cube loglik(J, nclasses, S, arma::fill::zeros);
    arma::vec logliks(S, arma::fill::zeros);
    arma::vec mid(S, arma::fill::zeros);
    arma::vec loglik_case(S, arma::fill::zeros);
    arma::mat pYXX(S, nclasses, arma::fill::zeros);
    // arma::mat posterior(S, nclasses, arma::fill::randu);
    // posterior = arma::normalise(posterior, 1);
    arma::mat latentloglik(S, nclasses, arma::fill::zeros);
    arma::mat jointlogp(S, nclasses, arma::fill::zeros);

    myestimator->Y = Y;
    myestimator->S = S;
    myestimator->J = J;
    myestimator->nclasses = nclasses;
    myestimator->n = n;
    myestimator->lclasses = classes;
    myestimator->constant_class = constant_class;
    myestimator->fix_logit = fix_logit;
    myestimator->conditionals = conditionals;
    myestimator->indices_classes = indices_classes;
    myestimator->indices_target_classes = indices_target_classes;
    myestimator->indices_conditionals = indices_conditionals;
    myestimator->indices_target_conditionals = indices_target_conditionals;
    myestimator->indices = indices;
    myestimator->indices_conditionals2 = indices_conditionals2;
    myestimator->indices_target_conditionals2 = indices_target_conditionals2;

    myestimator->gconditionals = conditionals;
    myestimator->loglik = loglik;
    myestimator->logliks = logliks;
    myestimator->mid = mid;
    myestimator->loglik_case = loglik_case;
    myestimator->pYXX = pYXX;
    // myestimator->posterior = posterior;
    myestimator->latentloglik = latentloglik;
    myestimator->jointlogp = jointlogp;
    myestimator->evalE = false;

    criterion = myestimator;

  } else if (estimator == "lca_gaussian") {

    lca_gaussian* myestimator = new lca_gaussian();

    arma::mat Y = estimator_setup["Y"];
    int S = estimator_setup["S"];
    int J = estimator_setup["J"];
    int nclasses = estimator_setup["nclasses"];
    arma::vec n = estimator_setup["n"];
    arma::vec classes = estimator_setup["classes"];
    double constant_class = estimator_setup["constant_class"];
    bool fix_logit = estimator_setup["fix_logit"];
    std::vector<arma::mat> conditionals = estimator_setup["conditionals"];
    arma::uvec indices_classes = estimator_setup["indices_classes"];
    arma::uvec indices_target_classes = estimator_setup["indices_target_classes"];
    std::vector<std::vector<arma::uvec>> indices_conditionals = estimator_setup["indices_conditionals"];
    std::vector<std::vector<arma::uvec>> indices_target_conditionals = estimator_setup["indices_target_conditionals"];
    arma::uvec indices = estimator_setup["indices"];
    arma::uvec indices_conditionals2 = estimator_setup["indices_conditionals2"];
    arma::uvec indices_target_conditionals2 = estimator_setup["indices_target_conditionals2"];

    arma::cube loglik(J, nclasses, S, arma::fill::zeros);
    arma::vec logliks(S, arma::fill::zeros);
    arma::vec mid(S, arma::fill::zeros);
    arma::vec loglik_case(S, arma::fill::zeros);
    arma::mat pYXX(S, nclasses, arma::fill::zeros);
    // arma::mat posterior(S, nclasses, arma::fill::randu);
    // posterior = arma::normalise(posterior, 1);
    arma::mat latentloglik(S, nclasses, arma::fill::zeros);
    arma::mat jointlogp(S, nclasses, arma::fill::zeros);

    myestimator->Y = Y;
    myestimator->S = S;
    myestimator->J = J;
    myestimator->nclasses = nclasses;
    myestimator->n = n;
    myestimator->classes = classes;
    myestimator->constant_class = constant_class;
    myestimator->fix_logit = fix_logit;
    myestimator->conditionals = conditionals;
    myestimator->indices_classes = indices_classes;
    myestimator->indices_target_classes = indices_target_classes;
    myestimator->indices_conditionals = indices_conditionals;
    myestimator->indices_target_conditionals = indices_target_conditionals;
    myestimator->indices = indices;
    myestimator->indices_conditionals2 = indices_conditionals2;
    myestimator->indices_target_conditionals2 = indices_target_conditionals2;

    myestimator->gconditionals = conditionals;
    myestimator->loglik = loglik;
    myestimator->logliks = logliks;
    myestimator->mid = mid;
    myestimator->loglik_case = loglik_case;
    myestimator->pYXX = pYXX;
    // myestimator->posterior = posterior;
    myestimator->latentloglik = latentloglik;
    myestimator->jointlogp = jointlogp;
    myestimator->evalE = true;

    criterion = myestimator;

  } else if (estimator == "latentloglik_combination") {

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
    myestimator->jointlogp = jointlogp;
    myestimator->evalE = false;

    criterion = myestimator;

  } else if (estimator == "latentloglik_combination_softmax") {

    latentloglik_combination_softmax* myestimator = new latentloglik_combination_softmax();

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
    myestimator->jointlogp = jointlogp;
    myestimator->evalE = false;

    criterion = myestimator;

  } else if (estimator == "lreg") {

    lreg* myestimator = new lreg();

    arma::uvec indices = estimator_setup["indices"];
    arma::mat Y = estimator_setup["Y"];
    arma::mat predictors = estimator_setup["predictors"];

    myestimator->indices = indices;
    myestimator->Y = Y;
    myestimator->predictors = predictors;

    criterion = myestimator;

  } else {

    Rcpp::stop("Unkown estimator");

  }

  return criterion;

}

// Product Estimator

class product_estimator {

public:

  void param(arguments_optim& x, std::vector<estimators*>& xestimators) {

    // double eps = std::numeric_limits<double>::epsilon();
    // x.parameters.elem(arma::find(arma::abs(x.parameters) < eps)).fill(eps);

    x.latentloglik.zeros();
    // x.latentpars = xestimators[1]->latentpars;

    for(int i=0; i < x.nestimators; ++i) {

      arma::uvec indices = xestimators[i]->indices;
      xestimators[i]->parameters = x.parameters.elem(indices);
      // xestimators[i]->latentpars = x.latentpars;
      xestimators[i]->latentloglik = x.latentloglik;

      xestimators[i]->param();

      // arma::vec v = xestimators[i]->latentpars;
      // for (arma::uword j = 0; j < v.n_elem; ++j) {
      //   Rprintf("%g \n", v[j]);
      // }

      x.latentloglik += xestimators[i]->latentloglik;
      x.latentpars = xestimators[i]->latentpars; // ADD INDICES

    }

    // Pick the final value of latentloglik:
    // int last = x.nestimators - 1L;
    // x.latentloglik = xestimators[last]->latentloglik;

  }

  void F(arguments_optim& x, std::vector<estimators*>& xestimators) {

    x.f = 0;

    for(int i=0; i < x.nestimators; ++i) {

      xestimators[i]->latentloglik = x.latentloglik;
      xestimators[i]->latentpars = x.latentpars; // ADD INDICES

      xestimators[i]->F();
      x.f += xestimators[i]->f;

    }

  }

  void G(arguments_optim& x, std::vector<estimators*>& xestimators) {

    x.g.set_size(x.parameters.n_elem); x.g.zeros();

    for(int i=0; i < x.nestimators; ++i) {

      xestimators[i]->latentloglik = x.latentloglik;
      xestimators[i]->latentpars = x.latentpars; // ADD INDICES

      xestimators[i]->G();

      arma::uvec indices = xestimators[i]->indices;

      // Rprintf("Indices:\n");
      // for (arma::uword i = 0; i < indices.n_elem; ++i) {
      //   Rprintf("%u ", indices[i]);
      // }
      //
      // for (arma::uword j = 0; j < xestimators[i]->g.n_elem; ++j) {
      //   Rprintf("%g ", xestimators[i]->g[j]);
      // }
      // Rprintf("\n\n");

      x.g.elem(indices) += xestimators[i]->g;

      // for (arma::uword j = 0; j < x.g.n_elem; ++j) {
      //   Rprintf("%g \n", x.g[j]);
      // }

    }

  }

  void dG(arguments_optim& x, std::vector<estimators*>& xestimators) {

    x.dg.set_size(x.parameters.n_elem); x.dg.zeros();

    for(int i=0; i < x.nestimators; ++i) {

      arma::uvec indices = xestimators[i]->indices;
      xestimators[i]->dparameters = x.dparameters.elem(indices);
      xestimators[i]->g = x.g.elem(indices); // Maybe delete PROBLEMS
      xestimators[i]->dG();
      x.dg.elem(indices) += xestimators[i]->dg;

    }

  }

  void H(arguments_optim& x, std::vector<estimators*>& xestimators) {

    int nhessian = x.parameters.n_elem;
    x.hessian.set_size(nhessian, nhessian); x.hessian.zeros();

    for(int i=0; i < x.nestimators; ++i) {

      xestimators[i]->H();
      arma::uvec indices = xestimators[i]->indices;
      x.hessian(indices, indices) += xestimators[i]->hessian;

    }

  }

  void E(arguments_optim& x, std::vector<estimators*>& xestimators) {

    x.latentloglik.zeros();

    for(int i=0; i < x.nestimators; ++i) {

      if(xestimators[i]->evalE) {

        xestimators[i]->posterior = x.posterior;
        arma::uvec indices = xestimators[i]->indices;
        xestimators[i]->E();
        x.parameters.elem(indices) = xestimators[i]->parameters;
        x.latentloglik += xestimators[i]->latentloglik;

      }

    }

  }

  void M(arguments_optim& x, std::vector<estimators*>& xestimators) {

    // x.loglik = 0.00;
    //
    // for(int i=0; i < x.nestimators; ++i) {
    //
    //   xestimators[i]->latentloglik = x.latentloglik;
    //
    //   xestimators[i]->M();
    //   x.loglik += xestimators[i]->loglik;
    //   x.posterior = xestimators[i]->posterior;
    //
    // }
    //
    // x.loglik /= x.nestimators;

    // For each response pattern, compute its joint and marginal probabilities:
    int S = x.latentloglik.n_rows;
    arma::vec mid(S);
    x.latentpars = xestimators[0]->latentpars;
    x.n = xestimators[0]->n;
    // for(int s=0; s < S; ++s) {
    //   x.posterior.row(s) = arma::trunc_exp(x.latentloglik.row(s) + x.latentpars.t()); // P(data |X = c) P(X = c)
    //   mid[s] = arma::accu(x.posterior.row(s)); // P(data)
    // }

    x.posterior = x.latentloglik;
    x.posterior.each_row() += x.latentpars.t();
    x.posterior = arma::trunc_exp(x.posterior);
    mid = arma::sum(x.posterior, 1);

    x.posterior.each_col() /= mid; // P(X = c | data) = P(data | X = c) P(X = c) / P(data)

    arma::vec logliks = x.n % arma::trunc_log(mid);
    x.loglik = -arma::accu(logliks);

  }

  void outcomes(arguments_optim& x, std::vector<estimators*>& xestimators) {

    x.modhessian.resize(x.nestimators);
    x.dparam_dS.resize(x.nestimators);
    x.doubles.resize(x.nestimators);
    x.vectors.resize(x.nestimators);
    x.matrices.resize(x.nestimators);
    x.list_matrices.resize(x.nestimators);

    for(int i=0; i < x.nestimators; ++i) {

      arma::uvec indices = xestimators[i]->indices;
      xestimators[i]->parameters = x.parameters.elem(indices);
      xestimators[i]->outcomes();
      x.modhessian[i] = xestimators[i]->modhessian;
      x.dparam_dS[i] = xestimators[i]->dparam_dS;
      x.doubles[i] = xestimators[i]->doubles;
      x.vectors[i] = xestimators[i]->vectors;
      x.matrices[i] = xestimators[i]->matrices;
      x.list_matrices[i] = xestimators[i]->list_matrices;

    }

  }

};

