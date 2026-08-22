/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 22/08/2026
 */

/*
 * Model-implied ordinal threshold transformation:
 *
 * tauhat = (kappas - meanshat)/sqrt(diag(Shat))
 *
 * Since each observed variable may have a different number of thresholds,
 * threshold_items maps every element of kappas to its observed variable.
 * threshold_items must use zero-based indices.
 */

class tau_param: public transformations {

public:

  int p, ntau;
  arma::uvec indices_kappas, indices_meanshat, indices_Shat,
    threshold_items;
  arma::vec kappas, meanshat, tauhat, dkappas, dmeanshat, dtauhat,
    variances, sqrt_variances, inv_sqrt_variances, grad_out;
  arma::mat Shat, dShat;

  void transform(arguments_optim& x) {

    kappas = x.transparameters.elem(indices_kappas);
    meanshat = x.transparameters.elem(indices_meanshat);
    Shat = arma::reshape(x.transparameters.elem(indices_Shat), p, p);

    variances = Shat.diag();

    if(!variances.is_finite() || arma::any(variances <= 0.0)) {
      Rcpp::stop("tau_param requires strictly positive finite diagonal elements of Shat.");
    }

    sqrt_variances = arma::sqrt(variances);
    inv_sqrt_variances = 1.0/sqrt_variances;

    tauhat.set_size(ntau);

    for(int k = 0; k < ntau; ++k) {

      const arma::uword j = threshold_items[k];
      tauhat[k] =
        (kappas[k]-meanshat[j])*inv_sqrt_variances[j];

    }

    x.transparameters.elem(indices_out) = tauhat;

  }

  void update_grad(arguments_optim& x) {

    grad_out = x.grad.elem(indices_out);

    arma::vec grad_kappas(ntau, arma::fill::zeros);
    arma::vec grad_meanshat(p, arma::fill::zeros);
    arma::mat grad_Shat(p, p, arma::fill::zeros);

    for(int k = 0; k < ntau; ++k) {

      const arma::uword j = threshold_items[k];
      const double a = kappas[k]-meanshat[j];
      const double inv_sd = inv_sqrt_variances[j];
      const double inv_var_sd = inv_sd/variances[j];
      const double g = grad_out[k];

      grad_kappas[k] += g*inv_sd;
      grad_meanshat[j] -= g*inv_sd;
      grad_Shat(j, j) -= 0.5*g*a*inv_var_sd;

    }

    x.grad.elem(indices_kappas) += grad_kappas;
    x.grad.elem(indices_meanshat) += grad_meanshat;
    x.grad.elem(indices_Shat) += arma::vectorise(grad_Shat);

  }

  void dtransform(arguments_optim& x) {

    dkappas = x.dtransparameters.elem(indices_kappas);
    dmeanshat = x.dtransparameters.elem(indices_meanshat);
    dShat = arma::reshape(x.dtransparameters.elem(indices_Shat), p, p);

    dtauhat.set_size(ntau);

    for(int k = 0; k < ntau; ++k) {

      const arma::uword j = threshold_items[k];
      const double a = kappas[k]-meanshat[j];
      const double da = dkappas[k]-dmeanshat[j];
      const double v = variances[j];
      const double dv = dShat(j, j);
      const double inv_sd = inv_sqrt_variances[j];
      const double inv_var_sd = inv_sd/v;

      dtauhat[k] =
        da*inv_sd-
        0.5*a*inv_var_sd*dv;

    }

    x.dtransparameters.elem(indices_out) = dtauhat;

  }

  void update_dgrad(arguments_optim& x) {

    arma::vec dgrad_out = x.dgrad.elem(indices_out);

    arma::vec dgrad_kappas(ntau, arma::fill::zeros);
    arma::vec dgrad_meanshat(p, arma::fill::zeros);
    arma::mat dgrad_Shat(p, p, arma::fill::zeros);

    for(int k = 0; k < ntau; ++k) {

      const arma::uword j = threshold_items[k];
      const double a = kappas[k]-meanshat[j];
      const double da = dkappas[k]-dmeanshat[j];
      const double v = variances[j];
      const double dv = dShat(j, j);
      const double inv_sd = inv_sqrt_variances[j];
      const double inv_v32 = inv_sd/v;
      const double inv_v52 = inv_v32/v;
      const double g = grad_out[k];
      const double dg = dgrad_out[k];
      const double d_inv_sd = -0.5*inv_v32*dv;
      const double h = -0.5*a*inv_v32;
      const double dh =
        -0.5*da*inv_v32+
        0.75*a*inv_v52*dv;

      dgrad_kappas[k] +=
        dg*inv_sd+
        g*d_inv_sd;

      dgrad_meanshat[j] -=
        dg*inv_sd+
        g*d_inv_sd;

      dgrad_Shat(j, j) +=
        dg*h+
        g*dh;

    }

    x.dgrad.elem(indices_kappas) += dgrad_kappas;
    x.dgrad.elem(indices_meanshat) += dgrad_meanshat;
    x.dgrad.elem(indices_Shat) += arma::vectorise(dgrad_Shat);

  }

  void jacobian(arguments_optim& x) {

    kappas = x.transparameters.elem(indices_kappas);
    meanshat = x.transparameters.elem(indices_meanshat);
    Shat = arma::reshape(x.transparameters.elem(indices_Shat), p, p);

    variances = Shat.diag();

    if(!variances.is_finite() || arma::any(variances <= 0.0)) {
      Rcpp::stop("tau_param requires strictly positive finite diagonal elements of Shat.");
    }

    sqrt_variances = arma::sqrt(variances);
    inv_sqrt_variances = 1.0/sqrt_variances;

    const arma::uword ninputs =
      indices_kappas.n_elem+
      indices_meanshat.n_elem+
      indices_Shat.n_elem;

    jacob.zeros(ntau, ninputs);

    const arma::uword offset_kappas = 0L;
    const arma::uword offset_meanshat = indices_kappas.n_elem;
    const arma::uword offset_Shat =
      indices_kappas.n_elem+indices_meanshat.n_elem;

    for(int k = 0; k < ntau; ++k) {

      const arma::uword j = threshold_items[k];
      const double a = kappas[k]-meanshat[j];
      const double inv_sd = inv_sqrt_variances[j];
      const double inv_var_sd = inv_sd/variances[j];
      const arma::uword jj = j+j*p;

      jacob(k, offset_kappas+k) = inv_sd;
      jacob(k, offset_meanshat+j) = -inv_sd;
      jacob(k, offset_Shat+jj) = -0.5*a*inv_var_sd;

    }

  }

  void outcomes(arguments_optim& x) {

    (void)x;

    vectors.resize(1);
    vectors[0] = tauhat;
    names_vectors.resize(1);
    names_vectors[0] = "tauhat";

    matrices.resize(1);
    matrices[0] = jacob;
    names_matrices.resize(1);
    names_matrices[0] = "jacobian";

  }

};

tau_param* choose_tau_param(const Rcpp::List& trans_setup) {

  tau_param* mytrans = new tau_param();

  std::vector<arma::uvec> indices_in = trans_setup["indices_in"];
  std::vector<arma::uvec> indices_out = trans_setup["indices_out"];
  arma::uvec threshold_items = trans_setup["threshold_items"];
  int p = trans_setup["p"];

  if(p < 1) {
    Rcpp::stop("tau_param requires a positive p dimension.");
  }

  if(indices_in.size() != 3L || indices_out.size() != 1L) {
    Rcpp::stop("tau_param requires kappas, meanshat, and Shat inputs and one tauhat output.");
  }

  arma::uvec indices_kappas = indices_in[0];
  arma::uvec indices_meanshat = indices_in[1];
  arma::uvec indices_Shat = indices_in[2];
  arma::uvec indices_tauhat = indices_out[0];

  const arma::uword ntau = indices_kappas.n_elem;

  if(indices_meanshat.n_elem != static_cast<arma::uword>(p)) {
    Rcpp::stop("tau_param requires one model-implied mean per observed variable.");
  }

  if(indices_Shat.n_elem != static_cast<arma::uword>(p*p)) {
    Rcpp::stop("tau_param requires Shat to be a p by p matrix.");
  }

  if(indices_tauhat.n_elem != ntau) {
    Rcpp::stop("tau_param requires tauhat and kappas to have the same length.");
  }

  if(threshold_items.n_elem != ntau) {
    Rcpp::stop("threshold_items must contain one variable index per threshold.");
  }

  if(ntau > 0L && arma::any(threshold_items >= static_cast<arma::uword>(p))) {
    Rcpp::stop("threshold_items contains an index outside the range of observed variables.");
  }

  mytrans->indices_kappas = indices_kappas;
  mytrans->indices_meanshat = indices_meanshat;
  mytrans->indices_Shat = indices_Shat;
  mytrans->indices_in = arma::join_cols(
    indices_kappas,
    arma::join_cols(indices_meanshat, indices_Shat)
  );
  mytrans->indices_out = indices_tauhat;
  mytrans->threshold_items = threshold_items;
  mytrans->p = p;
  mytrans->ntau = static_cast<int>(ntau);

  return mytrans;

}
