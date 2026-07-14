/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 14/07/2026
 */

class mvnormal: public transformations {

public:

  int S, J, I, n_in, n_out;
  arma::uvec indices_mu, indices_sigma, indices_loglik, indices_in, indices_out;
  arma::mat y, mu, dmu, loglik, jacob;
  arma::cube Sigma, dSigma;

  void transform(arguments_optim& x) {

    mu = arma::reshape(x.transparameters.elem(indices_mu), J, I);
    loglik = arma::reshape(x.transparameters.elem(indices_loglik), S, I);

    Sigma.set_size(J, J, I);
    for (int i = 0; i < I; ++i) {
      Sigma.slice(i) = arma::reshape(
        x.transparameters.elem(indices_sigma.subvec(i * J * J, (i + 1) * J * J - 1)),
        J, J
      );
      // For Sigma to be positive-definite:
      if(!Sigma.slice(i).is_sympd()) {
        arma::vec eigval;
        arma::mat eigvec;
        eig_sym(eigval, eigvec, Sigma.slice(i));
        arma::vec d = arma::clamp(eigval, 0.00001, eigval.max());
        Sigma.slice(i) = eigvec * arma::diagmat(d) * eigvec.t();
      }
    }

    for (int i = 0; i < I; ++i) {
      for (int s = 0; s < S; ++s) {

        arma::vec ys = y.row(s).t();
        arma::uvec obs = arma::find_finite(ys);
        const arma::uword m = obs.n_elem;
        if (m == 0) continue;

        arma::vec mui = mu.col(i);
        arma::vec d = ys.elem(obs) - mui.elem(obs);
        arma::mat Sig = Sigma.slice(i).submat(obs, obs);
        arma::mat K = arma::inv_sympd(Sig);

        double logdetSig = arma::log_det_sympd(Sig);
        double quad = arma::as_scalar(d.t() * K * d);

        loglik(s, i) += -0.5 * (static_cast<double>(m) * std::log(2.0 * M_PI) + logdetSig + quad);
      }
    }

    x.transparameters.elem(indices_out) = arma::vectorise(loglik);

  }

  void update_grad(arguments_optim& x) {

    arma::mat df_dmu(J, I, arma::fill::zeros);
    arma::cube df_dSigma(J, J, I, arma::fill::zeros);
    arma::vec grad_out = x.grad.elem(indices_out);

    int k = 0L;
    for (int i = 0; i < I; ++i) {
      for (int s = 0; s < S; ++s, ++k) {

        arma::vec ys = y.row(s).t();
        arma::uvec obs = arma::find_finite(ys);
        const arma::uword m = obs.n_elem;
        if (m == 0) continue;

        arma::vec mui = mu.col(i);
        arma::vec d = ys.elem(obs) - mui.elem(obs);
        arma::mat Sig = Sigma.slice(i).submat(obs, obs);
        arma::mat K = arma::inv_sympd(Sig);

        arma::vec gmu = K * d;
        arma::mat gSigma = 0.5 * (gmu * gmu.t() - K);

        const double gk = grad_out(k);

        for (arma::uword a = 0; a < m; ++a) {
          df_dmu(obs(a), i) += gk * gmu(a);
        }

        for (arma::uword a = 0; a < m; ++a) {
          for (arma::uword b = 0; b < m; ++b) {
            df_dSigma(obs(a), obs(b), i) += gk * gSigma(a, b);
          }
        }
      }
    }

    arma::vec grad_sigma(indices_sigma.n_elem, arma::fill::zeros);
    for (int i = 0; i < I; ++i) {
      grad_sigma.subvec(i * J * J, (i + 1) * J * J - 1) =
        arma::vectorise(df_dSigma.slice(i));
    }

    x.grad.elem(indices_mu) += arma::vectorise(df_dmu);
    x.grad.elem(indices_sigma) += grad_sigma;

  }

  void dtransform(arguments_optim& x) {

    dmu = arma::reshape(x.dtransparameters.elem(indices_mu), J, I);

    dSigma.set_size(J, J, I);
    for (int i = 0; i < I; ++i) {
      dSigma.slice(i) = arma::reshape(
        x.dtransparameters.elem(indices_sigma.subvec(i * J * J, (i + 1) * J * J - 1)),
        J, J
      );
    }

    arma::mat dloglik = arma::reshape(x.dtransparameters.elem(indices_loglik), S, I);

    for (int i = 0; i < I; ++i) {
      for (int s = 0; s < S; ++s) {

        arma::vec ys = y.row(s).t();
        arma::uvec obs = arma::find_finite(ys);
        const arma::uword m = obs.n_elem;
        if (m == 0) continue;

        arma::vec mui = mu.col(i);
        arma::vec dmui = dmu.col(i);
        arma::vec d = ys.elem(obs) - mui.elem(obs);
        arma::vec dmu_obs = dmui.elem(obs);
        arma::mat Sig = Sigma.slice(i).submat(obs, obs);
        arma::mat dSig = dSigma.slice(i).submat(obs, obs);
        arma::mat K = arma::inv_sympd(Sig);

        arma::vec gmu = K * d;
        arma::mat gSigma = 0.5 * (gmu * gmu.t() - K);

        dloglik(s, i) +=
          arma::dot(gmu, dmu_obs) +
          arma::accu(gSigma % dSig);
      }
    }

    x.dtransparameters.elem(indices_out) = arma::vectorise(dloglik);

  }

  void update_dgrad(arguments_optim& x) {

    arma::mat ddf_dmu(J, I, arma::fill::zeros);
    arma::cube ddf_dSigma(J, J, I, arma::fill::zeros);

    arma::vec grad_out = x.grad.elem(indices_out);
    arma::vec dgrad_out = x.dgrad.elem(indices_out);

    int k = 0L;
    for (int i = 0; i < I; ++i) {
      for (int s = 0; s < S; ++s, ++k) {

        arma::vec ys = y.row(s).t();
        arma::uvec obs = arma::find_finite(ys);
        const arma::uword m = obs.n_elem;
        if (m == 0) continue;

        arma::vec mui = mu.col(i);
        arma::vec dmui = dmu.col(i);
        arma::vec d = ys.elem(obs) - mui.elem(obs);
        arma::vec dmu_obs = dmui.elem(obs);
        arma::mat Sig = Sigma.slice(i).submat(obs, obs);
        arma::mat dSig = dSigma.slice(i).submat(obs, obs);
        arma::mat K = arma::inv_sympd(Sig);
        arma::mat dK = -K * dSig * K;

        arma::vec gmu = K * d;
        arma::mat gSigma = 0.5 * (gmu * gmu.t() - K);

        arma::vec dgmu = dK * d - K * dmu_obs;
        arma::mat dgSigma = 0.5 * (dgmu * gmu.t() + gmu * dgmu.t() - dK);

        const double gk = grad_out(k);
        const double dgk = dgrad_out(k);

        for (arma::uword a = 0; a < m; ++a) {
          ddf_dmu(obs(a), i) += dgk * gmu(a) + gk * dgmu(a);
        }

        for (arma::uword a = 0; a < m; ++a) {
          for (arma::uword b = 0; b < m; ++b) {
            ddf_dSigma(obs(a), obs(b), i) += dgk * gSigma(a, b) + gk * dgSigma(a, b);
          }
        }
      }
    }

    arma::vec dgrad_sigma(indices_sigma.n_elem, arma::fill::zeros);
    for (int i = 0; i < I; ++i) {
      dgrad_sigma.subvec(i * J * J, (i + 1) * J * J - 1) =
        arma::vectorise(ddf_dSigma.slice(i));
    }

    x.dgrad.elem(indices_mu) += arma::vectorise(ddf_dmu);
    x.dgrad.elem(indices_sigma) += dgrad_sigma;

  }

  void jacobian(arguments_optim& x) {

    jacob.set_size(n_out, n_in);
    jacob.zeros();

    const int offset_mu = 0;
    const int offset_sigma = J * I;
    const int offset_loglik = J * I + I * J * J;

    int k = 0L;
    for (int i = 0; i < I; ++i) {
      const int mu_offset = offset_mu + i * J;
      const int sigma_offset = offset_sigma + i * J * J;

      for (int s = 0; s < S; ++s, ++k) {

        jacob(k, offset_loglik + k) = 1.0;

        arma::vec ys = y.row(s).t();
        arma::uvec obs = arma::find_finite(ys);
        const arma::uword m = obs.n_elem;
        if (m == 0) continue;

        arma::vec mui = mu.col(i);
        arma::vec d = ys.elem(obs) - mui.elem(obs);
        arma::mat Sig = Sigma.slice(i).submat(obs, obs);
        arma::mat K = arma::inv_sympd(Sig);

        arma::vec gmu = K * d;
        arma::mat gSigma = 0.5 * (gmu * gmu.t() - K);

        for (arma::uword a = 0; a < m; ++a) {
          jacob(k, mu_offset + obs(a)) = gmu(a);
        }

        for (arma::uword a = 0; a < m; ++a) {
          for (arma::uword b = 0; b < m; ++b) {
            jacob(k, sigma_offset + obs(a) + obs(b) * J) = gSigma(a, b);
          }
        }
      }
    }

  }

  void update_vcov(arguments_optim& x) {

    indices_in = arma::join_cols(
      indices_mu,
      arma::join_cols(indices_sigma, indices_loglik)
    );

    x.vcov(indices_out, indices_out) =
      jacob * x.vcov(indices_in, indices_in) * jacob.t();

  }

  void dconstraints(arguments_optim& x) {

    constraints = false;

  }

  void outcomes(arguments_optim& x) {

    matrices.resize(1);
    matrices[0] = jacob;
    names_matrices.resize(1);
    names_matrices[0] = "jacobian";

  }

};

mvnormal* choose_mvnormal(const Rcpp::List& trans_setup) {

  mvnormal* mytrans = new mvnormal();

  std::vector<arma::uvec> indices_in = trans_setup["indices_in"];
  std::vector<arma::uvec> indices_out = trans_setup["indices_out"];
  arma::mat y = trans_setup["y"];
  int S = trans_setup["S"];
  int J = trans_setup["J"];
  int I = trans_setup["I"];

  arma::uvec indices_mu = indices_in[0];
  arma::uvec indices_sigma = indices_in[1];
  arma::uvec indices_loglik = indices_in[2];

  int n_in =
    indices_mu.n_elem +
    indices_sigma.n_elem +
    indices_loglik.n_elem;

  int n_out = indices_out[0].n_elem;

  arma::mat mu(J, I, arma::fill::zeros);
  arma::mat dmu(J, I, arma::fill::zeros);
  arma::mat loglik(S, I, arma::fill::zeros);
  arma::cube Sigma(J, J, I, arma::fill::zeros);
  arma::cube dSigma(J, J, I, arma::fill::zeros);

  mytrans->indices_mu = indices_mu;
  mytrans->indices_sigma = indices_sigma;
  mytrans->indices_loglik = indices_loglik;
  mytrans->indices_out = indices_out[0];
  mytrans->n_in = n_in;
  mytrans->n_out = n_out;
  mytrans->y = y;
  mytrans->S = S;
  mytrans->J = J;
  mytrans->I = I;
  mytrans->mu = mu;
  mytrans->dmu = dmu;
  mytrans->loglik = loglik;
  mytrans->Sigma = Sigma;
  mytrans->dSigma = dSigma;

  return mytrans;

}
