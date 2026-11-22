/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 11/08/2026
 */

#include <map>

class mvnormal: public transformations {

public:

  int S, J, I, n_in, n_out;
  arma::uvec indices_mu, indices_sigma, indices_loglik;
  // Each pattern holds cases with the same observed variables (one case per column).
  struct observation_pattern {
    arma::uvec obs, rows;
    arma::mat data;
    std::vector<arma::mat> precision, mean_score;
  };
  std::vector<observation_pattern> patterns;
  struct case_index {
    arma::uword pattern, column;
  };
  std::vector<case_index> cases;
  arma::vec last_mu, last_sigma, log_density;
  arma::mat dmu;
  arma::cube dSigma;

  void initialize_patterns(const arma::mat& y) {
    std::map<std::vector<arma::uword>, arma::uword> indices;
    std::vector<std::vector<arma::uword>> rows;
    for(int s = 0; s < S; ++s) {
      arma::uvec obs = arma::find_finite(y.row(s));
      if(obs.is_empty()) continue; // A fully missing row adds zero log density.
      std::vector<arma::uword> key(obs.begin(), obs.end());
      auto entry = indices.emplace(key, patterns.size());
      if(entry.second) {
        patterns.emplace_back();
        patterns.back().obs = obs;
        rows.emplace_back();
      }
      arma::uword p = entry.first->second;
      cases.push_back({p, rows[p].size()});
      rows[p].push_back(s);
    }
    for(arma::uword p = 0; p < patterns.size(); ++p) {
      auto& block = patterns[p];
      block.rows = arma::uvec(rows[p]);
      block.data = y.submat(block.rows, block.obs).t();
      block.precision.resize(I);
      block.mean_score.resize(I);
    }
    log_density.zeros(S * I);
  }

  void transform(arguments_optim& x) {
    arma::vec means = x.transparameters.elem(indices_mu);
    arma::vec covariances = x.transparameters.elem(indices_sigma);
    // In a structural fit the measurement parameters can stay fixed across iterations.
    if(last_mu.is_empty() || !arma::all(means == last_mu) ||
       !arma::all(covariances == last_sigma)) {
      last_mu.reset(); // A failed trial must not leave a valid cache key.
      for(int i = 0; i < I; ++i) {
        arma::vec mu = means.subvec(i*J, (i+1)*J-1);
        arma::mat Sigma = arma::reshape(covariances.subvec(i*J*J, (i+1)*J*J-1), J, J);
        if(!Sigma.is_sympd()) {
          arma::vec eigval;
          arma::mat eigvec;
          arma::eig_sym(eigval, eigvec, Sigma);
          arma::vec d = arma::clamp(eigval, 0.00001, eigval.max());
          Sigma = eigvec * arma::diagmat(d) * eigvec.t();
        }
        for(auto& block : patterns) {
          arma::mat Sig = Sigma.submat(block.obs, block.obs);
          arma::mat& K = block.precision[i];
          K = arma::inv_sympd(Sig);
          arma::mat residual = block.data.each_col() - mu.elem(block.obs);
          arma::mat& score = block.mean_score[i];
          score.set_size(residual.n_rows, residual.n_cols);
          double constant = block.obs.n_elem * std::log(2.0 * M_PI) + arma::log_det_sympd(Sig);
          for(arma::uword c = 0; c < residual.n_cols; ++c) {
            arma::vec d(residual.colptr(c), residual.n_rows, false, true);
            score.col(c) = K*d;
            log_density(block.rows(c) + i*S) = -0.5 * (constant + arma::as_scalar(d.t()*K*d));
          }
        }
      }
      last_mu = std::move(means);
      last_sigma = std::move(covariances);
    }
    x.transparameters.elem(indices_out) = x.transparameters.elem(indices_loglik) + log_density;
  }

  void update_grad(arguments_optim& x) {
    for(int i = 0; i < I; ++i) {
      arma::vec grad_mu(J, arma::fill::zeros);
      arma::mat grad_sigma(J, J, arma::fill::zeros);
      // Accumulate in input order: rounding can affect which tied random start wins.
      for(const auto& item : cases) {
        const auto& block = patterns[item.pattern];
        const double* score = block.mean_score[i].colptr(item.column);
        const double w = x.grad(indices_out(block.rows(item.column) + i*S));
        for(arma::uword b = 0; b < block.obs.n_elem; ++b) {
          grad_mu(block.obs(b)) += w*score[b];
          for(arma::uword a = 0; a < block.obs.n_elem; ++a) {
            double cov_score = 0.5 * (score[a]*score[b] - block.precision[i](a,b));
            grad_sigma(block.obs(a), block.obs(b)) += w*cov_score;
          }
        }
      }
      x.grad.elem(indices_mu.subvec(i*J, (i+1)*J-1)) += grad_mu;
      x.grad.elem(indices_sigma.subvec(i*J*J, (i+1)*J*J-1)) += arma::vectorise(grad_sigma);
    }
  }

  void dtransform(arguments_optim& x) {
    dmu = arma::reshape(x.dtransparameters.elem(indices_mu), J, I);
    dSigma.set_size(J, J, I);
    arma::vec dloglik = x.dtransparameters.elem(indices_loglik);
    for(int i = 0; i < I; ++i) {
      dSigma.slice(i) = arma::reshape(
        x.dtransparameters.elem(indices_sigma.subvec(i*J*J, (i+1)*J*J-1)), J, J);
      arma::vec dmui = dmu.col(i);
      for(const auto& block : patterns) {
        const arma::mat& score = block.mean_score[i];
        arma::mat dSig = dSigma.slice(i).submat(block.obs, block.obs);
        dloglik.elem(block.rows + i*S) += score.t() * dmui.elem(block.obs) +
          0.5 * (arma::sum(score % (dSig * score), 0).t() - arma::accu(dSig % block.precision[i]));
      }
    }
    x.dtransparameters.elem(indices_out) = dloglik;
  }

  void update_dgrad(arguments_optim& x) {
    for(int i = 0; i < I; ++i) {
      arma::vec grad_mu(J, arma::fill::zeros), dmui = dmu.col(i);
      arma::mat grad_sigma(J, J, arma::fill::zeros);
      for(const auto& block : patterns) {
        arma::uvec out = indices_out.elem(block.rows + i*S);
        arma::vec w = x.grad.elem(out), dw = x.dgrad.elem(out);
        const arma::mat& score = block.mean_score[i];
        const arma::mat& K = block.precision[i];
        arma::mat dSig = dSigma.slice(i).submat(block.obs, block.obs);
        arma::mat dK = -K * dSig * K;
        arma::mat shift = dSig * score;
        shift.each_col() += dmui.elem(block.obs);
        arma::mat dscore = -K * shift;
        grad_mu.elem(block.obs) += score * dw + dscore * w;
        arma::mat cross = (dscore.each_row() % w.t()) * score.t();
        grad_sigma.submat(block.obs, block.obs) += 0.5 *
          ((score.each_row() % dw.t()) * score.t() + cross + cross.t() -
           arma::accu(dw)*K - arma::accu(w)*dK);
      }
      x.dgrad.elem(indices_mu.subvec(i*J, (i+1)*J-1)) += grad_mu;
      x.dgrad.elem(indices_sigma.subvec(i*J*J, (i+1)*J*J-1)) += arma::vectorise(grad_sigma);
    }
  }

  void jacobian(arguments_optim& x) {
    jacob.zeros(n_out, n_in);
    for(int i = 0; i < I; ++i) {
      for(const auto& block : patterns) {
        arma::uvec rows = block.rows + i*S;
        const arma::mat& score = block.mean_score[i];
        jacob.submat(rows, block.obs + i*J) = score.t();
        for(arma::uword b = 0; b < block.obs.n_elem; ++b) {
          for(arma::uword a = 0; a < block.obs.n_elem; ++a) {
            arma::uword col = J*I + i*J*J + block.obs[a] + block.obs[b]*J;
            jacob.submat(rows, arma::uvec{col}) =
              0.5 * (score.row(a) % score.row(b) - block.precision[i](a,b)).t();
          }
        }
      }
    }
    jacob.cols(n_in-n_out, n_in-1).eye();
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

  mytrans->indices_mu = indices_mu;
  mytrans->indices_sigma = indices_sigma;
  mytrans->indices_loglik = indices_loglik;
  mytrans->indices_in = arma::join_cols(indices_mu,
                               arma::join_cols(indices_sigma, indices_loglik));
  mytrans->indices_out = indices_out[0];
  mytrans->n_in = n_in;
  mytrans->n_out = n_out;
  mytrans->S = S;
  mytrans->J = J;
  mytrans->I = I;
  mytrans->initialize_patterns(y);

  return mytrans;

}
