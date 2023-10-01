/*
 * Author: Marcos Jimenez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 15/09/2023
 *
 */

// Manifolds

class cfa_manifold {

public:

  virtual void param(arguments_cfa& x) = 0;

  virtual void dparam(arguments_cfa& x) = 0;

  virtual void grad(arguments_cfa& x) = 0;

  virtual void dgrad(arguments_cfa& x) = 0;

  virtual void proj(arguments_cfa& x) = 0;

  virtual void hess(arguments_cfa& x) = 0;

  virtual void retr(arguments_cfa& x) = 0;

};

// Identity manifold for CFA:

class cfa_id:public cfa_manifold {

public:

  void param(arguments_cfa& x) {

    // x.transformed = 0.5*(x.lower + x.upper) + 0.5*abs(x.upper - x.lower) % sin(x.parameters);
    // x.transformed = x.parameters;
    // x.transformed(x.target_positive) = x.parameters(x.target_positive) % x.parameters(x.target_positive);
    // x.lambda.elem(x.target_indexes) = x.transformed(x.lambda_indexes);
    // x.phi.elem(x.targetphi_indexes) = x.transformed(x.phi_indexes);
    // x.psi.elem(x.targetpsi_indexes) = x.transformed(x.psi_indexes);
    x.lambda.elem(x.target_indexes) = x.parameters(x.lambda_indexes);
    x.phi.elem(x.targetphi_indexes) = x.parameters(x.phi_indexes);
    x.psi.elem(x.targetpsi_indexes) = x.parameters(x.psi_indexes);
    x.phi = arma::symmatl(x.phi);
    x.psi = arma::symmatl(x.psi);

  }

  void dparam(arguments_cfa& x) {

    x.dlambda.elem(x.target_indexes) = x.dparameters(x.lambda_indexes);
    x.dphi.elem(x.targetphi_indexes) = x.dparameters(x.phi_indexes);
    x.dpsi.elem(x.targetpsi_indexes) = x.dparameters(x.psi_indexes);
    x.dphi = arma::symmatl(x.dphi);
    x.dpsi = arma::symmatl(x.dpsi);

  }

  void grad(arguments_cfa& x) {

    // x.g = x.gradient % (0.5*abs(x.upper - x.lower) % cos(x.parameters));
    x.g = x.gradient;
    // x.g(x.target_positive) %= 2*x.parameters(x.target_positive);

  }

  void dgrad(arguments_cfa& x) {

    // x.dg = x.dgradient % (0.5*abs(x.upper - x.lower) % cos(x.parameters)) +
    //   x.gradient % (0.5*abs(x.upper - x.lower) % -sin(x.dparameters));
    x.dg = x.dgradient;
    // x.dg(x.target_positive) *= 2;

  }

  void proj(arguments_cfa& x) {

    x.rg = x.g;

  }

  void hess(arguments_cfa& x) {

    x.dH = x.dg;

  }

  void retr(arguments_cfa& x) {

  }

};

// Positive-definite

class cfa_pos:public cfa_manifold {

public:

  void param(arguments_cfa& x) {

    x.lambda.elem(x.target_indexes) = x.parameters(x.lambda_indexes);
    x.T.elem(x.targetphi_indexes) = x.parameters(x.phi_indexes);
    x.phi = x.T.t() * x.T;
    x.psi.elem(x.targetpsi_indexes) = x.parameters(x.psi_indexes);
    x.psi = arma::symmatl(x.psi);
    // x.U.elem(x.targetpsi_indexes) = x.parameters(x.psi_indexes);
    // x.psi = x.U.t() * x.U;

  }

  void dparam(arguments_cfa& x) {

    x.dlambda.elem(x.target_indexes) = x.dparameters(x.lambda_indexes);
    x.dT.elem(x.targetphi_indexes) = x.dparameters(x.phi_indexes);
    x.dphi = x.T.t() * x.dT;
    x.dphi += x.dphi.t();
    // x.dU.elem(x.targetpsi_indexes) = x.dparameters(x.psi_indexes);
    // x.dpsi = x.U.t() * x.dU;
    // x.dpsi += x.dpsi.t();
    x.dpsi.elem(x.targetpsi_indexes) = x.dparameters(x.psi_indexes);
    x.dpsi = arma::symmatl(x.dpsi);

  }

  void grad(arguments_cfa& x) {

    x.g = x.gradient;
    x.gT = x.T * x.gphi;
    x.g(x.phi_indexes) = arma::vectorise(x.gT);
    // x.gU = x.U * x.gpsi;
    // x.g(x.psi_indexes) = arma::vectorise(x.gU);

  }

  void dgrad(arguments_cfa& x) {

    x.dg = x.dgradient;
    x.dgT = x.dT * x.gphi + x.T * x.dgphi;
    x.dg(x.phi_indexes) = arma::vectorise(x.dgT);
    // x.dgU = x.dU * x.gpsi + x.U * x.dgpsi;
    // x.dg(x.psi_indexes) = arma::vectorise(x.dgU);

  }

  void proj(arguments_cfa& x) {

    x.rg = x.g;

    arma::mat c1 = x.T.t() * x.gT;
    arma::mat X0 = c1 + c1.t();
    x.A = lyap_sym(x.phi, X0);
    x.A(x.oblq_indexes).zeros();
    arma::mat N = x.T * x.A;
    x.rg(x.phi_indexes) = arma::vectorise(x.gT - N);

    // c1 = x.U.t() * x.gU;
    // X0 = c1 + c1.t();
    // x.AU = lyap_sym(x.psi, X0);
    // x.AU(x.oblq_indexesU).zeros();
    // N = x.U * x.AU;
    // x.rg(x.psi_indexes) = arma::vectorise(x.gU - N);

  }

  void hess(arguments_cfa& x) {

    x.dH = x.dg;

    arma::mat dc1 = x.dT.t() * x.gT + x.T.t() * x.dgT; // Differential of c1
    arma::mat dX0 = dc1 + dc1.t(); // Differential of X0
    arma::mat c2 = x.A * x.dphi + x.dphi * x.A; // Differential of APhi + PhiA wrt Phi
    arma::mat Q = dX0 - c2;
    arma::mat dA = lyap_sym(x.phi, Q);
    dA(x.oblq_indexes).zeros();
    arma::mat drg = x.dgT - (x.dT * x.A + x.T * dA);
    // projection
    arma::mat c = x.T.t() * drg;
    arma::mat X0 = c + c.t();
    arma::mat A = lyap_sym(x.phi, X0);
    A(x.oblq_indexes).zeros();
    arma::mat N = x.T * A;
    x.dH(x.phi_indexes) = arma::vectorise(drg - N);

    // dc1 = x.dU.t() * x.gU + x.U.t() * x.dgU; // Differential of c1
    // dX0 = dc1 + dc1.t(); // Differential of X0
    // c2 = x.AU * x.dpsi + x.dpsi * x.AU; // Differential of APhi + PhiA wrt Phi
    // Q = dX0 - c2;
    // arma::mat dAU = lyap_sym(x.psi, Q);
    // dAU(x.oblq_indexesU).zeros();
    // drg = x.dgU - (x.dU * x.AU + x.U * dAU);
    // // projection
    // c = x.U.t() * drg;
    // X0 = c + c.t();
    // arma::mat AU = lyap_sym(x.psi, X0);
    // AU(x.oblq_indexesU).zeros();
    // N = x.U * AU;
    // x.dH(x.psi_indexes) = arma::vectorise(drg - N);

  }

  void retr(arguments_cfa& x) {

    int J = x.T.n_cols;

    for(int i=1; i < J; ++i) {

      arma::uvec indexes = consecutive(0, i-1L);
      arma::vec column = x.Phi_Target.col(i);
      arma::vec upper_column = column(indexes);
      arma::uvec zeros = arma::find(upper_column == 0);

      arma::mat Q;
      arma::mat R;
      qr_econ(Q, R, x.T.cols(zeros));
      arma::mat orthogonals = Q;

      x.T.col(i) = orthogonalize(orthogonals, x.T.col(i));

    }

    arma::vec z = arma::diagvec(x.T.t() * x.T);
    if(!x.free_indices_phi.is_empty()) z(x.free_indices_phi).ones();
    x.T *= arma::diagmat(1/sqrt(z));

    x.parameters(x.phi_indexes) = arma::vectorise(x.T);

  }

};

// Identity manifold for EFA:

class cfa_identity:public cfa_manifold {

public:

  void param(arguments_cfa& x) {

    x.psi2 = x.parameters;

  }

  void dparam(arguments_cfa& x) {

  }

  void grad(arguments_cfa& x) {

    x.g = x.g_psi2;

  }

  void dgrad(arguments_cfa& x) {

    Rcpp::stop("The differential of this estimator is not available yet.");

  }

  void proj(arguments_cfa& x) {

    x.rg = x.g;

  }

  void hess(arguments_cfa& x) {

  }

  void retr(arguments_cfa& x) {

  }

};

// Box-constraint manifold:

class cfa_box:public cfa_manifold {

public:

  void param(arguments_cfa& x) {

    x.psi2 = 0.5*(x.lower + x.upper) + 0.5*abs(x.upper - x.lower) % sin(x.parameters);

  }

  void dparam(arguments_cfa& x) {

  }

  void grad(arguments_cfa& x) {

    x.g = x.g_psi2 % (0.5*abs(x.upper - x.lower) % cos(x.parameters));

  }

  void dgrad(arguments_cfa& x) {

    Rcpp::stop("The differential of this estimator is not available yet.");

  }

  void proj(arguments_cfa& x) {

    x.rg = x.g;

  }

  void hess(arguments_cfa& x) {

  }

  void retr(arguments_cfa& x) {

  }

};

// DWLS manifold:

class efa_dwls_manifold:public cfa_manifold {

public:

  void param(arguments_cfa& x) {

    x.lambda.elem(x.lower_tri_ind) = x.parameters;

  }

  void dparam(arguments_cfa& x) {

  }

  void grad(arguments_cfa& x) {

    x.g = arma::vectorise(x.gL.elem(x.lower_tri_ind));

  }

  void dgrad(arguments_cfa& x) {

    // x.dg = x.G;

  }

  void proj(arguments_cfa& x) {

    x.rg = x.g;

  }

  void hess(arguments_cfa& x) {

    // x.dH = x.dg;

  }

  void retr(arguments_cfa& x) {

  }

};

// Choose the manifold:

cfa_manifold* choose_cfa_manifold(std::string projection) {

  cfa_manifold* manifold;
  if(projection == "id") {
    manifold = new cfa_id();
  } else if(projection == "positive") {
    manifold = new cfa_pos();
  } else if(projection == "identity") {
    manifold = new cfa_identity();
  } else if(projection == "dwls") {
    manifold = new efa_dwls_manifold();
  } else if(projection == "box") {
    manifold = new cfa_box();
  } else if(projection == "none") {

  } else {

    Rcpp::stop("Available projections: \n id, identity, box, and dwls_efa");

  }

  return manifold;

}

class cfa_manifold2 {

public:

  virtual void param(arguments_optim& x, std::vector<arguments_cfa>& structs) = 0;

  virtual void dparam(arguments_optim& x, std::vector<arguments_cfa>& structs) = 0;

  virtual void grad(arguments_optim& x, std::vector<arguments_cfa>& structs) = 0;

  virtual void dgrad(arguments_optim& x, std::vector<arguments_cfa>& structs) = 0;

  virtual void proj(arguments_optim& x, std::vector<arguments_cfa>& structs) = 0;

  virtual void hess(arguments_optim& x, std::vector<arguments_cfa>& structs) = 0;

  virtual void retr(arguments_optim& x, std::vector<arguments_cfa>& structs) = 0;

};

class ultimate_manifold: public cfa_manifold2 {

public:

  void param(arguments_optim& x, std::vector<arguments_cfa>& structs) {

    cfa_manifold* manifold;

    for(int i=0; i < x.nblocks; ++i) {

      structs[i].parameters = x.parameters;
      manifold = choose_cfa_manifold(structs[i].projection);
      manifold->param(structs[i]);

    }

  }

  void dparam(arguments_optim& x, std::vector<arguments_cfa>& structs) {

    cfa_manifold* manifold;

    for(int i=0; i < x.nblocks; ++i) {

      structs[i].dparameters = x.dparameters;
      manifold = choose_cfa_manifold(structs[i].projection);
      manifold->dparam(structs[i]);

    }

  }

  void grad(arguments_optim& x, std::vector<arguments_cfa>& structs) {

    cfa_manifold* manifold;
    x.g.set_size(x.parameters.n_elem); x.g.zeros();

    for(int i=0; i < x.nblocks; ++i) {

      manifold = choose_cfa_manifold(structs[i].projection);
      manifold->grad(structs[i]);
      x.g += structs[i].g;

    }

  }

  void dgrad(arguments_optim& x, std::vector<arguments_cfa>& structs) {

    cfa_manifold* manifold;
    x.dg.set_size(x.parameters.n_elem); x.dg.zeros();

    for(int i=0; i < x.nblocks; ++i) {

      manifold = choose_cfa_manifold(structs[i].projection);
      manifold->dgrad(structs[i]);
      x.dg += structs[i].dg;

    }

  }

  void proj(arguments_optim& x, std::vector<arguments_cfa>& structs) {

    cfa_manifold* manifold;
    x.rg.set_size(x.parameters.n_elem); x.rg.zeros();
    x.rg = x.g;

    for(int i=0; i < x.nblocks; ++i) {

      structs[i].g = x.g;
      manifold = choose_cfa_manifold(structs[i].projection);
      manifold->proj(structs[i]);
      x.rg(structs[i].phi_indexes) = structs[i].rg(structs[i].phi_indexes);

    }

    // x.rg = structs[0].rg;

  }

  void hess(arguments_optim& x, std::vector<arguments_cfa>& structs) {

    cfa_manifold* manifold;
    x.dH.set_size(x.parameters.n_elem); x.dH.zeros();

    for(int i=0; i < x.nblocks; ++i) {

      structs[i].dg = x.dg;
      manifold = choose_cfa_manifold(structs[i].projection);
      manifold->hess(structs[i]);

    }

    x.dH = structs[0].dH;

  }

  void retr(arguments_optim& x, std::vector<arguments_cfa>& structs) {

    cfa_manifold* manifold;

    // for(int i=0; i < x.nblocks; ++i) {
    //   structs[i].parameters = x.parameters;
    // }
    for(int i=0; i < x.nblocks; ++i) {

      structs[i].parameters = x.parameters;
      manifold = choose_cfa_manifold(structs[i].projection);
      manifold->param(structs[i]);
      manifold->retr(structs[i]);
      x.parameters(structs[i].phi_indexes) = structs[i].parameters(structs[i].phi_indexes);

    }
    // for(int i=0; i < x.nblocks; ++i) {
    //   x.parameters(structs[i].phi_indexes) = structs[i].parameters(structs[i].phi_indexes);
    // }

  }

};
