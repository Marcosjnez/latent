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

    x.g = x.gradient;

  }

  void dgrad(arguments_cfa& x) {

    x.dg = x.dgradient;

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
    // x.psi.elem(x.targetpsi_indexes) = x.parameters(x.psi_indexes);
    // x.psi = arma::symmatl(x.psi);
    x.U.elem(x.targetpsi_indexes) = x.parameters(x.psi_indexes);
    x.psi = x.U.t() * x.U;

  }

  void dparam(arguments_cfa& x) {

    x.dlambda.elem(x.target_indexes) = x.dparameters(x.lambda_indexes);
    x.dT.elem(x.targetphi_indexes) = x.dparameters(x.phi_indexes);
    x.dphi = x.T.t() * x.dT;
    x.dphi += x.dphi.t();
    x.dU.elem(x.targetpsi_indexes) = x.dparameters(x.psi_indexes);
    x.dpsi = x.U.t() * x.dU;
    x.dpsi += x.dpsi.t();
    // x.dpsi.elem(x.targetpsi_indexes) = x.dparameters(x.psi_indexes);
    // x.dpsi = arma::symmatl(x.dpsi);

  }

  void grad(arguments_cfa& x) {

    x.g = x.gradient;
    x.gT = x.T * x.gphi;
    x.g(x.phi_indexes) = arma::vectorise(x.gT);
    x.gU = x.U * x.gpsi;
    x.g(x.psi_indexes) = arma::vectorise(x.gU);

  }

  void dgrad(arguments_cfa& x) {

    x.dg = x.dgradient;
    x.dgT = x.dT * x.gphi + x.T * x.dgphi;
    x.dg(x.phi_indexes) = arma::vectorise(x.dgT);
    x.dgU = x.dU * x.gpsi + x.U * x.dgpsi;
    x.dg(x.psi_indexes) = arma::vectorise(x.dgU);

  }

  void proj(arguments_cfa& x) {

    x.rg = x.g;

    arma::mat c1 = x.T.t() * x.gT;
    arma::mat X0 = c1 + c1.t();
    x.A = lyap_sym(x.phi, X0);
    x.A(x.oblq_indexes).zeros();
    arma::mat N = x.T * x.A;
    x.rg(x.phi_indexes) = arma::vectorise(x.gT - N);

    c1 = x.U.t() * x.gU;
    X0 = c1 + c1.t();
    x.AU = lyap_sym(x.psi, X0);
    x.AU(x.psi_oblq_indexes).zeros();
    N = x.U * x.AU;
    x.rg(x.psi_indexes) = arma::vectorise(x.gU - N);

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

    dc1 = x.dU.t() * x.gU + x.U.t() * x.dgU; // Differential of c1
    dX0 = dc1 + dc1.t(); // Differential of X0
    c2 = x.AU * x.dpsi + x.dpsi * x.AU; // Differential of APhi + PhiA wrt Phi
    Q = dX0 - c2;
    arma::mat dAU = lyap_sym(x.psi, Q);
    dAU(x.psi_oblq_indexes).zeros();
    drg = x.dgU - (x.dU * x.AU + x.U * dAU);
    // projection
    c = x.U.t() * drg;
    X0 = c + c.t();
    arma::mat AU = lyap_sym(x.psi, X0);
    AU(x.psi_oblq_indexes).zeros();
    N = x.U * AU;
    x.dH(x.psi_indexes) = arma::vectorise(drg - N);

  }

  void retr(arguments_cfa& x) {

    // Retraction for phi:

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

    // Retraction for psi:

    J = x.U.n_cols;

    for(int i=1; i < J; ++i) {

      arma::uvec indexes = consecutive(0, i-1L);
      arma::vec column = x.Psi_Target.col(i);
      arma::vec upper_column = column(indexes);
      arma::uvec zeros = arma::find(upper_column == 0);

      arma::mat Q;
      arma::mat R;
      qr_econ(Q, R, x.U.cols(zeros));
      arma::mat orthogonals = Q;

      x.U.col(i) = orthogonalize(orthogonals, x.U.col(i));

    }

    z = arma::diagvec(x.U.t() * x.U);
    if(!x.free_indices_psi.is_empty()) z(x.free_indices_psi).ones();
    x.U *= arma::diagmat(1/sqrt(z));

    x.parameters(x.psi_indexes) = arma::vectorise(x.U);

  }

};

// Choose the manifold:

cfa_manifold* choose_cfa_manifold(std::string projection) {

  cfa_manifold* manifold;
  if(projection == "id") {
    manifold = new cfa_id();
  } else if(projection == "positive") {
    manifold = new cfa_pos();
  } else if(projection == "none") {

  } else {

    Rcpp::stop("Available projections: \n id and positive");

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
      x.rg(structs[i].psi_indexes) = structs[i].rg(structs[i].psi_indexes);

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
      x.parameters(structs[i].psi_indexes) = structs[i].parameters(structs[i].psi_indexes);

    }
    // for(int i=0; i < x.nblocks; ++i) {
    //   x.parameters(structs[i].phi_indexes) = structs[i].parameters(structs[i].phi_indexes);
    // }

  }

};

// // Manifolds for rotation
//
// class rot_manifold {
//
// public:
//
//   virtual void param(arguments_rotate& x) = 0;
//
//   virtual void dparam(arguments_rotate& x) = 0;
//
//   virtual void grad(arguments_rotate& x) = 0;
//
//   virtual void dgrad(arguments_rotate& x) = 0;
//
//   virtual void proj(arguments_rotate& x) = 0;
//
//   virtual void hess(arguments_rotate& x) = 0;
//
//   virtual void retr(arguments_rotate& x) = 0;
//
// };
//
// // Orthogonal manifold:
//
// class rot_orth:public rot_manifold {
//
// public:
//
//   void param(arguments_rotate& x) {
//
//     x.T = x.parameters.reshape(x.q, x.q);
//
//   }
//
//   void dparam(arguments_rotate& x) {
//
//     x.dT = x.dparameters.reshape(x.q, x.q);
//
//   }
//
//   void grad(arguments_rotate& x) {
//
//     x.g = x.gradient;
//
//   }
//
//   void dgrad(arguments_rotate& x) {
//
//     x.dg = x.dgradient;
//
//   }
//
//   void proj(arguments_rotate& x) {
//
//     x.rg = x.T * skew(x.T.t() * x.g);
//
//   }
//
//   void hess(arguments_rotate& x) {
//
//     arma::mat drg = x.dg - x.dT * symm(x.T.t() * x.g); // update this?
//     x.dH = x.T * skew(x.T.t() * drg);
//
//   }
//
//   void retr(arguments_rotate& x) {
//
//     arma::mat Q, R;
//     arma::qr_econ(Q, R, x.T);
//
//     x.T = Q;
//
//   }
//
// };
//
// // Oblique manifold:
//
// class rot_oblq:public rot_manifold {
//
// public:
//
//   void param(arguments_rotate& x) {
//
//     x.T = x.parameters.reshape(x.q, x.q);
//
//   }
//
//   void dparam(arguments_rotate& x) {
//
//     x.dT = x.dparameters.reshape(x.q, x.q);
//
//   }
//
//   void grad(arguments_rotate& x) {
//
//     x.g = x.gradient;
//
//   }
//
//   void dgrad(arguments_rotate& x) {
//
//     x.dg = x.dgradient;
//
//   }
//
//   void proj(arguments_rotate& x) {
//
//     x.rg = x.g - x.T * arma::diagmat( x.T.t() * x.g );
//
//   }
//
//   void hess(arguments_rotate& x) {
//
//     x.dH = x.dg - x.dT * arma::diagmat( x.T.t() * x.g) - x.T * arma::diagmat( x.T.t() * x.dg );
//     // arma::mat drg = x.dg - x.dT * arma::diagmat( x.T.t() * x.g) - x.T * arma::diagmat( x.dT.t() * x.g ) -
//     // x.T * arma::diagmat( x.T.t() * x.dg );
//     // x.dH = drg - x.T * arma::diagmat( x.T.t() * drg );
//
//   }
//
//   void retr(arguments_rotate& x) {
//
//     x.T *= arma::diagmat( 1 / sqrt(arma::sum(x.T % x.T, 0)) );
//
//   }
//
// };
//
// // Partially oblique manifold:
//
// class rot_poblq:public rot_manifold {
//
// public:
//
//   void param(arguments_rotate& x) {
//
//     x.T = x.parameters.reshape(x.q, x.q);
//
//   }
//
//   void dparam(arguments_rotate& x) {
//
//     x.dT = x.dparameters.reshape(x.q, x.q);
//
//   }
//
//   void grad(arguments_rotate& x) {
//
//     x.g = x.gradient;
//
//   }
//
//   void dgrad(arguments_rotate& x) {
//
//     x.dg = x.dgradient;
//
//   }
//
//   void proj(arguments_rotate& x) {
//
//     arma::mat c1 = x.T.t() * x.g;
//     arma::mat X0 = c1 + c1.t();
//     x.A = lyap_sym(x.Phi, X0);
//     x.A(x.oblq_indexes).zeros();
//     arma::mat N = x.T * x.A;
//     x.rg = x.g - N;
//
//   }
//
//   void hess(arguments_rotate& x) {
//
//     // Implicit differentiation of APhi + PhiA = X0
//     arma::mat dc1 = x.dT.t() * x.g + x.T.t() * x.dg; // Differential of c1
//     arma::mat dX0 = dc1 + dc1.t(); // Differential of X0
//     arma::mat c2 = x.A * x.dP + x.dP * x.A; // Differential of APhi + PhiA wrt Phi
//     arma::mat Q = dX0 - c2;
//     // dAPhi + PhidA = Q
//     arma::mat dA = lyap_sym(x.Phi, Q);
//     dA(x.oblq_indexes).zeros();
//     arma::mat drg = x.dg - (x.dT * x.A + x.T * dA);
//
//     // projection
//     arma::mat c = x.T.t() * drg;
//     arma::mat X0 = c + c.t();
//     arma::mat A = lyap_sym(x.Phi, X0);
//     A(x.oblq_indexes).zeros();
//     arma::mat N = x.T * A;
//     x.dH = drg - N;
//
//   }
//
//   void retr(arguments_rotate& x) {
//
//     if(x.Phi_Target.is_empty()) {
//
//       arma::mat Q;
//       arma::mat R;
//       qr_econ(Q, R, x.T);
//
//       arma::mat X2 = x.T;
//       x.T = Q;
//       x.T.cols(x.list_oblq_indexes[0]) = X2.cols(x.list_oblq_indexes[0]);
//
//       int J = x.T.n_cols;
//       int I = x.list_oblq_indexes.size()-1;
//
//       // In the following loop, the Gram-Schmidt process is performed between blocks:
//
//       for(int i=0; i < I; ++i) {
//
//         // Select the cumulative indexes of blocks:
//         std::vector<arma::uvec> indexes(&x.list_oblq_indexes[0], &x.list_oblq_indexes[i+1]);
//         arma::uvec cum_indexes = list_to_vector(indexes);
//
//         arma::mat orthogonals = Q.cols(cum_indexes);
//         int n = orthogonals.n_cols;
//
//         for(int j=0; j < x.list_oblq_indexes[i+1].size(); ++j) {
//
//           // orthogonalize every column of the following block:
//           int index = x.list_oblq_indexes[i+1][j];
//           x.T.col(index) = orthogonalize(orthogonals, X2.col(index));
//
//         }
//
//       }
//
//     } else {
//
//       int J = x.T.n_cols;
//
//       for(int i=1; i < J; ++i) {
//
//         arma::uvec indexes = consecutive(0, i);
//         arma::vec column = x.Phi_Target.col(i);
//         arma::vec upper_column = column(indexes);
//         arma::uvec zeros = arma::find(upper_column == 0);
//
//         arma::mat Q;
//         arma::mat R;
//         qr_econ(Q, R, x.T.cols(zeros));
//         arma::mat orthogonals = Q;
//         int n = orthogonals.n_cols;
//
//         x.T.col(i) = orthogonalize(orthogonals, x.T.col(i));
//
//       }
//
//     }
//
//     x.T *= arma::diagmat( 1 / sqrt(arma::diagvec( x.T.t() * x.T )) );
//
//   }
//
// };
//
// rot_manifold* choose_rot_manifold(std::string projection) {
//
//   rot_manifold* manifold;
//   if(projection == "orth") {
//     manifold = new rot_orth();
//   } else if(projection == "oblq") {
//     manifold = new rot_oblq();
//   } else if(projection == "poblq") {
//     manifold = new rot_poblq();
//   } else if(projection == "none") {
//
//   } else {
//
//     Rcpp::stop("Available projections: \n orth, oblq, poblq");
//
//   }
//
//   return manifold;
//
// }
//
// class rot_manifold2 {
//
// public:
//
//   virtual void param(arguments_optim& x, std::vector<arguments_rotate>& structs) = 0;
//
//   virtual void dparam(arguments_optim& x, std::vector<arguments_rotate>& structs) = 0;
//
//   virtual void grad(arguments_optim& x, std::vector<arguments_rotate>& structs) = 0;
//
//   virtual void dgrad(arguments_optim& x, std::vector<arguments_rotate>& structs) = 0;
//
//   virtual void proj(arguments_optim& x, std::vector<arguments_rotate>& structs) = 0;
//
//   virtual void hess(arguments_optim& x, std::vector<arguments_rotate>& structs) = 0;
//
//   virtual void retr(arguments_optim& x, std::vector<arguments_rotate>& structs) = 0;
//
// };
//
// class ultimate_rot_manifold: public rot_manifold2 {
//
// public:
//
//   void param(arguments_optim& x, std::vector<arguments_rotate>& structs) {
//
//     rot_manifold* manifold;
//
//     for(int i=0; i < x.nblocks; ++i) {
//
//       structs[i].parameters = x.parameters;
//       manifold = choose_rot_manifold(structs[i].projection);
//       manifold->param(structs[i]);
//
//     }
//
//   }
//
//   void dparam(arguments_optim& x, std::vector<arguments_rotate>& structs) {
//
//     rot_manifold* manifold;
//
//     for(int i=0; i < x.nblocks; ++i) {
//
//       structs[i].dparameters = x.dparameters;
//       manifold = choose_rot_manifold(structs[i].projection);
//       manifold->dparam(structs[i]);
//
//     }
//
//   }
//
//   void grad(arguments_optim& x, std::vector<arguments_rotate>& structs) {
//
//     rot_manifold* manifold;
//     x.g.set_size(x.parameters.n_elem); x.g.zeros();
//
//     for(int i=0; i < x.nblocks; ++i) {
//
//       manifold = choose_rot_manifold(structs[i].projection);
//       manifold->grad(structs[i]);
//       x.g += structs[i].g;
//
//     }
//
//   }
//
//   void dgrad(arguments_optim& x, std::vector<arguments_rotate>& structs) {
//
//     rot_manifold* manifold;
//     x.dg.set_size(x.parameters.n_elem); x.dg.zeros();
//
//     for(int i=0; i < x.nblocks; ++i) {
//
//       manifold = choose_rot_manifold(structs[i].projection);
//       manifold->dgrad(structs[i]);
//       x.dg += structs[i].dg;
//
//     }
//
//   }
//
//   void proj(arguments_optim& x, std::vector<arguments_rotate>& structs) {
//
//     rot_manifold* manifold;
//     x.rg.set_size(x.parameters.n_elem); x.rg.zeros();
//     x.rg = x.g;
//
//     for(int i=0; i < x.nblocks; ++i) {
//
//       structs[i].g = x.g;
//       manifold = choose_rot_manifold(structs[i].projection);
//       manifold->proj(structs[i]);
//       x.rg(structs[i].phi_indexes) = structs[i].rg(structs[i].phi_indexes);
//
//     }
//
//     // x.rg = structs[0].rg;
//
//   }
//
//   void hess(arguments_optim& x, std::vector<arguments_rotate>& structs) {
//
//     rot_manifold* manifold;
//     x.dH.set_size(x.parameters.n_elem); x.dH.zeros();
//
//     for(int i=0; i < x.nblocks; ++i) {
//
//       structs[i].dg = x.dg;
//       manifold = choose_rot_manifold(structs[i].projection);
//       manifold->hess(structs[i]);
//
//     }
//
//     x.dH = structs[0].dH;
//
//   }
//
//   void retr(arguments_optim& x, std::vector<arguments_rotate>& structs) {
//
//     rot_manifold* manifold;
//
//     // for(int i=0; i < x.nblocks; ++i) {
//     //   structs[i].parameters = x.parameters;
//     // }
//     for(int i=0; i < x.nblocks; ++i) {
//
//       structs[i].parameters = x.parameters;
//       manifold = choose_rot_manifold(structs[i].projection);
//       manifold->param(structs[i]);
//       manifold->retr(structs[i]);
//       x.parameters(structs[i].phi_indexes) = structs[i].parameters(structs[i].phi_indexes);
//
//     }
//     // for(int i=0; i < x.nblocks; ++i) {
//     //   x.parameters(structs[i].phi_indexes) = structs[i].parameters(structs[i].phi_indexes);
//     // }
//
//   }
//
// };
