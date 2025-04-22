// Partially Oblique manifold:

class poblq:public manifolds {

public:

  arma::mat X = arma::mat(q, q);
  arma::mat dX = arma::mat(q, q);
  arma::mat A, Phi, dP;
  arma::uvec oblq_indices;

  void param() {

    X = arma::reshape(parameters, q, q);
    Phi = X.t() * X;

  }

  void proj() {

    g.reshape(q, q);
    arma::mat c1 = X.t() * g;
    arma::mat X0 = c1 + c1.t();
    A = lyap_sym(Phi, X0);
    A(oblq_indices).zeros();
    arma::mat N = X * A;
    rg = g - N;

  }

  void hess() {

    g.reshape(q, q);
    dg.reshape(q, q);
    dX = arma::reshape(dparameters, q, q);
    dP = X.t() * dX;
    dP += dX.t();

    // Implicit differentiation of APhi + PhiA = X0
    arma::mat dc1 = dX.t() * g + X.t() * dg; // Differential of c1
    arma::mat dX0 = dc1 + dc1.t(); // Differential of X0
    arma::mat c2 = A * dP + dP * A; // Differential of APhi + PhiA wrt Phi
    arma::mat Q = dX0 - c2;
    // dAPhi + PhidA = Q
    arma::mat dA = lyap_sym(Phi, Q);
    dA(oblq_indices).zeros();
    arma::mat drg = dg - (dX * A + X * dA);

    // projection
    arma::mat c = X.t() * drg;
    arma::mat X0 = c + c.t();
    arma::mat A = lyap_sym(Phi, X0);
    A(oblq_indices).zeros();
    arma::mat N = X * A;
    dH = drg - N;

  }

  void retr() {

    int J = X.n_cols;

    for(int i=1; i < J; ++i) {

      arma::uvec indexes = consecutive(0, i);
      arma::vec column = PhiTarget.col(i);
      arma::vec upper_column = column(indexes);
      arma::uvec zeros = arma::find(upper_column == 0);

      arma::mat Q;
      arma::mat R;
      qr_econ(Q, R, X.cols(zeros));
      arma::mat orthogonals = Q;
      // int n = orthogonals.n_cols;

      X.col(i) = orthogonalize(orthogonals, X.col(i));

    }

    parameters = arma::vectorise(X * arma::diagmat(1 / sqrt(arma::diagvec(X.t() * X))));

  }

};

