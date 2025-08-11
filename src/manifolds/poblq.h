/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 14/07/2025
 */

// Partially Oblique manifold:

class poblq:public manifolds {

public:

  std::size_t q;
  arma::mat X = arma::mat(q, q);
  arma::mat dX = arma::mat(q, q);
  arma::mat A, Phi, dP;
  arma::uvec oblq_indices;
  arma::mat target;

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

    arma::vec indicator = target.diag();
    arma::mat target2 = target; // Do not modify the original target
    target2.diag().ones();
    int J = X.n_cols;

    for(int i=1; i < J; ++i) {

      arma::uvec indexes = consecutive(0, i);
      arma::vec column = target2.col(i);
      arma::vec upper_column = column(indexes);
      arma::uvec zeros = arma::find(upper_column == 0);

      arma::mat Q;
      arma::mat R;
      qr_econ(Q, R, X.cols(zeros));

      X.col(i) = orthogonalize(Q, X.col(i));

    }

    arma::mat Xstd = X * arma::diagmat(1 / sqrt(arma::diagvec(X.t() * X)));
    arma::uvec ones = arma::find(indicator == 0);
    X.cols(ones) = Xstd.cols(ones);
    parameters = arma::vectorise(X);
    // parameters = arma::vectorise(X * arma::diagmat(1 / sqrt(arma::diagvec(X.t() * X))));

  }

  void dconstraints() {

  }

  void outcomes() {

  }

};

poblq* choose_poblq(Rcpp::List manifold_setup) {

  poblq* mymanifold = new poblq();

  // Provide these:
  std::vector<arma::uvec> indices = manifold_setup["indices"];
  arma::mat target = manifold_setup["target"];

  mymanifold->indices = indices;
  mymanifold->q = target.n_cols;
  mymanifold->target = target;

  // target.diag() += 10;
  arma::uvec oblq_indices = arma::find(target == 1);
  mymanifold->oblq_indices = oblq_indices;

  return mymanifold;

}

arma::mat poblq(arma::mat X, arma::mat target) {

  arma::vec indicator = target.diag();
  target.diag().ones();
  int J = X.n_cols;

  for(int i=1; i < J; ++i) {

    arma::uvec indexes = consecutive(0, i);
    arma::vec column = target.col(i);
    arma::vec upper_column = column(indexes);
    arma::uvec zeros = arma::find(upper_column == 0);

    arma::mat Q;
    arma::mat R;
    qr_econ(Q, R, X.cols(zeros));

    X.col(i) = orthogonalize(Q, X.col(i));

  }

  arma::mat Xstd = X * arma::diagmat(1 / sqrt(arma::diagvec(X.t() * X)));
  arma::uvec ones = arma::find(indicator == 0);
  X.cols(ones) = Xstd.cols(ones);

  return X;

}

arma::mat rpoblq(int p, int q, arma::mat target) {

  arma::mat X(p, q, arma::fill::randn);
  X = poblq(X, target);

  return X;

}

