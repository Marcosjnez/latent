/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 24/08/2026
 */

// Oblique manifold:

class oblq:public manifolds {

public:

  arma::uvec indices;
  std::size_t p, q;
  arma::mat X, dX, g, dg;

  void param(arguments_optim& x) {

    if(indices.n_elem != p*q) {
      Rcpp::stop("The oblique manifold requires p*q parameters.");
    }

    X = arma::reshape(x.parameters(indices), p, q);

  }

  void proj(arguments_optim& x) {

    g = arma::reshape(x.g.elem(indices), p, q);
    x.rg.elem(indices) =
      arma::vectorise(g-X*arma::diagmat(X.t()*g));

  }

  void hess(arguments_optim& x) {

    g = arma::reshape(x.g.elem(indices), p, q);
    dg = arma::reshape(x.dg.elem(indices), p, q);
    dX = arma::reshape(x.dparameters.elem(indices), p, q);

    x.dH.elem(indices) =
      arma::vectorise(
        dg-dX*arma::diagmat(X.t()*g)-
        X*arma::diagmat(X.t()*dg)
      );

    // arma::mat drg = dg - dX * arma::diagmat( X.t() * g) - X * arma::diagmat(dX.t() * g) -
    // X * arma::diagmat(X.t() * dg);
    // dH = drg - X * arma::diagmat(X.t() * drg);

  }

  void retr(arguments_optim& x) {

    x.parameters(indices) =
      arma::vectorise(
        X*arma::diagmat(1/sqrt(arma::sum(X % X, 0)))
      );

  }

  void tangent_basis(arguments_optim& x) {

    X = arma::reshape(x.parameters(indices), p, q);

    const arma::uword dimension =
      q*(p > 0L ? p-1L : 0L);

    T.zeros(x.nparam, dimension);

    arma::uword column = 0L;

    for(arma::uword j=0L; j < q; ++j) {

      arma::mat T_column =
        tangent_complement(X.col(j));

      for(arma::uword k=0L; k < T_column.n_cols; ++k) {

        for(arma::uword i=0L; i < p; ++i) {
          arma::uword local_index = i+j*p;
          T(indices[local_index], column) =
            T_column(i, k);
        }

        column++;

      }

    }

  }

  void dconstraints(arguments_optim& x) {

    X = arma::reshape(x.parameters(indices), p, q);

    arma::uvec trans_indices =
      x.transparam2param.elem(indices);

    arma::sp_mat first_derivatives(
      x.transparameters.n_elem, q
    );
    std::vector<arma::sp_mat> second_derivatives;
    second_derivatives.reserve(q);

    for(arma::uword j=0L; j < q; ++j) {

      arma::sp_mat second_derivative(
        x.transparameters.n_elem,
        x.transparameters.n_elem
      );

      for(arma::uword i=0L; i < p; ++i) {

        arma::uword local_index = i+j*p;
        arma::uword trans_index =
          trans_indices[local_index];

        first_derivatives(trans_index, j) =
          2.00*X(i, j);
        second_derivative(trans_index,
                          trans_index) = 2.00;

      }

      second_derivatives.push_back(
        second_derivative
      );

    }

    append_constraint_derivatives(
      x,
      first_derivatives,
      second_derivatives
    );

  }

  void outcomes(arguments_optim& x) {

    tangent_basis(x);

    matrices.resize(1);
    matrices[0] = T;
    names_matrices.resize(1);
    names_matrices[0] = "tangent_basis";

  }

};

oblq* choose_oblq(Rcpp::List manifold_setup) {

  oblq* mymanifold = new oblq();

  arma::uvec indices = manifold_setup["indices"];
  std::size_t p = manifold_setup["p"];
  std::size_t q = manifold_setup["q"];

  mymanifold->indices = indices;
  mymanifold->p = p;
  mymanifold->q = q;

  return mymanifold;

}

arma::mat oblq(arma::mat X) {

  X *= arma::diagmat(1 / sqrt(arma::sum(X % X, 0)));

  return X;

}

arma::mat roblq(int p, int q) {

  arma::mat X(p, q, arma::fill::randn);
  X = oblq(X);

  return X;

}
