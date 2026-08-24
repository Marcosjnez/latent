/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 24/08/2026
 */

// Orthogonal manifold:

class orth:public manifolds {

public:

  arma::uvec indices;
  std::size_t p, q;
  arma::mat X, dX, g, dg;

  void param(arguments_optim& x) {

    if(indices.n_elem != p*q) {
      Rcpp::stop("The orthogonal manifold requires p*q parameters.");
    }

    X = arma::reshape(x.parameters(indices), p, q);

  }

  void proj(arguments_optim& x) {

    g = arma::reshape(x.g.elem(indices), p, q);
    x.rg.elem(indices) =
      arma::vectorise(g-X*symm(X.t()*g));

  }

  void hess(arguments_optim& x) {

    g = arma::reshape(x.g.elem(indices), p, q);
    dg = arma::reshape(x.dg.elem(indices), p, q);
    dX = arma::reshape(x.dparameters.elem(indices), p, q);

    arma::mat drg =
      dg-dX*symm(X.t()*g);

    x.dH.elem(indices) =
      arma::vectorise(drg-X*symm(X.t()*drg));

  }

  void retr(arguments_optim& x) {

    arma::mat Q, R;
    arma::qr_econ(Q, R, X);

    x.parameters(indices) = arma::vectorise(Q);

  }

  void tangent_basis(arguments_optim& x) {

    if(p < q) {
      Rcpp::stop("The orthogonal manifold requires p >= q.");
    }

    X = arma::reshape(x.parameters(indices), p, q);

    const arma::uword vertical_dimension =
      q*(q-1L)/2L;
    const arma::uword horizontal_dimension =
      (p-q)*q;
    const arma::uword dimension =
      vertical_dimension+horizontal_dimension;

    T.zeros(x.nparam, dimension);

    arma::uword column = 0L;
    const double inv_sqrt2 = 1.00/std::sqrt(2.00);

    for(arma::uword i=0L; i+1L < q; ++i) {
      for(arma::uword j=i+1L; j < q; ++j) {

        arma::mat Omega(q, q, arma::fill::zeros);
        Omega(i, j) = inv_sqrt2;
        Omega(j, i) = -inv_sqrt2;

        arma::vec direction =
          arma::vectorise(X*Omega);

        for(arma::uword k=0L; k < indices.n_elem; ++k) {
          T(indices[k], column) = direction[k];
        }

        column++;

      }
    }

    if(p > q) {

      arma::mat Q, R;
      arma::qr(Q, R, X);
      arma::mat X_perp = Q.cols(q, p-1L);

      for(arma::uword j=0L; j < q; ++j) {
        for(arma::uword i=0L; i < p-q; ++i) {

          arma::mat direction(p, q, arma::fill::zeros);
          direction.col(j) = X_perp.col(i);
          arma::vec direction_vec =
            arma::vectorise(direction);

          for(arma::uword k=0L; k < indices.n_elem; ++k) {
            T(indices[k], column) = direction_vec[k];
          }

          column++;

        }
      }

    }

  }

  void dconstraints(arguments_optim& x) {

    X = arma::reshape(x.parameters(indices), p, q);

    const arma::uword nconstraints =
      q*(q+1L)/2L;
    arma::uword first_column = x.dconstr.n_cols;
    x.dconstr.resize(x.transparameters.n_elem,
                     first_column+nconstraints);

    arma::uvec trans_indices =
      x.transparam2param.elem(indices);

    arma::uword constraint = 0L;

    for(arma::uword column=0L; column < q; ++column) {
      for(arma::uword row=column; row < q; ++row) {

        if(row == column) {

          for(arma::uword i=0L; i < p; ++i) {
            arma::uword local_index = i+column*p;
            x.dconstr(trans_indices[local_index],
                      first_column+constraint) =
              2.00*X(i, column);
          }

        } else {

          for(arma::uword i=0L; i < p; ++i) {

            arma::uword column_index = i+column*p;
            arma::uword row_index = i+row*p;

            x.dconstr(trans_indices[column_index],
                      first_column+constraint) =
              X(i, row);

            x.dconstr(trans_indices[row_index],
                      first_column+constraint) =
              X(i, column);

          }

        }

        constraint++;

      }
    }

  }

  void outcomes(arguments_optim& x) {

    tangent_basis(x);

    matrices.resize(1);
    matrices[0] = T;
    names_matrices.resize(1);
    names_matrices[0] = "tangent_basis";

  }

};

orth* choose_orth(Rcpp::List manifold_setup) {

  orth* mymanifold = new orth();

  arma::uvec indices = manifold_setup["indices"];
  std::size_t p = manifold_setup["p"];
  std::size_t q = manifold_setup["q"];

  mymanifold->indices = indices;
  mymanifold->p = p;
  mymanifold->q = q;

  return mymanifold;

}

arma::mat orth(arma::mat X) {

  arma::mat Q, R;
  arma::qr_econ(Q, R, X);

  return Q;

}

arma::mat rorth(int p, int q) {

  arma::mat X(p, q, arma::fill::randn);
  X = orth(X);

  return X;

}
