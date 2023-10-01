/*
 * Author: Marcos Jimenez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 28/05/2022
 *
 */

// #include "structures.h"
// #include "auxiliary_manifolds.h"

// Manifolds

class cor_manifold {

public:

  virtual void param(arguments_cor& x) = 0;

  virtual void dcor(arguments_cor& x) = 0;

  virtual void grad(arguments_cor& x) = 0;

  virtual void dgrad(arguments_cor& x) = 0;

  virtual void proj(arguments_cor& x) = 0;

  virtual void hess(arguments_cor& x) = 0;

  virtual void retr(arguments_cor& x) = 0;

};

// Oblique manifold:

class cor_oblq:public cor_manifold {

public:

  void param(arguments_cor& x) {

    x.correlation = x.T.t() * x.T;

  }

  void dcor(arguments_cor& x) {

    x.dcor = x.T.t() * x.dT;
    x.dcor += x.dcor.t();

  }

  void grad(arguments_cor& x) {

    x.g = x.T * x.gcor;

  }

  void dgrad(arguments_cor& x) {

    x.dg = x.dT * x.gcor + x.T * x.dgcor; // x.dgcor == (x.hcor % x.dcor) ¿?

  }

  void proj(arguments_cor& x) {

    x.rg = x.g - x.T * arma::diagmat( x.T.t() * x.g );

  }

  void hess(arguments_cor& x) {

    x.dH = x.dg - x.dT * arma::diagmat( x.T.t() * x.g) - x.T * arma::diagmat( x.T.t() * x.dg );
    // arma::mat drg = x.dg - x.dT * arma::diagmat( x.T.t() * x.g) - x.T * arma::diagmat( x.dT.t() * x.g ) -
      // x.T * arma::diagmat( x.T.t() * x.dg );
    // x.dH = drg - x.T * arma::diagmat( x.T.t() * drg );

  }

  void retr(arguments_cor& x) {

    x.T *= arma::diagmat( 1 / sqrt(arma::sum(x.T % x.T, 0)) );

  }

};

// Choose the manifold:

cor_manifold* choose_cor_manifold(std::string projection) {

  cor_manifold* manifold;
  if(projection == "oblq") {
    manifold = new cor_oblq();
  } else if(projection == "none") {

  } else {

    Rcpp::stop("Available projections: \n oblq, poblq");

  }

  return manifold;

}
