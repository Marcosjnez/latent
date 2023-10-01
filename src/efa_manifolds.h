/*
 * Author: Marcos Jimenez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 14/09/2023
 *
 */

// Manifolds

class efa_manifold {

public:

  virtual void param(arguments_efa& x) = 0;

  virtual void dparam(arguments_efa& x) = 0;

  virtual void grad(arguments_efa& x) = 0;

  virtual void dgrad(arguments_efa& x) = 0;

  virtual void proj(arguments_efa& x) = 0;

  virtual void hess(arguments_efa& x) = 0;

  virtual void retr(arguments_efa& x) = 0;

};

// Identity manifold:

class identity:public efa_manifold {

public:

  void param(arguments_efa& x) {

    x.psi2 = x.parameters;

  }

  void dparam(arguments_efa& x) {

  }

  void grad(arguments_efa& x) {

    x.g = x.g_psi2;

  }

  void dgrad(arguments_efa& x) {

    Rcpp::stop("The differential of this estimator is not available yet.");

  }

  void proj(arguments_efa& x) {

    x.rg = x.g;

  }

  void hess(arguments_efa& x) {

  }

  void retr(arguments_efa& x) {

  }

};

// Box-constraint manifold:

class box:public efa_manifold {

public:

  void param(arguments_efa& x) {

    x.psi2 = 0.5*(x.lower + x.upper) + 0.5*abs(x.upper - x.lower) % sin(x.parameters);

  }

  void dparam(arguments_efa& x) {

  }

  void grad(arguments_efa& x) {

    x.g = x.g_psi2 % (0.5*abs(x.upper - x.lower) % cos(x.parameters));

  }

  void dgrad(arguments_efa& x) {

    Rcpp::stop("The differential of this estimator is not available yet.");

  }

  void proj(arguments_efa& x) {

    x.rg = x.g;

  }

  void hess(arguments_efa& x) {

  }

  void retr(arguments_efa& x) {

  }

};

// DWLS manifold:

class dwls_efa:public efa_manifold {

public:

  void param(arguments_efa& x) {

    x.lambda.elem(x.lower_tri_ind) = x.parameters;

  }

  void dparam(arguments_efa& x) {

  }

  void grad(arguments_efa& x) {

    x.g = arma::vectorise(x.gL.elem(x.lower_tri_ind));

  }

  void dgrad(arguments_efa& x) {

    // x.dg = x.G;

  }

  void proj(arguments_efa& x) {

    x.rg = x.g;

  }

  void hess(arguments_efa& x) {

    // x.dH = x.dg;

  }

  void retr(arguments_efa& x) {

  }

};

// Choose the manifold:

efa_manifold* choose_efa_manifold(std::string mani) {

  efa_manifold* manifold;
  if(mani == "identity") {
    manifold = new identity();
  } else if(mani == "dwls") {
    manifold = new dwls_efa();
  } else if(mani == "box") {
    manifold = new box();
  } else {

    Rcpp::stop("Available manifolds for factor extraction: \n identity, dwls, box");

  }

  return manifold;

}
