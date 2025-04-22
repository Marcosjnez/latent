/*
 * Author: Marcos Jimenez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 14/09/2024
 */

// Transformations

class transformations {

public:

  arma::vec parameters, dparameters, gradient, dgradient, g, dg;
  arma::uvec indices;
  int q;

  virtual void trans() = 0;

  virtual void dtrans() = 0;

  virtual void grad() = 0;

  virtual void dgrad() = 0;

};

// No transformation:

class none:public transformations {

public:

  void trans() {

  }

  void dtrans() {

  }

  void grad() {

    g = arma::vectorise(gradient);

  }

  void dgrad() {

    dg = arma::vectorise(dgradient);

  }

};

// Choose the transformation:

transformations* choose_trans(Rcpp::List trans_setup, transformations* xtrans) {

  transformations* trans;
  std::string transform = trans_setup["trans"];
  std::vector<arma::uvec> trans_indices;

  if(transform == "none") {

    none* mytrans = new none();
    trans = mytrans;

  } else if(transform == "exp") {


  } else {

    Rcpp::stop("Available transformations: \n none and crossprod");

  }

  return trans;

}

// Product transformation:

class product_transform {

public:

  void trans(arguments_optim& x, std::vector<transformations*>& xtransformations) {

    for(int i=0; i < x.ntransformations; ++i) {

      arma::uvec indices = xtransformations[i]->indices;
      xtransformations[i]->parameters = x.parameters.elem(indices);
      xtransformations[i]->trans();
      x.parameters.elem(indices) = xtransformations[i]->parameters;

    }

  }

  void dtrans(arguments_optim& x, std::vector<transformations*>& xtransformations) {

    for(int i=0; i < x.ntransformations; ++i) {

      arma::uvec indices = xtransformations[i]->indices;
      xtransformations[i]->dparameters = x.dparameters.elem(indices);
      xtransformations[i]->dtrans();
      x.dparameters.elem(indices) = xtransformations[i]->dparameters;

    }

  }

  void grad(arguments_optim& x, std::vector<transformations*>& xtransformations) {

    x.g.set_size(x.parameters.n_elem); x.g.zeros();

    for(int i=0; i < x.ntransformations; ++i) {

      arma::uvec indices = xtransformations[i]->indices;
      xtransformations[i]->gradient = x.gradient(indices);
      xtransformations[i]->grad();
      x.g.elem(indices) += xtransformations[i]->gradient;

    }

  }

  void dgrad(arguments_optim& x, std::vector<transformations*>& xtransformations) {

    x.dg.set_size(x.parameters.n_elem); x.dg.zeros();

    for(int i=0; i < x.ntransformations; ++i) {

      arma::uvec indices = xtransformations[i]->indices;
      xtransformations[i]->dgradient = x.dgradient(indices);
      xtransformations[i]->dgrad();
      x.dg.elem(indices) += xtransformations[i]->dgradient;

    }

  }

};
