// Optimization algorithms for cfa

class cfa_optim {

public:

  virtual void optim(arguments_optim& x, std::vector<arguments_cfa>& structs) = 0;

};

// Riemannian gradient descent:

class cfa_RGD:public cfa_optim {

public:

  void optim(arguments_optim& x, std::vector<arguments_cfa>& structs) {

    gd(x, structs);

  }

};

// Riemannian Newton Trust-Region:

class cfa_RNTR:public cfa_optim {

public:

  void optim(arguments_optim& x, std::vector<arguments_cfa>& structs) {

    ntr(x, structs);

  }

};

// BFGS algorithm:

class cfa_BFGS:public cfa_optim {

public:

  void optim(arguments_optim& x, std::vector<arguments_cfa>& structs) {

    bfgs(x, structs);

  }

};

// L-BFGS algorithm:

class cfa_LBFGS:public cfa_optim {

public:

  void optim(arguments_optim& x, std::vector<arguments_cfa>& structs) {

    lbfgs(x, structs);

  }

};

// Select the optimization algorithm:

cfa_optim* choose_cfa_optim(std::string optim) {

  cfa_optim* algorithm;
  if(optim == "gradient") {
    algorithm = new cfa_RGD();
  } else if(optim == "newtonTR") {
    algorithm = new cfa_RNTR();
  } else if(optim == "BFGS") {
    algorithm = new cfa_BFGS();
  } else if(optim == "L-BFGS") {
    algorithm = new cfa_LBFGS();
  } else {

    Rcpp::stop("Available optimization rutines for cfa: \n 'gradient', 'BFGS', 'L-BFGS', and 'newtonTR'.");

  }

  return algorithm;

}
