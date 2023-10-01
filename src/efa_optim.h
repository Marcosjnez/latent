// Optimization algorithms for factor analysis

class efa_optim {

public:

  virtual efa_NTR optim(arguments_efa x,
                        efa_manifold *manifold,
                        efa_criterion *criterion) = 0;

};

// Gradient descent:

class efa_RGD:public efa_optim {

public:

  efa_NTR optim(arguments_efa x,
                efa_manifold *manifold,
                efa_criterion *criterion) {

    return gd(x, manifold, criterion);

  }

};

// L-BFGS algorithm:

class efa_LBFGS:public efa_optim {

public:

  efa_NTR optim(arguments_efa x,
                efa_manifold *manifold,
                efa_criterion *criterion) {

    return lbfgs(x, manifold, criterion);

  }

};

// Newton Trust-Region:

class efa_RNTR:public efa_optim {

public:

  efa_NTR optim(arguments_efa x,
                efa_manifold *manifold,
                efa_criterion *criterion) {

    return ntr(x, manifold, criterion);

  }

};

efa_optim* choose_efa_optim(std::string optim) {

  efa_optim* algorithm;
  if(optim == "gradient") {
    algorithm = new efa_RGD();
  } else if(optim == "L-BFGS") {
    algorithm = new efa_LBFGS();
  } else if(optim == "newtonTR") {
    algorithm = new efa_RNTR();
  } else {

    Rcpp::stop("Available optimization rutines for factor extraction: \n gradient, L-BFGS, newtonTR");

  }

  return algorithm;

}
