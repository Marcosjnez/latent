// Optimization algorithms for rotation

class rotation_optim {

public:

  virtual NTR optim(arguments_rotate x, rotation_manifold *manifold,
                     rotation_criterion *criterion) = 0;

};

// Riemannian gradient descent:

class RGD:public rotation_optim {

public:

  NTR optim(arguments_rotate x, rotation_manifold *manifold,
             rotation_criterion *criterion) {

    return gd(x, manifold, criterion);

  }

};

// Riemannian Newton Trust-Region:

class RNTR:public rotation_optim {

public:

  NTR optim(arguments_rotate x, rotation_manifold *manifold,
             rotation_criterion *criterion) {

    return ntr(x, manifold, criterion);

  }

};

// BFGS algorithm:

class BFGS:public rotation_optim {

public:

  NTR optim(arguments_rotate x, rotation_manifold *manifold,
            rotation_criterion *criterion) {

    return bfgs(x, manifold, criterion);

  }

};

// L-BFGS algorithm:

class LBFGS:public rotation_optim {

public:

  NTR optim(arguments_rotate x, rotation_manifold *manifold,
            rotation_criterion *criterion) {

    return lbfgs(x, manifold, criterion);

  }

};

// Select the optimization algorithm:

rotation_optim* choose_optim(std::string optim) {

  rotation_optim* algorithm;
  if(optim == "gradient") {
    algorithm = new RGD();
  } else if(optim == "newtonTR") {
    algorithm = new RNTR();
  } else if(optim == "BFGS") {
    algorithm = new BFGS();
  } else if(optim == "L-BFGS") {
    algorithm = new LBFGS();
  } else {

    Rcpp::stop("Available optimization rutines for rotation: \n 'gradient', 'BFGS', 'L-BFGS', 'newtonTR'. The default method is 'newtonTR'.");

  }

  return algorithm;

}
