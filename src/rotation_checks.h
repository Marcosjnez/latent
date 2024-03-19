/*
 * Author: Marcos Jimenez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 18/03/2022
 *
 */

// #include "structures.h"
// #include "auxiliary_checks.h"
// #include "auxiliary_manifolds.h"

void check_rotate(arguments_rotate& x, int random_starts, int cores) {

  if (x.nullable_PhiTarget.isNotNull()) {
    x.Phi_Target = Rcpp::as<arma::mat>(x.nullable_PhiTarget);
  }

  // Create a list of column indexes for each block of factors:

  if(x.nullable_blocks.isNotNull()) {

    // If the main input was a list of indexes for each block of factors,
    // adjust the indexes to start from 0:

    x.blocks = Rcpp::as<std::vector<std::vector<arma::uvec>>>(x.nullable_blocks);;
    x.rows_list = x.blocks[0];
    x.cols_list = x.blocks[1];
    arma::uvec v_rows = list_to_vector(x.rows_list);
    arma::uvec v_cols = list_to_vector(x.cols_list);

    if(v_rows.max() > x.p || v_rows.min() < 1) {
      Rcpp::stop("blocks must contain values between 1 and the number of variables");
    }
    if(v_cols.max() > x.q || v_cols.min() < 1) {
      Rcpp::stop("blocks must contain values between 1 and the number of factors");
    }

    // Number of blocks:
    x.n_blocks = x.rows_list.size();
    if(x.n_blocks != x.cols_list.size()) {
      Rcpp::stop("blocks must be a list containing two lists of the same size: one for the items and another for the factors");
    }

    for(int i=0; i < x.n_blocks; ++i) {

      x.rows_list[i] -= 1;
      x.cols_list[i] -= 1;

    }

    // Weights for each block of factors:

    if(x.nullable_block_weights.isNull()) {

      // If not specified, set the blocks to 1:
      x.block_weights.set_size(x.n_blocks);
      x.block_weights.ones();

    } else {

      x.block_weights = Rcpp::as<arma::vec>(x.nullable_block_weights);

      if(x.block_weights.size() != x.n_blocks) {
        Rcpp::stop("The number of weights must equal the number of blocks");
      }

    }

    // For mixed criteria, specify the size of the std::vector objects:
    x.qi.resize(x.n_blocks), x.pi.resize(x.n_blocks), x.Li.resize(x.n_blocks), x.Li2.resize(x.n_blocks),
    x.Mi.resize(x.n_blocks), x.Ni.resize(x.n_blocks), x.HLi2.resize(x.n_blocks), x.LoLi2.resize(x.n_blocks),
    x.termi.resize(x.n_blocks), x.IgCL2Ni.resize(x.n_blocks), x.f1i.resize(x.n_blocks),
    x.f2i.resize(x.n_blocks), x.I_gamma_Ci.resize(x.n_blocks),
    x.V.resize(x.n_blocks), x.vari.resize(x.n_blocks), x.varqi.resize(x.n_blocks),
    x.prodvari.resize(x.n_blocks), x.prodvarqi.resize(x.n_blocks),
    x.dvarqdLi.resize(x.n_blocks), x.dmudLi.resize(x.n_blocks), x.dc2dLi.resize(x.n_blocks),
    x.mui.resize(x.n_blocks), x.c2i.resize(x.n_blocks), x.dc2dmui.resize(x.n_blocks),
    x.prodc2i.resize(x.n_blocks), x.Hi.resize(x.n_blocks), x.dxtLi.resize(x.n_blocks),
    x.Ii.resize(x.n_blocks), x.dc2dLi.resize(x.n_blocks), x.dmudPi.resize(x.n_blocks),
    x.dc2dPi.resize(x.n_blocks), x.LtLxIi.resize(x.n_blocks), x.expmmui.resize(x.n_blocks),
    x.Phii.resize(x.n_blocks), x.dxtPi.resize(x.n_blocks), x.HL2i.resize(x.n_blocks);
    x.largeri.resize(x.n_blocks); x.loweri.resize(x.n_blocks);

    // Number of factors in each block (x.qi):
    for(int i=0; i < x.n_blocks; ++i) {
      x.qi[i] = x.cols_list[i].size();
      x.pi[i] = x.rows_list[i].size();
    }

    // Number of specified rotation criteria:
    x.n_rotations = x.rotations.size();

    if(x.n_rotations == 1) {

      // Resize the vector of rotations of length one to match the length of blocks:
      x.rotations.resize(x.n_blocks, x.rotations[0]);

    } else {

      if(x.n_rotations != x.n_blocks) Rcpp::stop("The number of rotation criteria and blocks of factors does not match. \n Provide either one rotation criteria for all the blocks or one for each block");

    }

  } else {

    if(x.n_rotations > 1) Rcpp::stop("Multiple rotation criteria specified but no block of factors was specified \n Provide a block for each rotation criteria");

  }

  // Check inputs for partially oblique projection:

  if(x.projection == "poblq") {

    if(x.nullable_oblq_factors.isNull()) {

      Rcpp::stop("Please, provide a vector with the number of factors in each oblique factor in the oblq_factors argument");

    } else {

      // Create indexes to fix to zero the entries corresponding to the nonconstrained correlations:

      x.oblq_indexes = Rcpp::as<arma::uvec>(x.nullable_oblq_factors);
      if(arma::accu(x.oblq_indexes) > x.q) Rcpp::stop("To many factors declared in oblq_factors");
      x.list_oblq_indexes = vector_to_list2(x.oblq_indexes);

      arma::mat X(x.q, x.q, arma::fill::ones);
      arma::mat Q = zeros(X, x.list_oblq_indexes);
      // Select the indexes for the duplicated oblique entries
      x.oblq_indexes = arma::find(Q == 0);
      // Select the indexes for the lower diagonal (nonducaplicated) orthogonal entries
      arma::uvec set_to_zero = arma::trimatu_ind(arma::size(X), 0);
      Q.elem(set_to_zero).zeros();
      x.orth_indexes = arma::find(Q == 1);
      // Select the indexes for the lower diagonal (nonducaplicated) oblique entries
      Q.elem(set_to_zero).ones();
      x.loblq_indexes = arma::find(Q == 0);

      if(!x.Phi_Target.is_empty()) {
        // Select the indexes for the duplicated oblique entries
        arma::mat Phi_Target = x.Phi_Target;
        Phi_Target.diag() += 10;
        x.oblq_indexes = arma::find(Phi_Target == 1);
        Phi_Target.elem(set_to_zero).ones();
        x.orth_indexes = arma::find(Phi_Target == 0);
        Q.elem(set_to_zero).zeros();
        x.loblq_indexes = arma::find(Phi_Target == 1);
      }

      // These indexes will be used to set to zero the oblique entries in the
      // A matrix in poblq projection and to extract the rotation constraints

    }

  }

  // Check and build defaults for each rotation criteria:

  std::vector<std::string> all_rotations = {"cf", "oblimin", "geomin", "target",
                                            "xtarget", "varimax", "varimin",
                                            "equavar", "simplix", "clfl", "invar",
                                            "geomin", "none"};

  // Check for invalid rotations:
  for (auto i: x.rotations) {
    if (std::find(all_rotations.begin(), all_rotations.end(), i) == all_rotations.end()) {
      Rcpp::stop("Available rotations: \n cf, oblimin, geomin, varimax, varimin, target, xtarget, equavar, simplix, clfl, invar");
    }
  }

  // TARGET ROTATION
  if(std::find(x.rotations.begin(), x.rotations.end(), "target") != x.rotations.end()) {

    if (x.nullable_Target.isNotNull()) {
      x.Target = Rcpp::as<arma::mat>(x.nullable_Target);
    } else {
      Rcpp::stop("Provide a Target for target rotation");
    }
    if (x.nullable_Weight.isNotNull()) {
      x.Weight = Rcpp::as<arma::mat>(x.nullable_Weight);
    } else {
      // If missing, Weight = 1 - Target (Partially specified target rotation)
      x.Weight = 1 - x.Target;
    }

    // Check the dimensions of the target:
    arma::mat dimcheck(x.p, x.q);
    if(x.n_blocks > 1) {
      // Find the position of the block corresponding to the target criteria and
      // resize dimcheck according to the number of factors specified in the block:
      int pos = std::find(x.rotations.begin(), x.rotations.end(), "target") - x.rotations.begin();
      dimcheck.set_size(x.pi[pos], x.qi[pos]);
    }

    if(arma::size(x.Target) != arma::size(dimcheck) ||
       arma::size(x.Weight) != arma::size(dimcheck)) {

      Rcpp::stop("Incompatible Target and Weight dimensions");

    }

    x.Weight2 = x.Weight % x.Weight;

  }

  // XTARGET ROTATION
  if(std::find(x.rotations.begin(), x.rotations.end(), "xtarget") != x.rotations.end()) {

    if(x.w < 0) Rcpp::stop("w must be nonnegative");

    if (x.nullable_Target.isNotNull()) {
      x.Target = Rcpp::as<arma::mat>(x.nullable_Target);
    } else {
      Rcpp::stop("Provide a Target for xtarget rotation");
    }
    if (x.nullable_PhiTarget.isNotNull()) {
      x.Phi_Target = Rcpp::as<arma::mat>(x.nullable_PhiTarget);
    } else {
      Rcpp::stop("Provide a PhiTarget for xtarget rotation");
    }
    if (x.nullable_Weight.isNotNull()) {
      x.Weight = Rcpp::as<arma::mat>(x.nullable_Weight);
    } else {
      // If missing, Weight = 1 - Target (Partially specified target rotation)
      x.Weight = 1-x.Target;
    }
    if (x.nullable_PhiWeight.isNotNull()) {
      x.Phi_Weight = Rcpp::as<arma::mat>(x.nullable_PhiWeight);
    } else {
      // If missing, PhiWeight = 1 - PhiTarget (Partially specified target rotation)
      x.Phi_Weight = 1-x.Phi_Target;
    }

    // Check the dimensions of the Target and PhiTarget
    arma::mat dimcheck1(x.p, x.q);
    arma::mat dimcheck2(x.q, x.q);
    if(x.n_blocks > 1) {
      // Find the position of the block corresponding to the target criteria and
      // resize dimcheck according to the number of factors specified in the target:
      int pos = std::find(x.rotations.begin(), x.rotations.end(), "xtarget") - x.rotations.begin();
      dimcheck1.set_size(x.pi[pos], x.qi[pos]);
      dimcheck2.set_size(x.qi[pos], x.qi[pos]);
    }

    if(arma::size(x.Target) != arma::size(dimcheck1) ||
       arma::size(x.Weight) != arma::size(dimcheck1) ||
       arma::size(x.Phi_Target) != arma::size(dimcheck2) ||
       arma::size(x.Phi_Weight) != arma::size(dimcheck2)) {

      Rcpp::stop("Incompatible dimensions between Target, PhiTarget, Weight or PhiWeight");

    }

    x.Phi_Weight.diag().zeros(); // This should have no effect, actually
    x.Weight2 = x.Weight % x.Weight;
    x.Phi_Weight2 = x.Phi_Weight % x.Phi_Weight;

  }

  // Linear CLF ROTATION
  if(std::find(x.rotations.begin(), x.rotations.end(), "clfl") != x.rotations.end()) {

    if(x.clf_epsilon.min() < 0) {

      Rcpp::stop("clf_epsilon must be nonnegative");

    }

    // Resize the constant vector so that it has the same length as the number of blocks.
    // Put each element of the original vector in the position corresponding to the CLF block:
    x.clf_epsilon = new_k(x.rotations, "clfl", x.clf_epsilon);

  }

  // GEOMIN ROTATION
  if(std::find(x.rotations.begin(), x.rotations.end(), "geomin") != x.rotations.end()) {

    if(x.epsilon.min() <= 0) {

      Rcpp::stop("epsilon must be positive");

    }

    // Resize the constant vector so that it has the same length as the number of blocks.
    // Put each element of the original vector in the position corresponding to the geomin block:
    x.epsilon = new_k(x.rotations, "geomin", x.epsilon);

  }

  // GEOMIN 2.0 ROTATION
  if(std::find(x.rotations.begin(), x.rotations.end(), "geomin2") != x.rotations.end()) {

    x.epsilones = random_orth(x.q, 1);

  }

  // OBLIMIN ROTATION
  if(std::find(x.rotations.begin(), x.rotations.end(), "oblimin") != x.rotations.end()) {

    if(x.gamma.min() < 0) {

      Rcpp::stop("gamma must be nonnegative");

    }

    x.N.set_size(x.q, x.q); x.N.ones();
    x.N.diag(0).zeros();
    // Resize the constant vector so that it has the same length as the number of blocks.
    // Put each element of the original vector in the position corresponding to the oblimin block:
    x.gamma = new_k(x.rotations, "oblimin", x.gamma);

    arma::mat I(x.p, x.p, arma::fill::eye), gamma_C(x.p, x.p, arma::fill::ones);

    if(!x.cols_list.empty()) {
      for(int i=0; i < x.n_blocks; ++i) { // if(!x.gamma[i].empty())
        arma::mat Ii(x.pi[i], x.pi[i], arma::fill::eye);
        arma::mat gamma_Ci(x.pi[i], x.pi[i], arma::fill::ones);
        gamma_Ci *= (x.gamma[i]/x.pi[i]);
        x.I_gamma_Ci[i] = (Ii - gamma_Ci);
      }
    } else {
      gamma_C *= (x.gamma[0]/x.p);
      x.I_gamma_C = (I - gamma_C);
    }

  }

  // CF ROTATION
  if(std::find(x.rotations.begin(), x.rotations.end(), "cf") != x.rotations.end()) {

    x.N.set_size(x.q, x.q); x.N.ones();
    x.N.diag(0).zeros();
    x.M.set_size(x.p, x.p); x.M.ones();
    x.M.diag(0).zeros();

    // Resize the constant vector so that it has the same length as the number of blocks.
    // Put each element of the original vector in the position corresponding to the CF block:
    x.k = new_k(x.rotations, "cf", x.k);

  }

  // VARIMAX ROTATION
  if(std::find(x.rotations.begin(), x.rotations.end(), "varimax") != x.rotations.end()) {

    if((x.projection == "oblq") | (x.projection == "poblq")) {
      Rcpp::warning("Usually, the varimax criterion does not converge with (partially) oblique projection. \n Consider using cf with k = 1/(number of items), which is equivalent to varimax for orthogonal rotation but also converges with (partially) oblique projection.");
    }

    // if(!x.blocks_list.empty()) {
    //   for(int i=0; i < x.n_blocks; ++i) {
    //     arma::vec v(x.p, arma::fill::ones);
    //     arma::mat I(x.p, x.p, arma::fill::eye);
    //     x.Hi[i] = I - v * v.t() / (x.p + 0.0); // Centering matrix
    //   }
    // } else {
      arma::vec v(x.p, arma::fill::ones);
      arma::mat I(x.p, x.p, arma::fill::eye);
      x.H = I - v * v.t() / (x.p + 0.0); // Centering matrix
    // }

  }

  // VARIMIN ROTATION
  if(std::find(x.rotations.begin(), x.rotations.end(), "varimin") != x.rotations.end()) {

    arma::vec v(x.p, arma::fill::ones);
    arma::mat I(x.p, x.p, arma::fill::eye);
    x.H = I - v * v.t() / (x.p + 0.0); // Centering matrix

  }

  // Check rotation control parameters:

  Rcpp::List rot_control;


  if (x.nullable_rot_control.isNotNull()) {

    rot_control = Rcpp::as<Rcpp::List>(x.nullable_rot_control);

  }

  if(rot_control.containsElementNamed("maxit")) {

    x.maxit = rot_control["maxit"];

  }

  if(rot_control.containsElementNamed("eps")) {

    x.eps = rot_control["eps"];

  }

  if(rot_control.containsElementNamed("optim")) {

    std::string optim = rot_control["optim"];
    x.optim = optim;

  }

  if(rot_control.containsElementNamed("indexes1")) {

    arma::mat indexes1 = rot_control["indexes1"];
    x.indexes1 = indexes1;

  }

  if(rot_control.containsElementNamed("indexes2")) {

    arma::mat indexes2 = rot_control["indexes2"];
    x.indexes2 = indexes2;

  }

  // INVARIANCE ROTATION
  // if(std::find(x.rotations.begin(), x.rotations.end(), "invar") != x.rotations.end()) {
  //
  //   if (x.nullable_indexes1.isNotNull()) {
  //     x.indexes1 = Rcpp::as<arma::mat>(x.nullable_indexes1);
  //   } else {
  //     Rcpp::stop("Provide a matrix with the indexes of the fixed loadings");
  //   }
  //
  //   if (x.nullable_indexes2.isNotNull()) {
  //     x.indexes2 = Rcpp::as<arma::mat>(x.nullable_indexes2);
  //   } else {
  //     Rcpp::stop("Provide a matrix with the indexes of the fixed loadings");
  //   }
  //
  // }

  // Check parallelization setup:

  if(random_starts < 1) Rcpp::stop("random_starts must be a positive integer");
  if(cores < 1) Rcpp::stop("The number of cores must be a positive integer");

}
