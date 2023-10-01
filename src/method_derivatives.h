/*
 * Author: Marcos Jimenez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 18/03/2022
 *
 */

// #include "auxiliary_manifolds.h"
// #include "asymptotic_cov.h"


/*
 * Hessian matrix: d2f/(dtheta dtheta')
 */

/*
 * Derivatives for Rhat wrt Lambda, Phi, and Psi
 */

// arma::mat gLRhat(arma::mat Lambda, arma::mat Phi) {
//
//   int p = Lambda.n_rows;
//   int q = Lambda.n_cols;
//   arma::mat I(p, p, arma::fill::eye);
//   arma::mat LP = Lambda * Phi;
//   arma::mat g1 = arma::kron(LP, I);
//   arma::mat g21 = arma::kron(I, LP);
//   arma::mat g2 = g21 * dxt(p, q);
//   arma::mat g = g1 + g2;
//
//   return g;
//
// }
//
// arma::mat gPRhat(arma::mat Lambda, arma::mat Phi, arma::uvec indexes) {
//
//   int q = Phi.n_cols;
//   // arma::uvec indexes = trimatl_ind(arma::size(Phi), -1);
//   arma::mat g1 = arma::kron(Lambda, Lambda);
//   arma::mat g2 = g1 * dxt(q, q);
//   arma::mat g_temp = g1 + g2;
//   arma::mat g = g_temp.cols(indexes);
//
//   return g;
// }

arma::mat guRhat(int p) {

  arma::mat gu(p*p, p, arma::fill::zeros);

  for(int i=0; i < p; ++i) {

    int index = i*p + i;
    gu(index, i) = 1;

  }

  return gu;

}

/*
 * Derivatives for dwls:
 * d2f/(dtheta dtheta') and d2f/(dtheta ds)
 */

arma::mat hessian_Lambda(arma::mat S, arma::mat Lambda,
                         arma::mat Phi, arma::vec psi, arma::mat W) {

  int p = Lambda.n_rows;
  int q = Phi.n_rows;
  int pq = p*q;

  arma::vec w = arma::vectorise(W);
  arma::mat Rhat = Lambda * Phi * Lambda.t() + arma::diagmat(psi);
  arma::mat residuals = (S - Rhat) % W;
  arma::mat I(p, p, arma::fill::eye);
  arma::mat g0 = gLRhat(Lambda, Phi);
  g0.each_col() %= w;
  arma::mat h1 = 2*arma::kron((Lambda * Phi).t(), I) * g0;
  arma::mat h2 = -2*arma::kron(Phi, residuals);

  arma::mat hessian_Lambda = h1 + h2; // second derivatives for Lambda

  return hessian_Lambda;

}

arma::mat hessian_Phi(arma::mat Lambda, arma::mat Phi, arma::uvec indexes,
                      arma::mat W) {

  // arma::uvec indexes = trimatl_ind(arma::size(Phi), -1);
  arma::vec w = arma::vectorise(W);
  arma::mat LxL = 2*arma::kron(Lambda, Lambda);
  arma::mat LxLt = LxL.t();
  arma::mat g0 = gPRhat(Lambda, Phi, indexes);
  g0.each_col() %= w;

  arma::mat hessian_Phi = LxLt.rows(indexes) * g0;

  return hessian_Phi;

}

arma::mat hessian_psi(arma::vec psi, arma::mat W) {

  int p = psi.size();
  arma::mat I(p, p, arma::fill::zeros);
  I.diag() = W.diag();

  return I;

}

arma::mat hessian_LambdaPhi(arma::mat S, arma::mat Lambda,
                            arma::mat Phi, arma::vec psi,
                            arma::uvec indexes, arma::mat W) {

  int p = Lambda.n_rows;
  int q = Phi.n_rows;
  arma::vec w = arma::vectorise(W);

  arma::mat LambdaPhi = Lambda * Phi;
  arma::mat Rhat = LambdaPhi * Lambda.t() + arma::diagmat(psi);
  arma::mat residuals = W % (S - Rhat);
  arma::mat I1(p, p, arma::fill::eye);
  arma::mat I2(q, q, arma::fill::eye);

  // arma::uvec indexes = trimatl_ind(arma::size(Phi), -1);
  arma::mat g0 = gPRhat(Lambda, Phi, indexes);
  g0.each_col() %= w;
  arma::mat h1 = 2*arma::kron(LambdaPhi.t(), I1) * g0;
  arma::mat h21 = -2*arma::kron(I2, residuals * Lambda);
  arma::mat h22 = h21 + h21 * dxt(q, q);
  arma::mat h2 = h22.cols(indexes);
  arma::mat h = h1 + h2;

  return h;

}

arma::mat hessian_Lambdapsi(arma::mat Lambda, arma::mat Phi, arma::mat W) {

  int p = Lambda.n_rows;
  int q = Phi.n_rows;
  arma::rowvec w = arma::diagvec(W);

  arma::mat LP = 2*Lambda * Phi;
  arma::mat h_Lambdapsi = arma::diagmat(LP.col(0));

  for(int i=1; i < q; ++i) {

    arma::mat aug = arma::diagmat(LP.col(i));
    h_Lambdapsi = arma::join_cols(h_Lambdapsi, aug);

  }

  h_Lambdapsi.each_row() %= w;

  return h_Lambdapsi;

}

arma::mat hessian_Phipsi(arma::mat Lambda, arma::mat Phi, arma::uvec indexes) {

  int p = Lambda.n_rows;
  // arma::vec w = arma::vectorise(W);

  arma::mat gPR = gPRhat(Lambda, Phi, indexes);
  // gPR.each_col() %= w;
  arma::uvec indexes_psi(p);
  for(int i=0; i < p; ++i) indexes_psi[i] = i*p + i;
  arma::mat h_Phipsi = gPR.rows(indexes_psi);

  return h_Phipsi.t();

}

arma::mat hessian_uls(arma::mat S, arma::mat Lambda, arma::mat Phi,
                         std::string projection, arma::uvec indexes) {

  /*
   * Compute all the functions above and create the hessian matrix
   */

  /*
   * Second derivatives for Lambda
   */

  int p = Lambda.n_rows;
  int q = Phi.n_rows;
  int pq = p*q;

  arma::mat LambdaPhi = Lambda * Phi;
  arma::mat Rhat = Lambda * Phi * Lambda.t();
  Rhat.diag().ones();
  arma::mat residuals = S - Rhat;
  arma::mat I1(p, p, arma::fill::eye);
  arma::mat gLR = gLRhat(Lambda, Phi);
  arma::mat LambdaPhitxI1 = 2*arma::kron(LambdaPhi.t(), I1);
  arma::mat h1 = LambdaPhitxI1 * gLR;
  arma::mat h2 = -2*arma::kron(Phi, residuals);

  arma::mat h_Lambda = h1 + h2;

  arma::mat gPR;
  // arma::uvec indexes;
  arma::mat h_Phi;

  if(projection == "oblq" | projection == "poblq") {

    /*
     * Second derivatives for Phi
     */

    // indexes = trimatl_ind(arma::size(Phi), -1);
    arma::mat LxL = 2*arma::kron(Lambda, Lambda);
    arma::mat LxLt = LxL.t();
    gPR = gPRhat(Lambda, Phi, indexes);
    h_Phi = LxLt.rows(indexes) * gPR;

  }
  /*
   * Second derivatives for psi
   */

  arma::mat h_psi = I1;
  arma::mat h_LambdaPhi;

  if(projection == "oblq" | projection == "poblq") {

    /*
     * Second derivatives for Lambda and Phi
     */

    arma::mat I2(q, q, arma::fill::eye);

    h1 = LambdaPhitxI1 * gPR;
    arma::mat h21 = -2*arma::kron(I2, residuals * Lambda);
    arma::mat h22 = h21 + h21 * dxt(q, q);
    h2 = h22.cols(indexes);
    h_LambdaPhi = h1 + h2;

  }

  /*
   * Second derivatives for Lambda and Psi
   */

  arma::mat LP = 2*LambdaPhi;
  arma::mat h_Lambdapsi = arma::diagmat(LP.col(0));

  for(int i=1; i < q; ++i) {

    arma::mat aug = arma::diagmat(LP.col(i));
    h_Lambdapsi = arma::join_cols(h_Lambdapsi, aug);

  }

  arma::mat h_Phipsi;

  if(projection == "oblq" | projection == "poblq") {

    /*
     * Second derivatives for Phi and Psi
     */

    arma::uvec indexes2(p);
    for(int i=0; i < p; ++i) indexes2[i] = i*p + i;
    h_Phipsi = gPR.rows(indexes2).t();

  }

  arma::mat hessian;

  if(projection == "oblq" | projection == "poblq") {

    /*
     * Join all the derivatives such that
     * h_Lambda         h_LambdaPhi   h_Lambdapsi
     * h_LambdaPhi.t()  h_Phi         h_Phipsi
     * h_Lambdapsi.t()  h_Phipsi.t()  h_psi
     */

    int q_cor;
    if(projection == "oblq") {
      q_cor = q*(q-1)/2;
    }
    if(projection == "poblq") {
      q_cor = indexes.size();
    }

    // insert columnwise:
    arma::mat hessian1 = h_Lambda; // pq x pq
    hessian1.insert_cols(pq, h_LambdaPhi); // pq x q_cor
    hessian1.insert_cols(pq + q_cor, h_Lambdapsi); // pq x p

    // insert columnwise:
    arma::mat hessian2 = h_LambdaPhi.t(); // q_cor x pq
    hessian2.insert_cols(pq, h_Phi); // q_cor x q_cor
    hessian2.insert_cols(pq + q_cor, h_Phipsi); // q_cor x p

    // insert columnwise:
    arma::mat hessian3 = h_Lambdapsi.t(); // p x pq
    hessian3.insert_cols(pq, h_Phipsi.t()); // p x pq
    hessian3.insert_cols(pq + q_cor, h_psi); // p x p

    // stack the blocks:
    hessian = arma::join_cols(hessian1, hessian2, hessian3);

  } else if(projection == "orth") {

    /*
     * Join all the derivatives such that
     * h_Lambda        h_Lambdapsi
     * h_Lambdapsi.t()    h_psi
     */

    // insert columnwise:
    arma::mat hessian1 = h_Lambda; // pq x pq
    hessian1.insert_cols(pq, h_Lambdapsi); // pq x p

    // insert columnwise:
    arma::mat hessian2 = h_Lambdapsi.t(); // p x pq
    hessian2.insert_cols(pq, h_psi); // p x p

    // stack the blocks:
    hessian = arma::join_cols(hessian1, hessian2);

  }

  return hessian;

}

/*
 * d2f/(dtheta ds) for uls
 */

arma::mat gLS_uls(arma::mat S, arma::mat Lambda, arma::mat Phi) {

  int p = Lambda.n_rows;
  arma::uvec indexes = trimatl_ind(arma::size(S), -1);

  arma::mat LambdaPhi = Lambda * Phi;
  arma::mat I(p, p, arma::fill::eye);
  arma::mat g1 = -2*arma::kron(LambdaPhi.t(), I);
  arma::mat g2 = g1 * dxt(p, p);
  arma::mat g_temp = g1 + g2;
  arma::mat g = g_temp.cols(indexes);

  return g;

}

arma::mat gPS_uls(arma::mat S, arma::mat Lambda, arma::mat Phi,
                     arma::uvec indexes1) {

  int p = S.n_rows;
  // arma::uvec indexes1 = trimatl_ind(arma::size(Phi), -1);
  arma::uvec indexes2 = trimatl_ind(arma::size(S), -1);

  arma::mat g1 = -2*arma::kron(Lambda.t(), Lambda.t());
  arma::mat g2 = g1 * dxt(p, p);
  arma::mat g_temp = g1 + g2;
  arma::mat g = g_temp(indexes1, indexes2);

  return g;

}

arma::mat gLPS_uls(arma::mat S, arma::mat Lambda, arma::mat Phi,
                      std::string projection, arma::uvec indexes1) {

  /*
   * Compute d2f/(dtheta ds)
   */

  int p = Lambda.n_rows;
  int q = Lambda.n_cols;

  // arma::uvec indexes1 = trimatl_ind(arma::size(Phi), -1);
  arma::uvec indexes2 = trimatl_ind(arma::size(S), -1);

  arma::mat LambdaPhi = Lambda * Phi;
  arma::mat I(p, p, arma::fill::eye);
  arma::mat g1 = -2*arma::kron(LambdaPhi.t(), I);
  arma::mat g2 = g1 * dxt(p, p);
  arma::mat g_temp = g1 + g2;
  arma::mat g = g_temp.cols(indexes2);

  arma::mat gd;
  int k;

  if(projection == "orth") {

    gd = g;
    k = p*q;

  } else if(projection == "oblq" | projection == "poblq") {

    arma::mat d1 = -2*arma::kron(Lambda.t(), Lambda.t());
    arma::mat d2 = d1 * dxt(p, p);
    arma::mat d_temp = d1 + d2;
    arma::mat d = d_temp(indexes1, indexes2);
    gd = arma::join_cols(g, d);

    if(projection == "oblq") k = p*q + q*(q-1)/2;
    if(projection == "poblq") k = p*q + indexes1.size();

  }

  /*
   * fill p rows with zeros
   */

  gd.insert_rows(k, p);

  return gd;

}

/*
 * Hessian for ml
 */

arma::mat hessian_ml(arma::mat S, arma::mat Lambda, arma::mat Phi,
                     std::string projection, arma::uvec indexes) {

  /*
   * Compute all the functions above and create the hessian matrix
   */

  /*
   * Second derivatives for Lambda
   */

  int p = Lambda.n_rows;
  int q = Phi.n_rows;
  int pp = p*p;
  int pq = p*q;
  arma::mat I1(p, p, arma::fill::eye);

  arma::mat LambdaPhi = Lambda * Phi;
  arma::mat Rhat = LambdaPhi * Lambda.t();
  Rhat.diag().ones();
  arma::mat Rhat_inv = arma::inv_sympd(Rhat);
  arma::mat residuals = Rhat - S;
  arma::mat Ri_res_Ri = 2*Rhat_inv * residuals * Rhat_inv;
  arma::mat h1 = arma::kron(Phi, Ri_res_Ri);

  arma::mat dRhat_dL = gLRhat(Lambda, Phi);
  arma::mat dRi_res_Ri_dRhat = 2*arma::kron(Rhat_inv, Rhat_inv) -
    arma::kron(Ri_res_Ri, Rhat_inv) - arma::kron(Rhat_inv, Ri_res_Ri);
  arma::mat dRi_res_Ri_dL = dRi_res_Ri_dRhat * dRhat_dL;
  arma::mat LPtxI = arma::kron(LambdaPhi.t(), I1);
  arma::mat h2 = LPtxI * dRi_res_Ri_dL;

  arma::mat h_Lambda = h1 + h2;

  arma::mat dRhat_dP;
  // arma::uvec indexes;
  arma::mat dRi_res_Ri_dP;
  arma::mat h_Phi;

  if(projection == "oblq" | projection == "poblq") {

    /*
     * Second derivatives for Phi
     */

    // indexes = trimatl_ind(arma::size(Phi), -1);
    arma::mat LxL = 2*arma::kron(Lambda, Lambda);
    arma::mat LxLt = LxL.t();
    dRhat_dP = gPRhat(Lambda, Phi, indexes);
    dRi_res_Ri_dP = dRi_res_Ri_dRhat * dRhat_dP;
    arma::mat LtxLt = arma::kron(Lambda.t(), Lambda.t());

    h_Phi = LtxLt.rows(indexes) * dRi_res_Ri_dP;

  }
  /*
   * Second derivatives for psi
   */

  arma::uvec indexes2 = arma::linspace<arma::uvec>(0, pp-1, p);
  arma::mat dRhat_du = guRhat(p);
  arma::mat dRi_res_Ri_du = dRi_res_Ri_dRhat * dRhat_du;
  arma::mat h_psi = 0.5*dRi_res_Ri_du.rows(indexes2);

  arma::mat h_LambdaPhi;

  if(projection == "oblq" | projection == "poblq") {

    /*
     * Second derivatives for Lambda and Phi
     */

    arma::mat I2(q, q, arma::fill::eye);

    arma::mat dxtPhi = dxt(q, q);
    h1 = LPtxI * dRi_res_Ri_dP;
    arma::mat h21 = arma::kron(I2, Ri_res_Ri * Lambda);
    arma::mat h22 = h21 * dxtPhi;
    h2 = h21 + h22;

    h_LambdaPhi = h1 + h2.cols(indexes);

  }

  /*
   * Second derivatives for Lambda and Psi
   */

  arma::mat h_Lambdapsi = LPtxI * dRi_res_Ri_du;

  arma::mat h_Phipsi;

  if(projection == "oblq" | projection == "poblq") {

    /*
     * Second derivatives for Phi and Psi
     */

    h_Phipsi = 0.5*dRi_res_Ri_dP.rows(indexes2).t();

  }

  arma::mat hessian;

  if(projection == "oblq" | projection == "poblq") {

    /*
     * Join all the derivatives such that
     * h_Lambda         h_LambdaPhi   h_Lambdapsi
     * h_LambdaPhi.t()  h_Phi         h_Phipsi
     * h_Lambdapsi.t()  h_Phipsi.t()  h_psi
     */

    int q_cor;
    if(projection == "oblq") q_cor = q*(q-1)/2;
    if(projection == "poblq") q_cor = indexes.size();

    // insert columnwise:
    arma::mat hessian1 = h_Lambda; // pq x pq
    hessian1.insert_cols(pq, h_LambdaPhi); // pq x q_cor
    hessian1.insert_cols(pq + q_cor, h_Lambdapsi); // pq x p

    // insert columnwise:
    arma::mat hessian2 = h_LambdaPhi.t(); // q_cor x pq
    hessian2.insert_cols(pq, h_Phi); // q_cor x q_cor
    hessian2.insert_cols(pq + q_cor, h_Phipsi); // q_cor x p

    // insert columnwise:
    arma::mat hessian3 = h_Lambdapsi.t(); // p x pq
    hessian3.insert_cols(pq, h_Phipsi.t()); // p x pq
    hessian3.insert_cols(pq + q_cor, h_psi); // p x p

    // stack the blocks:
    hessian = arma::join_cols(hessian1, hessian2, hessian3);

  } else if(projection == "orth") {

    /*
     * Join all the derivatives such that
     * h_Lambda        h_Lambdapsi
     * h_Lambdapsi.t()    h_psi
     */

    // insert columnwise:
    arma::mat hessian1 = h_Lambda; // pq x pq
    hessian1.insert_cols(pq, h_Lambdapsi); // pq x p

    // insert columnwise:
    arma::mat hessian2 = h_Lambdapsi.t(); // p x pq
    hessian2.insert_cols(pq, h_psi); // p x p

    // stack the blocks:
    hessian = arma::join_cols(hessian1, hessian2);

  }

  return hessian;

}

/*
 * d2f/(dtheta ds) for ml
 */

arma::mat gLPS_ml(arma::mat S, arma::mat Lambda, arma::mat Phi,
                  std::string projection, arma::uvec indexes1) {

  /*
   * Compute d2f/(dtheta ds)
   */

  int p = Lambda.n_rows;
  int q = Lambda.n_cols;
  int pp = p*p;

  // arma::uvec indexes1 = trimatl_ind(arma::size(Phi), -1);
  arma::uvec indexes2 = trimatl_ind(arma::size(S), -1);
  arma::uvec indexes3 = arma::linspace<arma::uvec>(0, pp-1, p);
  arma::mat I1(p, p, arma::fill::eye);

  /*
   * for Lambda
   */

  arma::mat LambdaPhi = Lambda * Phi;
  arma::mat Rhat = LambdaPhi * Lambda.t();
  Rhat.diag().ones();
  arma::mat Rhat_inv = arma::inv_sympd(Rhat);
  arma::mat dRi_res_Ri_dS = -2*arma::kron(Rhat_inv, Rhat_inv);
  arma::mat dxtS = dxt(p, p);
  dRi_res_Ri_dS += dRi_res_Ri_dS * dxtS;
  arma::mat g_temp = arma::kron(LambdaPhi.t(), I1) * dRi_res_Ri_dS;
  arma::mat g = g_temp.cols(indexes2);

  arma::mat gd;
  int k;

  /*
   * for Phi
   */

  if(projection == "orth") {

    gd = g;
    k = p*q;

  } else if(projection == "oblq" | projection == "poblq") {

    arma::mat d1 = arma::kron(Lambda.t(), Lambda.t());
    arma::mat d_temp = d1 * dRi_res_Ri_dS;

    arma::mat d = d_temp(indexes1, indexes2);
    gd = arma::join_cols(g, d);
    if(projection == "oblq") k = p*q + q*(q-1)/2;
    if(projection == "poblq") k = p*q + indexes1.size();

  }

  /*
   * for psi
   */

  arma::mat h = 0.5*dRi_res_Ri_dS(indexes3, indexes2);
  arma::mat gdh = arma::join_cols(gd, h);

  return gdh;

}


// B matrix for sandwich estimator

arma::mat B(arma::mat S, arma::mat Lambda, arma::mat Phi,
            arma::uvec loblq_indexes,
            Rcpp::Nullable<arma::mat> nullable_X = R_NilValue,
            std::string method = "uls", std::string projection = "oblq",
            std::string type = "continuous", double eta = 1) {

  /*
   * B matrix from the sandwich covariance matrix estimator
   */

  arma::uvec indexes = trimatl_ind(arma::size(S), -1);
  arma::mat gLPS;

  if(method == "uls") {

    gLPS = gLPS_uls(S, Lambda, Phi, projection, loblq_indexes);

  } else if(method == "ml") {

    gLPS = gLPS_ml(S, Lambda, Phi, projection, loblq_indexes);

  }

  arma::mat asymp = asymp_cov(S, nullable_X, eta, type);
  arma::mat asymp2 = asymp(indexes, indexes);

  arma::mat B = gLPS * asymp2 * gLPS.t();

  return B;

}



