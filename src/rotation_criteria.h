/*
 * Author: Marcos Jimenez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 28/05/2022
 *
 */

// #include "auxiliary_criteria.h"

// Criteria for rotation

class rotation_criterion {

public:

  virtual void F(arguments_rotate& x) = 0;

  virtual void gLP(arguments_rotate& x) = 0;

  virtual void hLP(arguments_rotate& x) = 0;

  virtual void dgLP(arguments_rotate& x) = 0;

};

/*
 * none
 */

class none: public rotation_criterion {

public:

  void F(arguments_rotate& x) {}

  void gLP(arguments_rotate& x) {}

  void hLP(arguments_rotate& x) {}

  void dgLP(arguments_rotate& x) {}

};

/*
 * Crawford-Ferguson family
 */

class cf: public rotation_criterion {

public:

  void F(arguments_rotate& x){

    x.L2 = x.L % x.L;
    x.L2N = x.L2 * x.N;
    x.ML2 = x.M * x.L2;
    double ff1 = (1-x.k[0]) * arma::accu(x.L2 % x.L2N) / 4;
    double ff2 = x.k[0] * arma::accu(x.L2 % x.ML2) / 4;

    x.f = ff1 + ff2;

  }

  void gLP(arguments_rotate& x){

    x.f1 = (1-x.k[0]) * x.L % x.L2N;
    x.f2 = x.k[0] * x.L % x.ML2;
    x.gL = x.f1 + x.f2;

  }

  void hLP(arguments_rotate& x) {

    arma::colvec L_vector = arma::vectorise(x.L);

    arma::mat Ip(x.p, x.p, arma::fill::eye);
    arma::mat c1 = arma::kron(x.N.t(), Ip) * arma::diagmat(2*L_vector);
    c1.each_col() %= L_vector;
    arma::mat c2 = arma::diagmat(arma::vectorise(x.L2 * x.N));
    arma::mat gf1 = (1-x.k[0]) * (c1 + c2);

    arma::mat Iq(x.q, x.q, arma::fill::eye);
    arma::mat c3 = arma::kron(Iq, x.M) * arma::diagmat(2*L_vector);
    c3.each_col() %= L_vector;
    arma::mat c4 = arma::diagmat(arma::vectorise(x.M * x.L2));
    arma::mat gf2 = x.k[0] * (c3 + c4);
    x.hL = gf1 + gf2;

  }

  void dgLP(arguments_rotate& x) {

    arma::mat dL2 = 2 * x.dL % x.L;
    arma::mat df1 = (1-x.k[0]) * x.dL % x.L2N + (1-x.k[0]) * x.L % (dL2 * x.N);
    arma::mat df2 = x.k[0] * x.dL % x.ML2 + x.k[0] * x.L % (x.M * dL2);

    x.dgL = df1 + df2;

  }

};

/*
 * Varimax
 */

class varimax: public rotation_criterion {

public:

  void F(arguments_rotate& x) {

    x.L2 = x.L % x.L;
    x.HL2 = x.H * x.L2;

    x.f = -arma::trace(x.HL2.t() * x.HL2) / 4;

  }

  void gLP(arguments_rotate& x) {

    x.gL = -x.L % x.HL2;

  }

  void hLP(arguments_rotate& x) {

    arma::mat Iq(x.q, x.q, arma::fill::eye);
    arma::mat c1 = arma::diagmat(arma::vectorise(x.HL2));
    arma::colvec L_vector = arma::vectorise(x.L);
    arma::mat c2 = arma::kron(Iq, x.H) * arma::diagmat(2*L_vector);
    c2.each_col() %= L_vector;
    x.hL = -(c1 + c2);

  }

  void dgLP(arguments_rotate& x) {

    arma::mat dL2 = 2 * x.dL % x.L;
    x.dgL = -x.dL % x.HL2 - x.L % (x.H * dL2);

  }

};

/*
 * Varimin
 */

class varimin: public rotation_criterion {

public:

  void F(arguments_rotate& x) {

    x.L2 = x.L % x.L;
    x.HL2 = x.H * x.L2;

    x.f = arma::trace(x.HL2.t() * x.HL2) / 4;

  }

  void gLP(arguments_rotate& x) {

    x.gL = x.L % x.HL2;

  }

  void hLP(arguments_rotate& x) {

    arma::mat Iq(x.q, x.q, arma::fill::eye);
    arma::mat c1 = arma::diagmat(arma::vectorise(x.HL2));
    arma::colvec L_vector = arma::vectorise(x.L);
    arma::mat c2 = arma::kron(Iq, x.H) * arma::diagmat(2*L_vector);
    c2.each_col() %= L_vector;
    x.hL = c1 + c2;

  }

  void dgLP(arguments_rotate& x) {

    arma::mat dL2 = 2 * x.dL % x.L;
    x.dgL = x.dL % x.HL2 + x.L % (x.H * dL2);

  }

};

/*
 * Oblimin
 */

class oblimin: public rotation_criterion {

public:

  void F(arguments_rotate& x) {

    x.L2 = x.L % x.L;
    x.IgCL2N = x.I_gamma_C * x.L2 * x.N;

    x.f = arma::accu(x.L2 % x.IgCL2N) / 4;

  }

  void gLP(arguments_rotate& x) {

    x.gL = x.L % x.IgCL2N;

  }

  void hLP(arguments_rotate& x) {

    arma::mat c1 = arma::diagmat(arma::vectorise(x.IgCL2N));
    arma::mat diagL = arma::diagmat(arma::vectorise(x.L));
    arma::mat diag2L = 2*diagL;
    arma::mat c2 = diagL * arma::kron(x.N, x.I_gamma_C) * diag2L;
    x.hL = c1 + c2;

  }

  void dgLP(arguments_rotate& x) {

    x.dgL = x.dL % x.IgCL2N + x.L % (x.I_gamma_C * (2*x.dL % x.L) * x.N);

  }

};

/*
 * Geomin
 */

class geomin: public rotation_criterion {

public:

  void F(arguments_rotate& x) {

    x.q2 = 2/(x.q + 0.0);
    x.L2 = x.L % x.L;
    x.L2 += x.epsilon[0];
    x.term = arma::trunc_exp(arma::sum(arma::trunc_log(x.L2), 1) / x.q);

    x.f = arma::accu(x.term);

  }

  void gLP(arguments_rotate& x) {

    x.LoL2 = x.L / x.L2;
    x.gL = x.LoL2 * x.q2;
    x.gL.each_col() %= x.term;

  }

  void hLP(arguments_rotate& x) {

    arma::mat cx = x.q2 * x.LoL2;

    arma::mat c1 = x.q2*(arma::vectorise(x.L2) - arma::vectorise(2*x.L % x.L)) /
      arma::vectorise(x.L2 % x.L2);
    arma::mat gcx = arma::diagmat(c1);
    arma::mat c2 = (1/x.L2) % (2*x.L) / x.q;
    c2.each_col() %= x.term;
    arma::mat gterm = cbind_diag(c2);
    arma::mat v = gterm.t() * cx;
    x.hL = gcx;
    arma::mat term2 = x.term;
    for(int i=0; i < (x.q-1); ++i) term2 = arma::join_cols(term2, x.term);
    x.hL.each_col() %= term2;
    x.hL += kdiag(v);

  }

  void dgLP(arguments_rotate& x) {

    arma::mat c1 = (x.epsilon[0] - x.L % x.L) / (x.L2 % x.L2) % x.dL;
    c1.each_col() %= x.term;
    arma::mat c2 = x.LoL2;
    arma::vec term2 = x.q2 * x.term % arma::sum(x.LoL2 % x.dL, 1);
    c2.each_col() %= term2;

    x.dgL = x.q2 * (c1 + c2);

  }

};

/*
 * Target
 */

class target: public rotation_criterion {

public:
  void F(arguments_rotate& x) {

    x.f1 = x.Weight % (x.L - x.Target);

    x.f = 0.5*arma::accu(x.f1 % x.f1);

  }

  void gLP(arguments_rotate& x) {

    x.gL = x.Weight % x.f1;

  }

  void hLP(arguments_rotate& x) {

    arma::mat W2 = x.Weight % x.Weight;
    x.hL = arma::diagmat(arma::vectorise(W2));

  }

  void dgLP(arguments_rotate& x) {

    x.dgL = x.Weight2 % x.dL;

  }

};

/*
 * xTarget
 */

class xtarget: public rotation_criterion {

public:

  void F(arguments_rotate& x) {

    x.f1 = x.Weight % (x.L - x.Target);
    x.f2 = x.Phi_Weight % (x.Phi - x.Phi_Target);

    x.f = 0.5*arma::accu(x.f1 % x.f1) + 0.25*x.w*arma::accu(x.f2 % x.f2);

  }

  void gLP(arguments_rotate& x) {

    x.gL = x.Weight % x.f1;
    x.gP = x.w * x.Phi_Weight % x.f2;

  }

  void hLP(arguments_rotate& x) {

    x.hL = arma::diagmat(arma::vectorise(x.Weight2));

    arma::mat wPW2 = arma::diagmat(arma::vectorise(x.w*x.Phi_Weight2));
    arma::mat hP = wPW2 + wPW2 * dxt(x.q, x.q);
    arma::uvec lower_indices = arma::trimatl_ind(arma::size(x.Phi), -1);
    x.hP = hP;
    // x.hP = hP.cols(lower_indices);


  }

  void dgLP(arguments_rotate& x) {

    x.dgL = x.Weight2 % x.dL;
    x.dgP = x.w * x.Phi_Weight2 % x.dP;

  }

};

/*
 * equavar (orthogonal)
 */

class equavar: public rotation_criterion {

public:

  void F(arguments_rotate& x) {

    double q = x.q + 0.0;
    x.var = arma::diagvec(x.L.t() * x.L);
    x.varq = arma::pow(x.var, 1/q);
    x.prodvarq = arma::prod(x.varq);

    x.f = x.prodvarq * -1;

  }

  void gLP(arguments_rotate& x) {

    double q = x.q + 0.0;
    x.dvarqdL = 2*x.L;
    x.dvarqdL.each_row() %= arma::pow(x.var, 1/q-1)/q;
    x.gL = x.dvarqdL;
    x.gL.each_row() %= (1/x.varq);
    x.gL *= x.prodvarq * -1;

  }

  void hLP(arguments_rotate& x) {

    Rcpp::stop("Standard errors not implement yet for equavar");

  }

  void dgLP(arguments_rotate& x) {

    double q = x.q + 0.0;
    arma::mat dvardL = 2*x.L;
    arma::mat DdvardL = 2*x.dL;
    arma::rowvec Dvar = arma::diagvec(x.dL.t() * x.L + x.L.t() * x.dL);
    arma::rowvec Dvarq = arma::pow(x.var, 1/q-1)/q % Dvar;
    arma::mat dvarqdL = dvardL;
    dvarqdL.each_row() %= arma::pow(x.var, 1/q-1)/q;
    arma::mat DdvarqdL1 = DdvardL;
    DdvarqdL1.each_row() %= arma::pow(x.var, 1/q-1)/q;
    arma::mat DdvarqdL2 = dvardL;
    DdvarqdL2.each_row() %= arma::pow(x.var, 1/q-2)/q * (1/q-1) % Dvar;
    arma::mat DdvarqdL = DdvarqdL1 + DdvarqdL2;
    double Dprodvarq = arma::accu(x.prodvarq / x.varq % Dvarq);

    arma::mat dgL1 = Dprodvarq * dvarqdL + x.prodvarq * DdvarqdL;
    dgL1.each_row() %= x.varq;
    arma::mat dgL2 = x.prodvarq * dvarqdL;
    dgL2.each_row() %= Dvarq;
    x.dgL = dgL1 - dgL2;
    x.dgL.each_row() /= (x.varq % x.varq) * -1;

  }

};

/*
 * Linear CLF
 */

class clf_l: public rotation_criterion {

public:

  void F(arguments_rotate& x) {

    double b = 1 / (2*x.clf_epsilon[0]);
    double a = x.clf_epsilon[0] - b*x.clf_epsilon[0]*x.clf_epsilon[0];

    arma::mat absL = arma::abs(x.L);
    x.lower = arma::find(absL <= x.clf_epsilon[0]);
    x.larger = arma::find(absL > x.clf_epsilon[0]);
    double f1 = arma::accu(a + b*absL.elem(x.lower) % absL.elem(x.lower));
    double f2 = arma::accu(absL.elem(x.larger));

    x.f = f1 + f2;

  }

  void gLP(arguments_rotate& x) {

    double b = 1 / (2*x.clf_epsilon[0]);
    x.gL.set_size(arma::size(x.L));
    x.gL.elem(x.lower) = 2*b*x.L.elem(x.lower);
    x.gL.elem(x.larger) = arma::sign(x.L.elem(x.larger));

  }

  void hLP(arguments_rotate& x) {

    int p = x.L.n_rows;
    int q = x.L.n_cols;
    arma::mat hL(p*q, p*q, arma::fill::zeros);

    arma::mat absL = arma::abs(x.L);
    x.lower = arma::find(absL <= x.clf_epsilon[0]);
    double d = 1 / x.clf_epsilon[0];

    // for(int i=0; i < x.lower.size(); ++i) {
    //   // x.hL(arma::span(x.lower[i]), arma::span(x.lower[i])) = d;
    //   x.hL(x.lower[i], x.lower[i]) = d;
    // }

  }

  void dgLP(arguments_rotate& x) {

    double b = 1 / (2*x.clf_epsilon[0]);
    x.dgL.set_size(arma::size(x.L));
    x.dgL.zeros();
    x.dgL.elem(x.lower) = 2*b*x.dL.elem(x.lower);

  }

};

/*
 * Invariance
 */

class invar: public rotation_criterion {

public:

  void F(arguments_rotate& x) {

    x.f = 0;
    int k = x.indexes1.n_rows;
    for(int i=0; i < k; ++i) {

      double dif = x.L(x.indexes1(i, 0), x.indexes1(i, 1)) -
        x.L(x.indexes2(i, 0), x.indexes2(i, 1));
      x.f += dif*dif;

    }

  }

  void gLP(arguments_rotate& x) {

    int k = x.indexes1.n_rows;
    x.gL.set_size(arma::size(x.L));
    x.gL.zeros();

    for(int i=0; i < k; ++i) {

      double dif = x.L(x.indexes1(i, 0), x.indexes1(i, 1)) -
        x.L(x.indexes2(i, 0), x.indexes2(i, 1));
      x.gL(x.indexes1(i, 0), x.indexes1(i, 1)) = 2*dif;
      x.gL(x.indexes2(i, 0), x.indexes2(i, 1)) = -x.gL(x.indexes1(i, 0), x.indexes1(i, 1));

    }

  }

  void hLP(arguments_rotate& x) {

    Rcpp::warning("Standard errors not implement yet for invar");

  }

  void dgLP(arguments_rotate& x) {

    x.dgL.set_size(arma::size(x.L));
    x.dgL.zeros();
    int k = x.indexes1.n_rows;

    for(int i=0; i < k; ++i) {

      x.dgL(x.indexes1(i, 0), x.indexes1(i, 1)) = 2*(x.dL(x.indexes1(i, 0), x.indexes1(i, 1))) -
        x.dL(x.indexes2(i, 0), x.indexes2(i, 1));
      x.dgL(x.indexes2(i, 0), x.indexes2(i, 1)) = -x.dgL(x.indexes1(i, 0), x.indexes1(i, 1));

    }

  }

};

/*
 * Repeated Crawford-Ferguson family
 */

class rep_cf: public rotation_criterion {

public:

  void F(arguments_rotate& x){

    arma::uvec rows_indexes = x.rows_list[x.i];
    arma::uvec cols_indexes = x.cols_list[x.i];

    x.Li[x.i] = x.L(rows_indexes, cols_indexes);
    x.Li2[x.i] = x.Li[x.i] % x.Li[x.i];
    x.Ni[x.i] = x.N(cols_indexes, cols_indexes);
    x.Mi[x.i] = x.M(rows_indexes, rows_indexes);
    double ff1 = (1-x.k[x.i]) * arma::accu(x.Li2[x.i] % (x.Li2[x.i] * x.Ni[x.i])) / 4;
    double ff2 = x.k[x.i] * arma::accu(x.Li2[x.i] % (x.Mi[x.i] * x.Li2[x.i])) / 4;

    x.f += (ff1 + ff2) * x.block_weights[x.i];

  }

  void gLP(arguments_rotate& x) {

    arma::uvec rows_indexes = x.rows_list[x.i];
    arma::uvec cols_indexes = x.cols_list[x.i];

    arma::mat f1 = (1-x.k[x.i]) * x.Li[x.i] % (x.Li2[x.i] * x.Ni[x.i]);
    arma::mat f2 = x.k[x.i] * x.Li[x.i] % (x.Mi[x.i] * x.Li2[x.i]);
    x.gL(rows_indexes, cols_indexes) += (f1 + f2) * x.block_weights[x.i];

  }

  void hLP(arguments_rotate& x) {

    Rcpp::stop("Standard errors not implement yet for block criteria");

    arma::colvec L_vector = arma::vectorise(x.Li[x.i]);

    arma::mat Ip(x.p, x.p, arma::fill::eye);
    arma::mat c1 = arma::kron(x.Ni[x.i].t(), Ip) * arma::diagmat(2*L_vector);
    c1.each_col() %= L_vector;
    arma::mat c2 = arma::diagmat(arma::vectorise(x.Li2[x.i] * x.Ni[x.i]));
    arma::mat gf1 = (1-x.k[x.i]) * (c1 + c2);

    arma::mat Iq(x.qi[x.i], x.qi[x.i], arma::fill::eye);
    arma::mat c3 = arma::kron(Iq, x.M) * arma::diagmat(2*L_vector);
    c3.each_col() %= L_vector;
    arma::mat c4 = arma::diagmat(arma::vectorise(x.M * x.Li2[x.i]));
    arma::mat gf2 = x.k[x.i] * (c3 + c4);

    x.hL += (gf1 + gf2) * x.block_weights[x.i];

  }

  void dgLP(arguments_rotate& x) {

    arma::uvec rows_indexes = x.rows_list[x.i];
    arma::uvec cols_indexes = x.cols_list[x.i];

    arma::mat dLi = x.dL(rows_indexes, cols_indexes);
    arma::mat dLi2 = 2 * dLi % x.Li[x.i];
    arma::mat df1 = (1-x.k[x.i]) * dLi % (x.Li2[x.i] * x.Ni[x.i]) +
      (1-x.k[x.i]) * x.Li[x.i] % (dLi2 * x.Ni[x.i]);
    arma::mat df2 = x.k[x.i] * dLi % (x.Mi[x.i] * x.Li2[x.i]) +
      x.k[x.i] * x.Li[x.i] % (x.Mi[x.i] * dLi2);

    x.dgL(rows_indexes, cols_indexes) += (df1 + df2) * x.block_weights[x.i];

  }

};

/*
 * Repeated Varimax
 */

class rep_varimax: public rotation_criterion {

public:

  void F(arguments_rotate& x) {

    arma::uvec rows_indexes = x.rows_list[x.i];
    arma::uvec cols_indexes = x.cols_list[x.i];

    x.Li[x.i] = x.L(rows_indexes, cols_indexes);
    x.Li2[x.i] = x.Li[x.i] % x.Li[x.i];
    x.Hi[x.i] = x.H(rows_indexes, rows_indexes);
    x.HLi2[x.i] = x.Hi[x.i] * x.Li2[x.i];

    x.f -= arma::trace(x.HLi2[x.i].t() * x.HLi2[x.i]) / 4 * x.block_weights[x.i];

  }

  void gLP(arguments_rotate& x) {

    arma::uvec rows_indexes = x.rows_list[x.i];
    arma::uvec cols_indexes = x.cols_list[x.i];

    x.gL(rows_indexes, cols_indexes) -= (x.Li[x.i] % x.HLi2[x.i]) * x.block_weights[x.i];

  }

  void hLP(arguments_rotate& x) {

    Rcpp::stop("Standard errors not implement yet for block criteria");

    arma::mat Iq(x.qi[x.i], x.qi[x.i], arma::fill::eye);
    arma::mat c1 = arma::diagmat(arma::vectorise(x.HL2i[x.i]));
    arma::colvec L_vector = arma::vectorise(x.Li[x.i]);
    arma::mat c2 = arma::kron(Iq, x.Hi[x.i]) * arma::diagmat(2*L_vector);
    c2.each_col() %= L_vector;

    x.hL -= (c1 + c2) * x.block_weights[x.i];

  }

  void dgLP(arguments_rotate& x) {

    arma::uvec rows_indexes = x.rows_list[x.i];
    arma::uvec cols_indexes = x.cols_list[x.i];

    arma::mat dLi = x.dL(rows_indexes, cols_indexes);
    arma::mat dLi2 = 2 * dLi % x.Li[x.i];
    x.dgL(rows_indexes, cols_indexes) -= (dLi % x.HLi2[x.i] - x.Li[x.i] % (x.Hi[x.i] * dLi2)) * x.block_weights[x.i];

  }

};

/*
 * Repeated Varimin
 */

class rep_varimin: public rotation_criterion {

public:

  void F(arguments_rotate& x) {

    arma::uvec rows_indexes = x.rows_list[x.i];
    arma::uvec cols_indexes = x.cols_list[x.i];

    x.Li[x.i] = x.L(rows_indexes, cols_indexes);
    x.Li2[x.i] = x.Li[x.i] % x.Li[x.i];
    x.Hi[x.i] = x.H(rows_indexes, rows_indexes);
    x.HLi2[x.i] = x.Hi[x.i] * x.Li2[x.i];

    x.f += arma::trace(x.HLi2[x.i].t() * x.HLi2[x.i]) / 4 * x.block_weights[x.i];

  }

  void gLP(arguments_rotate& x) {

    arma::uvec rows_indexes = x.rows_list[x.i];
    arma::uvec cols_indexes = x.cols_list[x.i];

    x.gL(rows_indexes, cols_indexes) += (x.Li[x.i] % x.HLi2[x.i]) * x.block_weights[x.i];

  }

  void hLP(arguments_rotate& x) {

    Rcpp::stop("Standard errors not implement yet for block criteria");

    arma::mat Iq(x.qi[x.i], x.qi[x.i], arma::fill::eye);
    arma::mat c1 = arma::diagmat(arma::vectorise(x.HL2i[x.i]));
    arma::colvec L_vector = arma::vectorise(x.Li[x.i]);
    arma::mat c2 = arma::kron(Iq, x.Hi[x.i]) * arma::diagmat(2*L_vector);
    c2.each_col() %= L_vector;

    x.hL += (c1 + c2) * x.block_weights[x.i];

  }

  void dgLP(arguments_rotate& x) {

    arma::uvec rows_indexes = x.rows_list[x.i];
    arma::uvec cols_indexes = x.cols_list[x.i];

    arma::mat dLi = x.dL(rows_indexes, cols_indexes);
    arma::mat dLi2 = 2 * dLi % x.Li[x.i];
    x.dgL(rows_indexes, cols_indexes) -= (dLi % x.HLi2[x.i] - x.Li[x.i] % (x.Hi[x.i] * dLi2)) * x.block_weights[x.i];

  }

};

/*
 * Repeated Oblimin
 */

class rep_oblimin: public rotation_criterion {

public:

  void F(arguments_rotate& x) {

    arma::uvec rows_indexes = x.rows_list[x.i];
    arma::uvec cols_indexes = x.cols_list[x.i];

    x.Li[x.i] = x.L(rows_indexes, cols_indexes);
    x.Ni[x.i] = x.N(cols_indexes, cols_indexes);
    x.Li2[x.i] = x.Li[x.i] % x.Li[x.i];
    x.IgCL2Ni[x.i] = x.I_gamma_Ci[x.i] * x.Li2[x.i] * x.Ni[x.i];

    x.f += arma::trace(x.Li2[x.i].t() * x.IgCL2Ni[x.i]) / 4 * x.block_weights[x.i];

  }

  void gLP(arguments_rotate& x) {

    arma::uvec rows_indexes = x.rows_list[x.i];
    arma::uvec cols_indexes = x.cols_list[x.i];

    x.gL(rows_indexes, cols_indexes) += (x.Li[x.i] % x.IgCL2Ni[x.i]) * x.block_weights[x.i];

  }

  void hLP(arguments_rotate& x) {

    Rcpp::stop("Standard errors not implement yet for block criteria");

    arma::mat c1 = arma::diagmat(arma::vectorise(x.IgCL2Ni[x.i]));
    arma::mat diagL = arma::diagmat(arma::vectorise(x.Li[x.i]));
    arma::mat diag2L = 2*diagL;
    arma::mat c2 = diagL * arma::kron(x.Ni[x.i], x.I_gamma_Ci[x.i]) * diag2L;
    x.hL += (c1 + c2) * x.block_weights[x.i];

  }

  void dgLP(arguments_rotate& x) {

    arma::uvec rows_indexes = x.rows_list[x.i];
    arma::uvec cols_indexes = x.cols_list[x.i];

    arma::mat dLi = x.dL(rows_indexes, cols_indexes);
    x.dgL(rows_indexes, cols_indexes) += (dLi % x.IgCL2Ni[x.i] +
      x.Li[x.i] % (x.I_gamma_Ci[x.i] * (2*dLi % x.Li[x.i]) * x.Ni[x.i])) * x.block_weights[x.i];

  }

};

/*
 * Repeated Geomin
 */

class rep_geomin: public rotation_criterion {

public:

  void F(arguments_rotate& x) {

    arma::uvec rows_indexes = x.rows_list[x.i];
    arma::uvec cols_indexes = x.cols_list[x.i];

    int q = x.qi[x.i];
    x.q2 = 2/(q + 0.0);

    x.Li[x.i] = x.L(rows_indexes, cols_indexes);
    x.Li2[x.i] = x.Li[x.i] % x.Li[x.i];
    x.Li2[x.i] += x.epsilon[x.i];
    x.termi[x.i] = arma::exp(arma::sum(arma::log(x.Li2[x.i]), 1) / q);

    x.f += arma::accu(x.termi[x.i]) * x.block_weights[x.i];

  }

  void gLP(arguments_rotate& x) {

    arma::uvec rows_indexes = x.rows_list[x.i];
    arma::uvec cols_indexes = x.cols_list[x.i];

    int q = x.qi[x.i];
    x.q2 = 2/(q + 0.0);

    x.LoLi2[x.i] = x.Li[x.i] / x.Li2[x.i];
    x.termi[x.i] = arma::exp(arma::sum(arma::log(x.Li2[x.i]), 1) / q);
    arma::mat gLi = x.LoLi2[x.i] * x.q2;
    gLi.each_col() %= x.termi[x.i];
    x.gL(rows_indexes, cols_indexes) += gLi * x.block_weights[x.i];

  }

  void hLP(arguments_rotate& x) {

    Rcpp::stop("Standard errors not implement yet for block criteria");

    arma::mat cx = x.q2 * x.LoLi2[x.i];
    arma::mat c1 = x.q2*(arma::vectorise(x.Li2[x.i]) -
      arma::vectorise(2*x.Li[x.i] % x.Li[x.i])) /
      arma::vectorise(x.Li2[x.i] % x.Li2[x.i]);
    arma::mat gcx = arma::diagmat(c1);
    arma::mat c2 = (1/x.Li2[x.i]) % (2*x.Li[x.i]) / x.qi[x.i];
    c2.each_col() %= x.termi[x.i];
    arma::mat gterm = cbind_diag(c2);
    arma::mat v = gterm.t() * cx;
    arma::mat hL = gcx;
    arma::mat term2 = x.termi[x.i];
    for(int i=0; i < (x.qi[x.i]-1); ++i) term2 = arma::join_cols(term2, x.termi[x.i]);
    hL.each_col() %= term2;
    hL += kdiag(v);

    x.hL += hL * x.block_weights[x.i];

  }

  void dgLP(arguments_rotate& x) {

    arma::uvec rows_indexes = x.rows_list[x.i];
    arma::uvec cols_indexes = x.cols_list[x.i];

    int q = x.qi[x.i];
    x.q2 = 2/(q + 0.0);

    arma::mat dLi = x.dL(rows_indexes, cols_indexes);
    arma::mat c1 = (x.epsilon[x.i] - x.Li[x.i] % x.Li[x.i]) / (x.Li2[x.i] % x.Li2[x.i]) % dLi;
    c1.each_col() %= x.termi[x.i];
    arma::mat c2 = x.LoLi2[x.i];
    arma::vec termi2 = x.q2 * x.termi[x.i] % arma::sum(x.LoLi2[x.i] % dLi, 1);
    c2.each_col() %= termi2;

    x.dgL(rows_indexes, cols_indexes) += x.q2 * (c1 + c2) * x.block_weights[x.i];

  }

};

/*
 * Repeated Target
 */

class rep_target: public rotation_criterion {

public:

  void F(arguments_rotate& x) {

    arma::uvec rows_indexes = x.rows_list[x.i];
    arma::uvec cols_indexes = x.cols_list[x.i];

    x.Li[x.i] = x.L(rows_indexes, cols_indexes);
    x.f1i[x.i] = x.Weight % (x.Li[x.i] - x.Target);
    x.f += 0.5*arma::accu(x.f1i[x.i] % x.f1i[x.i]) * x.block_weights[x.i];

  }

  void gLP(arguments_rotate& x) {

    arma::uvec rows_indexes = x.rows_list[x.i];
    arma::uvec cols_indexes = x.cols_list[x.i];

    x.gL(rows_indexes, cols_indexes) += x.Weight % x.f1i[x.i] * x.block_weights[x.i];

  }

  void hLP(arguments_rotate& x) {

    Rcpp::stop("Standard errors not implement yet for block criteria");

    arma::mat W2 = x.Weight % x.Weight;
    x.hL += arma::diagmat(arma::vectorise(W2));

  }

  void dgLP(arguments_rotate& x) {

    arma::uvec rows_indexes = x.rows_list[x.i];
    arma::uvec cols_indexes = x.cols_list[x.i];

    arma::mat dLi = x.dL(rows_indexes, cols_indexes);

    x.dgL(rows_indexes, cols_indexes) += x.Weight % x.Weight % dLi * x.block_weights[x.i];

  }

};

/*
 * Repeated xTarget
 */

class rep_xtarget: public rotation_criterion {

public:

  void F(arguments_rotate& x) {

    arma::uvec rows_indexes = x.rows_list[x.i];
    arma::uvec cols_indexes = x.cols_list[x.i];

    x.Li[x.i] = x.L(rows_indexes, cols_indexes);
    x.Phii[x.i] = x.Phi(cols_indexes, cols_indexes);

    x.f1i[x.i] = x.Weight % (x.Li[x.i] - x.Target);
    x.f2i[x.i] = x.Phi_Weight % (x.Phii[x.i] - x.Phi_Target);

    x.f += (0.5*arma::accu(x.f1i[x.i] % x.f1i[x.i]) +
      0.25*x.w*arma::accu(x.f2i[x.i] % x.f2i[x.i])) * x.block_weights[x.i];

  }

  void gLP(arguments_rotate& x) {

    arma::uvec rows_indexes = x.rows_list[x.i];
    arma::uvec cols_indexes = x.cols_list[x.i];

    x.gL(rows_indexes, cols_indexes) += x.Weight % x.f1i[x.i] * x.block_weights[x.i];
    x.gP(cols_indexes, cols_indexes) += x.w * x.Phi_Weight % x.f2i[x.i] * x.block_weights[x.i];

  }

  void hLP(arguments_rotate& x) {

    Rcpp::stop("Standard errors not implement yet for block criteria");

    x.hL += arma::diagmat(arma::vectorise(x.Weight2));

    arma::mat Phi_t = dxt(x.qi[x.i], x.qi[x.i]);
    arma::mat diag_PW2 = arma::diagmat(arma::vectorise(x.Phi_Weight2));
    x.hP += diag_PW2 + diag_PW2 * Phi_t; // w¿?

  }

  void dgLP(arguments_rotate& x) {

    arma::uvec rows_indexes = x.rows_list[x.i];
    arma::uvec cols_indexes = x.cols_list[x.i];

    arma::mat dPi = x.dP(cols_indexes, cols_indexes);

    x.dgL(rows_indexes, cols_indexes) += x.Weight2 % x.dL * x.block_weights[x.i];
    x.dgP(cols_indexes, cols_indexes) += x.w * x.Phi_Weight2 % dPi * x.block_weights[x.i];

  }

};

/*
 * Repeated equavar (orthogonal)
 */

class rep_equavar: public rotation_criterion {

public:

  void F(arguments_rotate& x) {

    arma::uvec rows_indexes = x.rows_list[x.i];
    arma::uvec cols_indexes = x.cols_list[x.i];

    x.Li[x.i] = x.L(rows_indexes, cols_indexes);

    double q = x.qi[x.i] + 0.0;
    x.vari[x.i] = arma::diagvec(x.Li[x.i].t() * x.Li[x.i]);
    x.varqi[x.i] = arma::pow(x.vari[x.i], 1/q);
    x.prodvarqi[x.i] = arma::prod(x.varqi[x.i]);

    x.f += x.prodvarqi[x.i] * -1 * x.block_weights[x.i];

  }

  void gLP(arguments_rotate& x) {

    arma::uvec rows_indexes = x.rows_list[x.i];
    arma::uvec cols_indexes = x.cols_list[x.i];

    double q = x.qi[x.i] + 0.0;
    x.dvarqdLi[x.i] = 2*x.Li[x.i];
    x.dvarqdLi[x.i].each_row() %= arma::pow(x.vari[x.i], 1/q-1)/q;
    arma::mat gL = x.dvarqdLi[x.i];
    gL.each_row() %= (1/x.varqi[x.i]);
    gL *= x.prodvarqi[x.i];
    x.gL(rows_indexes, cols_indexes) += gL * -1 * x.block_weights[x.i];

  }

  void hLP(arguments_rotate& x) {

    Rcpp::stop("Standard errors not implement yet for block criteria");

  }

  void dgLP(arguments_rotate& x) {

    arma::uvec rows_indexes = x.rows_list[x.i];
    arma::uvec cols_indexes = x.cols_list[x.i];

    arma::mat dLi = x.dL(rows_indexes, cols_indexes);

    double q = x.qi[x.i] + 0.0;
    arma::mat dvardL = 2*x.Li[x.i];
    arma::mat DdvardL = 2*dLi;
    arma::rowvec Dvar = arma::diagvec(dLi.t() * x.Li[x.i] + x.Li[x.i].t() * dLi);
    arma::rowvec Dvarq = arma::pow(x.vari[x.i], 1/q-1)/q % Dvar;
    arma::mat dvarqdL = dvardL;
    dvarqdL.each_row() %= arma::pow(x.vari[x.i], 1/q-1)/q;
    arma::mat DdvarqdL1 = DdvardL;
    DdvarqdL1.each_row() %= arma::pow(x.vari[x.i], 1/q-1)/q;
    arma::mat DdvarqdL2 = dvardL;
    DdvarqdL2.each_row() %= arma::pow(x.vari[x.i], 1/q-2)/q * (1/q-1) % Dvar;
    arma::mat DdvarqdL = DdvarqdL1 + DdvarqdL2;
    double Dprodvarq = arma::accu(x.prodvarqi[x.i] / x.varqi[x.i] % Dvarq);

    arma::mat dgL1 = Dprodvarq * dvarqdL + x.prodvarqi[x.i] * DdvarqdL;
    dgL1.each_row() %= x.varqi[x.i];
    arma::mat dgL2 = x.prodvarqi[x.i] * dvarqdL;
    dgL2.each_row() %= Dvarq;
    arma::mat dgL = dgL1 - dgL2;
    dgL /= (x.varqi[x.i] % x.varqi[x.i]);

    x.dgL(rows_indexes, cols_indexes) += dgL * -1 * x.block_weights[x.i];

  }

};

/*
 * Repeated Linear CLF
 */

class rep_clf_l: public rotation_criterion {

public:

  void F(arguments_rotate& x) {

    arma::uvec rows_indexes = x.rows_list[x.i];
    arma::uvec cols_indexes = x.cols_list[x.i];

    x.Li[x.i] = x.L(rows_indexes, cols_indexes);

    double b = 1 / (2*x.clf_epsilon[x.i]);
    double a = x.clf_epsilon[x.i] - b*x.clf_epsilon[x.i]*x.clf_epsilon[x.i];

    arma::mat absLi = arma::abs(x.Li[x.i]);
    x.loweri[x.i] = arma::find(absLi <= x.clf_epsilon[x.i]);
    x.largeri[x.i] = arma::find(absLi > x.clf_epsilon[x.i]);
    double f1 = arma::accu(a + b*absLi.elem(x.loweri[x.i]) % absLi.elem(x.loweri[x.i]));
    double f2 = arma::accu(absLi.elem(x.largeri[x.i]));

    x.f += (f1 + f2) * x.block_weights[x.i];

  }

  void gLP(arguments_rotate& x) {

    arma::uvec rows_indexes = x.rows_list[x.i];
    arma::uvec cols_indexes = x.cols_list[x.i];

    double b = 1 / (2*x.clf_epsilon[x.i]);
    arma::mat gLi(arma::size(x.Li[x.i]));
    gLi.elem(x.loweri[x.i]) = 2*b*x.Li[x.i].elem(x.loweri[x.i]);
    gLi.elem(x.largeri[x.i]) = arma::sign(x.Li[x.i].elem(x.largeri[x.i]));

    x.gL(rows_indexes, cols_indexes) += gLi * x.block_weights[x.i];

  }

  void hLP(arguments_rotate& x) {

    Rcpp::stop("Standard errors not implement yet for block criteria");

  }

  void dgLP(arguments_rotate& x) {

    arma::uvec rows_indexes = x.rows_list[x.i];
    arma::uvec cols_indexes = x.cols_list[x.i];

    arma::mat dLi = x.dL(rows_indexes, cols_indexes);

    double b = 1 / (2*x.clf_epsilon[x.i]);
    arma::mat dgLi(arma::size(x.Li[x.i]), arma::fill::zeros);
    dgLi.elem(x.loweri[x.i]) = 2*b*dLi.elem(x.loweri[x.i]);

    x.dgL(rows_indexes, cols_indexes) += dgLi * x.block_weights[x.i];

  }

};

/*
 * Repeated Invariance
 */

class rep_invar: public rotation_criterion {

public:

  void F(arguments_rotate& x) {

    arma::uvec rows_indexes = x.rows_list[x.i];
    arma::uvec cols_indexes = x.cols_list[x.i];

    x.Li[x.i] = x.L(rows_indexes, cols_indexes);

    int k = x.indexes1.n_rows;
    for(int i=0; i < k; ++i) {

      double dif = x.Li[x.i](x.indexes1(i, 0), x.indexes1(i, 1)) -
        x.Li[x.i](x.indexes2(i, 0), x.indexes2(i, 1));
      x.f += dif*dif*x.block_weights[x.i];

    }

  }

  void gLP(arguments_rotate& x) {

    arma::uvec rows_indexes = x.rows_list[x.i];
    arma::uvec cols_indexes = x.cols_list[x.i];

    x.Li[x.i] = x.L(rows_indexes, cols_indexes);

    int k = x.indexes1.n_rows;
    arma::mat gLi(arma::size(x.Li[x.i]), arma::fill::zeros);

    for(int i=0; i < k; ++i) {

      double dif = x.Li[x.i](x.indexes1(i, 0), x.indexes1(i, 1)) -
        x.Li[x.i](x.indexes2(i, 0), x.indexes2(i, 1));
      gLi(x.indexes1(i, 0), x.indexes1(i, 1)) = 2*dif;
      gLi(x.indexes2(i, 0), x.indexes2(i, 1)) = -gLi(x.indexes1(i, 0), x.indexes1(i, 1));

    }

    x.gL(rows_indexes, cols_indexes) += gLi * x.block_weights[x.i];

  }

  void hLP(arguments_rotate& x) {

    Rcpp::stop("Standard errors not implement yet for block criteria");

  }

  void dgLP(arguments_rotate& x) {

    arma::uvec rows_indexes = x.rows_list[x.i];
    arma::uvec cols_indexes = x.cols_list[x.i];

    arma::mat dLi = x.dL(rows_indexes, cols_indexes);

    int k = x.indexes1.n_rows;
    arma::mat dgLi(arma::size(x.Li[x.i]), arma::fill::zeros);

    for(int i=0; i < k; ++i) {

      dgLi(x.indexes1(i, 0), x.indexes1(i, 1)) = 2*(dLi(x.indexes1(i, 0), x.indexes1(i, 1))) -
        dLi(x.indexes2(i, 0), x.indexes2(i, 1));
      dgLi(x.indexes2(i, 0), x.indexes2(i, 1)) = -dgLi(x.indexes1(i, 0), x.indexes1(i, 1));

    }

    x.dgL(rows_indexes, cols_indexes) += dgLi * x.block_weights[x.i];

  }

};

/*
 * Mixed criteria (multiple criteria simultaneously allowed)
 */

// Choose rotation criteria for mixed criteria:

rotation_criterion* choose_rep_criterion(std::string rotation, std::string projection) {

  rotation_criterion* criterion;

  if (rotation == "target") {

    criterion = new rep_target();

  }  else if(rotation == "xtarget") {

    if(projection == "orth") {
      criterion = new rep_target();
    } else {
      criterion = new rep_xtarget();
    }

  } else if(rotation == "cf") {

    criterion = new rep_cf();

  } else if(rotation == "oblimin") {

    criterion = new rep_oblimin();

  } else if(rotation == "geomin") {

    criterion = new rep_geomin();

  } else if(rotation == "varimax") {

    criterion = new rep_varimax();

  } else if(rotation == "varimin") {

    criterion = new rep_varimin();

  } else if(rotation == "equavar") {

    criterion = new rep_equavar();

  } else if(rotation == "clfl") {

    criterion = new rep_clf_l();

  } else if(rotation == "invar") {

    criterion = new rep_invar();

  } else if(rotation == "none") {

    criterion = new none();

  } else {

    Rcpp::stop("Available rotations: \n cf, oblimin, geomin, varimax, varimin, target, xtarget, equavar, clfl, invar");

  }

  return criterion;

}

class mixed: public rotation_criterion {

public:

  void F(arguments_rotate& x) {

    rotation_criterion* criterion;
    x.f = 0;

    for(int i=0; i < x.n_blocks; ++i) {

      x.i = i;
      // int report = x.cols_list[x.i].size();
      // Rcpp::Rcout << report << std::endl;

      criterion = choose_rep_criterion(x.rotations[i], x.projection);
      criterion->F(x);

    }

  }

  void gLP(arguments_rotate& x) {

    rotation_criterion* criterion;
    x.gL.set_size(arma::size(x.L));
    x.gL.zeros();
    x.gP.set_size(x.q, x.q);
    x.gP.zeros();

    for(int i=0; i < x.n_blocks; ++i) {

      x.i = i;
      criterion = choose_rep_criterion(x.rotations[i], x.projection);
      criterion->gLP(x);

    }

  }

  void hLP(arguments_rotate& x) {

    rotation_criterion* criterion;
    int pq = x.p*x.q;
    int qq = x.q*x.q;
    x.hL.set_size(pq, pq);
    x.hL.zeros();
    x.hP.set_size(qq, qq);
    x.hP.zeros();

    for(int i=0; i < x.n_blocks; ++i) {

      x.i = i;
      criterion = choose_rep_criterion(x.rotations[i], x.projection);
      criterion->hLP(x);

    }

  }

  void dgLP(arguments_rotate& x) {

    rotation_criterion* criterion;
    x.dgL.set_size(arma::size(x.L));
    x.dgL.zeros();
    x.dgP.set_size(x.q, x.q);
    x.dgP.zeros();

    for(int i=0; i < x.n_blocks; ++i) {

      x.i = i;
      criterion = choose_rep_criterion(x.rotations[i], x.projection);
      criterion->dgLP(x);

    }

  }

};

// Choose the rotation criteria:

rotation_criterion* choose_criterion(std::vector<std::string> rotations, std::string projection,
                                 std::vector<arma::uvec> cols_list) {

  rotation_criterion *criterion;

  if(!cols_list.empty()) {

    // Rcpp::stop("Mixed rotation criteria not supported yet");
    criterion = new mixed();

  } else if (rotations[0] == "target") {

    criterion = new target();

  }  else if(rotations[0] == "xtarget") {

    if(projection == "orth") {
      criterion = new target();
    } else {
      criterion = new xtarget();
    }

  } else if(rotations[0] == "cf") {

    criterion = new cf();

  } else if(rotations[0] == "oblimin") {

    criterion = new oblimin();

  } else if(rotations[0] == "geomin") {

    criterion = new geomin();

  } else if(rotations[0] == "varimax") {

    criterion = new varimax();

  } else if(rotations[0] == "varimin") {

    criterion = new varimin();

  } else if(rotations[0] == "equavar") {

    criterion = new equavar();

  } else if(rotations[0] == "clfl") {

    criterion = new clf_l();

  } else if(rotations[0] == "invar") {

    criterion = new invar();

  } else if(rotations[0] == "none") {

  } else {

    Rcpp::stop("Available rotations: \n cf, oblimin, geomin, varimax, varimin, target, xtarget, equavar, clfl, invar");

  }

  return criterion;

}
