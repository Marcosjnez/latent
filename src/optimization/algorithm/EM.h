/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 13/07/2026
 */

// Expectation Maximization:

// optim_result em(arguments_optim x,
//                 std::vector<transformations*>& xtransforms,
//                 std::vector<manifolds*>& xmanifolds,
//                 std::vector<estimators*>& xestimators) {
//
//   product_transform* final_transform;
//   product_estimator* final_estimator;
//
//   final_transform->transform(x, xtransforms);
//
//   x.convergence = false;
//   arma::vec loglik(x.maxit); loglik[0] = x.loglik;
//   x.iterations = 0L;
//
//   do {
//
//     ++x.iterations;
//
//     final_estimator->param(x, xestimators);
//     final_estimator->F(x, xestimators); // Store f for outcomes()
//     final_estimator->E(x, xestimators);
//     // Rprintf("loglik = %g \n", x.loglik);
//     loglik[x.iterations] = x.loglik;
//
//     double diff = loglik[x.iterations]-loglik[x.iterations-1];
//     x.ng = sqrt(arma::accu(diff*diff));
//     if (x.ng < x.eps) {
//       x.convergence = true;
//       break;
//     }
//
//     // final_transform->M(x, xtransforms);
//
//     // Rf_error("617");
//
//   } while(x.iterations < x.maxit);
//
//   // Rf_error("652");
//   // final_transform->transform(x, xtransforms);
//   final_estimator->param(x, xestimators);
//   final_estimator->F(x, xestimators); // Store f for outcomes()
//   // x.f = -x.loglik;
//
//   optim_result result = std::make_tuple(x.parameters,
//                                         x.transparameters,
//                                         x.f,
//                                         x.iterations,
//                                         x.convergence,
//                                         x.ng, x.rg, x.g, x.dir);
//
//   return result;
//
// }

// class EM:public optim {
//
// public:
//
//   optim_result optimize(arguments_optim x,
//                         std::vector<transformations*>& xtransforms,
//                         std::vector<manifolds*>& xmanifolds,
//                         std::vector<estimators*>& xestimators) {
//
//     return em(x, xtransforms, xmanifolds, xestimators);
//
//   }
//
// };
