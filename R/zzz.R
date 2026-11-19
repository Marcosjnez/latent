# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 06/10/2025

.onAttach <- function(libname, pkgname) {

  path <- system.file("extdata", "standard.png", package = pkgname)
  if (nzchar(path)) {
    if (interactive() && requireNamespace("png", quietly = TRUE) && requireNamespace("grid", quietly = TRUE)) {
      img <- png::readPNG(path)
      grid::grid.raster(img)
    }
  }

  packageStartupMessage("\U0001f47b Welcome to latent!")

  version <- as.character(utils::packageVersion(pkgname))
  msg <- sprintf("
      __      __             __
     / /_____/ /____  ____  / /_
    / / __  / __/ _ \\/ __ \\/ __/
   / / /_/ / /_/  __/ / / / /_
  /_/\\____/\\__/\\___/_/ /_/\\__/   version %s

Type 'citation(\"%s\")' for citing this package in publications.
", version, pkgname)

  packageStartupMessage(msg)
  packageStartupMessage("Report bugs at m.j.jimenezhenriquez@vu.nl")
  packageStartupMessage("For tutorials, visit marcosjnez.github.io/latent/")
  # packageStartupMessage("🌈 Happy pride month!")

}

# Model builders unpack named lists with list2env(). Some names are also
# evaluated as columns in data.table or subset() expressions.
utils::globalVariables(c(
  ".", ".weights", "A_group", "alpha_group", "any_gaussian", "any_multinomial",
  "any_mvgaussian", "any_mvmultinomial", "B_group", "control_estimator",
  "cov_params", "cov_params_labels", "cov_patterns", "data", "data_param",
  "delta_group", "design", "eta_hat_list", "full2short", "gaussian",
  "gaussian_names", "group", "group1", "indep_pairs", "indep_pairs_list",
  "indicators_names", "init_logvar", "init_mean", "init_param", "init_sigma",
  "init_var", "joint_orderings", "kappa_group", "lambda_group",
  "latent_cov_group", "latent_means_group", "loginter_names", "M_group", "mat",
  "means_group", "means_params", "means_params_labels", "meanshat_group",
  "model_group", "multinomial", "multinomial_factor_lengths",
  "multinomial_factor_levels", "multinomial_log_rows", "multinomial_names",
  "multinomial_prob_rows", "multinomial_rownames", "mvgaussian",
  "mvgaussian_names", "mvmultinomial", "new_J", "new_K", "new_levels",
  "ngaussian", "njoints", "nmultinomial", "nmvgaussian", "nobs", "nobs_ij",
  "npairs", "nparam", "npatterns", "npossible_patterns", "ntrans", "nu_group",
  "NVCOV_cov", "NVCOV_means", "pairs", "param", "parameters",
  "parameters_labels", "pattern_means_group", "pattern_names",
  "pattern_S_group", "pattern_weights", "patterns_gaussian",
  "patterns_multinomial", "patterns_mvgaussian",
  "patterns_mvmultinomial_recoded", "pdesign", "pi_hat_list",
  "previous_models", "psi_group", "S_group", "S_matrix", "short2full",
  "sigma_names", "sqrtvars_group", "structural_group", "target", "target_psi",
  "target_theta", "tauhat_group", "taus_group", "taus_item", "theta_group",
  "trajectory_weights", "trans", "trans2param", "transparameters",
  "transparameters_labels", "ualpha_group", "ulambda_group", "upsi_group",
  "variable_type", "variables", "X_group", "X_matrix", "Xinv_group",
  "xpsi_group", "xtheta_group"
))
