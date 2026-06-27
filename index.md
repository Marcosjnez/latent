# latent: Latent Class and Factor Analysis Models

Fit measurement models with discrete or continuous latent variables.

# Installation in Windows and Linux

``` R
devtools::install_github("marcosjnez/latent", force = TRUE)
```

# Installation in macOS

Install the `macrtools` R package of James Balamuta:

``` R
# install.packages("remotes")
remotes::install_github("coatless-mac/macrtools")
```

Install Command Line Tools and the R Compilation Toolchain (this will
take some minutes):

``` R
macrtools::macos_rtools_install()
```

Get OpenMP support:

``` R
macrtools::openmp_install()
```

If you have difficulties during the installation, check the following
resources:

- <https://mac.thecoatlessprofessor.com/macrtools/index.html>
- <https://mac.thecoatlessprofessor.com/macrtools/reference/openmp.html>

Finally,

``` R
devtools::install_github("marcosjnez/latent", force = TRUE)
```

## Funding

The package development is supported by two projects:

- The “DYNANSE: Righting the Wrongs. A Life Course Dynamics Approach for
  Non-Standard Employment” project, which has received funding from the
  European Research Council (ERC) under the European Union’s Horizon
  2020 research and innovation program (grant agreement No 864471).
  September 2024 – February 2026

- NWO-Open Science NL, Open Science Infrastructure. Garnier-Villarreal,
  M.(PI) & Jiménez Henríquez, MJ (Grant Agreement No 500.040.2469). Fast
  estimation of classical and new latent variable models in R. March
  2026 – February 2028. Budget: €249.873
