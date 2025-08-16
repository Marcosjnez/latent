# latent: Latent Class and Factor Analysis Models

<div style = "padding-top:1em; padding-bottom: 0.5em;">
<img src="man/figures/standard.png" width = 135 align="right" />
</div>

Fit measurement models with discrete or continuous latent variables.

# Installation in Windows and Linux

    devtools::install_github("marcosjnez/latent")

# Installation in macOS

Install the `macrtools` R package of James Balamuta:

    # install.packages("remotes")
    remotes::install_github("coatless-mac/macrtools")

Install Command Line Tools and the R Compilation Toolchain (this will take some minutes):

    macrtools::macos_rtools_install()

Get OpenMP support:

    macrtools::openmp_install()

If you have difficulties during the installation, check the following resources:

* https://mac.thecoatlessprofessor.com/macrtools/index.html
* https://mac.thecoatlessprofessor.com/macrtools/reference/openmp.html
    
Finally,

    devtools::install_github("marcosjnez/latent")

## Funding
The package development is supported by the “DYNANSE: Righting the Wrongs. A Life Course Dynamics Approach for Non-Standard Employment” project, which has received funding from the European Research Council (ERC) under the European Union’s Horizon 2020 research and innovation programme (grant agreement No 864471).
