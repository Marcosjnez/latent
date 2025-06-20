# latent: Latent Class and Factor Analysis Models

<img src="man/figures/standard.png" align="right" />

Fit measurement models with discrete or continuous latent variables.

# Installation in Windows and Linux

    devtools::install_github("marcosjnez/latent")

# Installation in macOS

Install Command Line Tools:

    xcode-select --install

Install Homebrew:
 
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
    
Install OpenMP support:

    brew update && brew install libomp
    
Update Makevars

    add the following lines to ~/.R/Makevars.

On Intel Macs:

    CPPFLAGS += -Xclang -fopenmp
    LDFLAGS += -lomp

On Apple Silicon Macs (make sure this goes at the end of the file):

    LDFLAGS += -L/opt/homebrew/opt/libomp/lib -lomp
    CPPFLAGS += -I/opt/homebrew/opt/libomp/include -Xclang -fopenmp

Finally,

    devtools::install_github("marcosjnez/latent")
