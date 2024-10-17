# latent: Confirmatory and Exploratory Factor Analysis Models

Fit measurement models.

# Authors

* Marcos Jiménez (Creator; Author)

* Mauricio Garnier-Villarreal (Author)

* Vithor R. Franco (Author)

* Francisco J. Abad (Author)

* Eduardo Garcia-Garzon (Author)

* Luis E. Garrido (Author)

# Installation in Linux and Windows

Using the `remotes` package:

    install.packages("remotes")
    remotes::install_github("marcosjnez/latent")
    
Using the `devtools` package:
    
    install.packages("devtools")
    devtools::install_github("marcosjnez/latent")

# Installation in macOS

Delete the clang4,6,7 binary and the prior version of gfortran installed:

    sudo rm -rf /usr/local/clang{4,6,7}
    sudo rm -rf /usr/local/gfortran
    sudo rm -rf /usr/local/bin/gfortran

Remove the gfortran install receipts:

    sudo rm /private/var/db/receipts/com.gnu.gfortran.bom
    sudo rm /private/var/db/receipts/com.gnu.gfortran.plist

Remove the clang4 installer receipt:

    sudo rm /private/var/db/receipts/com.rbinaries.clang4.bom
    sudo rm /private/var/db/receipts/com.rbinaries.clang4.plist

Remove the Makevars file and the Renviron file:
   
    rm ~/.R/Makevars
    rm ~/.Renviron

Install Homebrew from the terminal:

    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"                                
    echo 'eval "$(/opt/homebrew/bin/brew shellenv)"' >> ~/.zprofile
    eval "$(/opt/homebrew/bin/brew shellenv)"

Verify the installation with:

    brew doctor

Install `xcode`, if not installed:

    xcode-select --install

Install `libomp`, `llvm` and `gettext` from the terminal:

    brew install libomp
    brew install gettext
    brew install llvm

Download and install `gfortran` from https://github.com/fxcoudert/gfortran-for-macOS/releases for your macOS version.

Create the directory ~/.R/, if it does not exist:

    sudo mkdir ~/.R

Create a file named Makevars in this directory:

    sudo touch ~/.R/Makevars

Open the Makevars file and paste the following lines in it:

    LOC = /usr/local/gfortran/
    CC=$(LOC)/bin/gcc -fopenmp
    CXX=$(LOC)/bin/g++ -fopenmp
    # -O3 should be faster than -O2 (default) level optimisation ..
    CFLAGS=-g -O3 -Wall -pedantic -std=gnu99 -mtune=native -pipe
    CXXFLAGS=-g -O3 -Wall -pedantic -std=c++11 -mtune=native -pipe
    LDFLAGS=-L/usr/local/opt/gettext/lib -L$(LOC)/lib -Wl,-rpath,$(LOC)/lib
    CPPFLAGS=-I/usr/local/opt/gettext/include -I$(LOC)/include
    -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include

The Makevars file can be opened and edited from the R console:

    file.edit("~/.R/Makevars")

Install `Rcpp` and `RcppArmadillo` from the R console:

    install.packages(c('Rcpp', 'RcppArmadillo'))

Install `latent`:

    devtools::install_github("marcosjnez/latent")
