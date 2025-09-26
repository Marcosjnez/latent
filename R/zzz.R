# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 22/05/2025

.onAttach <- function(libname, pkgname) {

  path <- system.file("extdata", "standard.png", package = pkgname)
  if (nzchar(path)) {
    if (interactive() && requireNamespace("png", quietly = TRUE) && requireNamespace("grid", quietly = TRUE)) {
      img <- png::readPNG(path)
      grid::grid.raster(img)
    }
  }

  packageStartupMessage("👻 Welcome to latent!")

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
