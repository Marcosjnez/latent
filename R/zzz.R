.onAttach <- function(libname, pkgname) {

  path <- system.file("extdata", "latent_mascot.png", package = pkgname)
  if (nzchar(path)) {
    if (interactive() && requireNamespace("png", quietly = TRUE) && requireNamespace("grid", quietly = TRUE)) {
      img <- png::readPNG(path)
      grid::grid.raster(img)
    }
    packageStartupMessage("👻 Welcome to latent!")
  }

  version <- as.character(utils::packageVersion(pkgname))
  msg <- sprintf("
      __      __             __
     / /_____/ /____  ____  / /_
    / / __  / __/ _ \\/ __ \\/ __/
   / / /_/ / /_/  __/ / / / /_
  /_/\\____/\\__/\\___/_/ /_/\\__/   version %s

Type 'citation(\"%s\")' for how to cite this package in publications.
", version, pkgname)

  packageStartupMessage(msg)

}
