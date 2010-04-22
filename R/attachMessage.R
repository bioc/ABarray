.onAttach <- function(libname, pkgname) {
    packageStartupMessage(
        paste("\n\nWelcome to Applied Biosystems AB1700\n",
              "    This package performs analysis for AB1700\n",
              "    gene expression data\n", sep = ""))
}
