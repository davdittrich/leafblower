# R_init_leafblower() in r_bridge.cpp is called automatically by R when the
# shared library is loaded (R calls R_init_<pkgname> on dlopen). No .onLoad()
# needed. This empty .onLoad is present to suppress R CMD check NOTE about
# an absent zzz.R on some R versions.
.onLoad <- function(libname, pkgname) invisible(NULL)
