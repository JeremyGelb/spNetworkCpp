#' @useDynLib spNetworkCpp
#' @importFrom Rcpp sourceCpp
NULL


.onUnload <- function (libpath) {
  library.dynam.unload("spnetwork", libpath)
}
