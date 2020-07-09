// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// continuous_nkde_cpp_arma_sparse
DataFrame continuous_nkde_cpp_arma_sparse(List neighbour_list, NumericVector events, NumericVector weights, DataFrame samples, NumericVector bws, std::string kernel_name, DataFrame nodes, DataFrame line_list, int max_depth, bool verbose);
RcppExport SEXP _spNetworkCpp_continuous_nkde_cpp_arma_sparse(SEXP neighbour_listSEXP, SEXP eventsSEXP, SEXP weightsSEXP, SEXP samplesSEXP, SEXP bwsSEXP, SEXP kernel_nameSEXP, SEXP nodesSEXP, SEXP line_listSEXP, SEXP max_depthSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type neighbour_list(neighbour_listSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type events(eventsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type samples(samplesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bws(bwsSEXP);
    Rcpp::traits::input_parameter< std::string >::type kernel_name(kernel_nameSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type nodes(nodesSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type line_list(line_listSEXP);
    Rcpp::traits::input_parameter< int >::type max_depth(max_depthSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(continuous_nkde_cpp_arma_sparse(neighbour_list, events, weights, samples, bws, kernel_name, nodes, line_list, max_depth, verbose));
    return rcpp_result_gen;
END_RCPP
}
// continuous_nkde_cpp_arma
DataFrame continuous_nkde_cpp_arma(List neighbour_list, NumericVector events, NumericVector weights, DataFrame samples, NumericVector bws, std::string kernel_name, DataFrame nodes, DataFrame line_list, int max_depth, bool verbose);
RcppExport SEXP _spNetworkCpp_continuous_nkde_cpp_arma(SEXP neighbour_listSEXP, SEXP eventsSEXP, SEXP weightsSEXP, SEXP samplesSEXP, SEXP bwsSEXP, SEXP kernel_nameSEXP, SEXP nodesSEXP, SEXP line_listSEXP, SEXP max_depthSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type neighbour_list(neighbour_listSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type events(eventsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type samples(samplesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bws(bwsSEXP);
    Rcpp::traits::input_parameter< std::string >::type kernel_name(kernel_nameSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type nodes(nodesSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type line_list(line_listSEXP);
    Rcpp::traits::input_parameter< int >::type max_depth(max_depthSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(continuous_nkde_cpp_arma(neighbour_list, events, weights, samples, bws, kernel_name, nodes, line_list, max_depth, verbose));
    return rcpp_result_gen;
END_RCPP
}
// discontinuous_nkde_cpp_arma_sparse
DataFrame discontinuous_nkde_cpp_arma_sparse(List neighbour_list, NumericVector events, NumericVector weights, DataFrame samples, NumericVector bws, std::string kernel_name, DataFrame nodes, DataFrame line_list, int max_depth, bool verbose);
RcppExport SEXP _spNetworkCpp_discontinuous_nkde_cpp_arma_sparse(SEXP neighbour_listSEXP, SEXP eventsSEXP, SEXP weightsSEXP, SEXP samplesSEXP, SEXP bwsSEXP, SEXP kernel_nameSEXP, SEXP nodesSEXP, SEXP line_listSEXP, SEXP max_depthSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type neighbour_list(neighbour_listSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type events(eventsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type samples(samplesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bws(bwsSEXP);
    Rcpp::traits::input_parameter< std::string >::type kernel_name(kernel_nameSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type nodes(nodesSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type line_list(line_listSEXP);
    Rcpp::traits::input_parameter< int >::type max_depth(max_depthSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(discontinuous_nkde_cpp_arma_sparse(neighbour_list, events, weights, samples, bws, kernel_name, nodes, line_list, max_depth, verbose));
    return rcpp_result_gen;
END_RCPP
}
// discontinuous_nkde_cpp_arma
DataFrame discontinuous_nkde_cpp_arma(List neighbour_list, NumericVector events, NumericVector weights, DataFrame samples, NumericVector bws, std::string kernel_name, DataFrame nodes, DataFrame line_list, int max_depth, bool verbose);
RcppExport SEXP _spNetworkCpp_discontinuous_nkde_cpp_arma(SEXP neighbour_listSEXP, SEXP eventsSEXP, SEXP weightsSEXP, SEXP samplesSEXP, SEXP bwsSEXP, SEXP kernel_nameSEXP, SEXP nodesSEXP, SEXP line_listSEXP, SEXP max_depthSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type neighbour_list(neighbour_listSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type events(eventsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type samples(samplesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bws(bwsSEXP);
    Rcpp::traits::input_parameter< std::string >::type kernel_name(kernel_nameSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type nodes(nodesSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type line_list(line_listSEXP);
    Rcpp::traits::input_parameter< int >::type max_depth(max_depthSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(discontinuous_nkde_cpp_arma(neighbour_list, events, weights, samples, bws, kernel_name, nodes, line_list, max_depth, verbose));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_spNetworkCpp_continuous_nkde_cpp_arma_sparse", (DL_FUNC) &_spNetworkCpp_continuous_nkde_cpp_arma_sparse, 10},
    {"_spNetworkCpp_continuous_nkde_cpp_arma", (DL_FUNC) &_spNetworkCpp_continuous_nkde_cpp_arma, 10},
    {"_spNetworkCpp_discontinuous_nkde_cpp_arma_sparse", (DL_FUNC) &_spNetworkCpp_discontinuous_nkde_cpp_arma_sparse, 10},
    {"_spNetworkCpp_discontinuous_nkde_cpp_arma", (DL_FUNC) &_spNetworkCpp_discontinuous_nkde_cpp_arma, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_spNetworkCpp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
