// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// genoEM
Rcpp::NumericVector genoEM(Rcpp::NumericVector refCounts, Rcpp::NumericVector altCounts, int ploidy, Rcpp::NumericVector h, Rcpp::NumericVector eps, int mrep, double mdiff, bool returnAll);
RcppExport SEXP _tripsAndDipR_genoEM(SEXP refCountsSEXP, SEXP altCountsSEXP, SEXP ploidySEXP, SEXP hSEXP, SEXP epsSEXP, SEXP mrepSEXP, SEXP mdiffSEXP, SEXP returnAllSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type refCounts(refCountsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type altCounts(altCountsSEXP);
    Rcpp::traits::input_parameter< int >::type ploidy(ploidySEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type h(hSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type mrep(mrepSEXP);
    Rcpp::traits::input_parameter< double >::type mdiff(mdiffSEXP);
    Rcpp::traits::input_parameter< bool >::type returnAll(returnAllSEXP);
    rcpp_result_gen = Rcpp::wrap(genoEM(refCounts, altCounts, ploidy, h, eps, mrep, mdiff, returnAll));
    return rcpp_result_gen;
END_RCPP
}
// genoEM_noise
Rcpp::NumericVector genoEM_noise(Rcpp::NumericVector refCounts, Rcpp::NumericVector altCounts, int ploidy, Rcpp::NumericVector h, Rcpp::NumericVector eps, int mrep, double mdiff, bool returnAll);
RcppExport SEXP _tripsAndDipR_genoEM_noise(SEXP refCountsSEXP, SEXP altCountsSEXP, SEXP ploidySEXP, SEXP hSEXP, SEXP epsSEXP, SEXP mrepSEXP, SEXP mdiffSEXP, SEXP returnAllSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type refCounts(refCountsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type altCounts(altCountsSEXP);
    Rcpp::traits::input_parameter< int >::type ploidy(ploidySEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type h(hSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type mrep(mrepSEXP);
    Rcpp::traits::input_parameter< double >::type mdiff(mdiffSEXP);
    Rcpp::traits::input_parameter< bool >::type returnAll(returnAllSEXP);
    rcpp_result_gen = Rcpp::wrap(genoEM_noise(refCounts, altCounts, ploidy, h, eps, mrep, mdiff, returnAll));
    return rcpp_result_gen;
END_RCPP
}
// llh_calc_BB_noise
Rcpp::NumericVector llh_calc_BB_noise(Rcpp::NumericVector refCounts, Rcpp::NumericVector altCounts, Rcpp::NumericVector tau, Rcpp::NumericVector mixWeights, int ploidy, Rcpp::NumericVector h, Rcpp::NumericVector eps);
RcppExport SEXP _tripsAndDipR_llh_calc_BB_noise(SEXP refCountsSEXP, SEXP altCountsSEXP, SEXP tauSEXP, SEXP mixWeightsSEXP, SEXP ploidySEXP, SEXP hSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type refCounts(refCountsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type altCounts(altCountsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mixWeights(mixWeightsSEXP);
    Rcpp::traits::input_parameter< int >::type ploidy(ploidySEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type h(hSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(llh_calc_BB_noise(refCounts, altCounts, tau, mixWeights, ploidy, h, eps));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_tripsAndDipR_genoEM", (DL_FUNC) &_tripsAndDipR_genoEM, 8},
    {"_tripsAndDipR_genoEM_noise", (DL_FUNC) &_tripsAndDipR_genoEM_noise, 8},
    {"_tripsAndDipR_llh_calc_BB_noise", (DL_FUNC) &_tripsAndDipR_llh_calc_BB_noise, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_tripsAndDipR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
