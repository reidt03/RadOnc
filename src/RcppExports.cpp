// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// INFORMATION
#define EMD_RUBNER_MAIN 1 #define R_NO_REMAP 1  #ifndef  DEBUG_LEVEL #define DEBUG_LEVEL 0 #endif /* DEBUG_LEVEL: 0 = NO MESSAGES 1 = PRINT THE NUMBER OF ITERATIONS AND THE FINAL RESULT 2 = PRINT THE RESULT AFTER EVERY ITERATION 3 = PRINT ALSO THE FLOW AFTER EVERY ITERATION 4 = PRINT A LOT OF INFORMATION(PROBABLY USEFUL ONLY FOR THE AUTHOR) */   #define MAX_SIG_SIZE1 (MAX_SIG_SIZE+1);
RcppExport SEXP _RadOnc_INFORMATION(SEXP (MAX_SIG_SIZE+1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< PROBABLY USEFUL ONLY FOR THE AUTHOR) */   #define MAX_SIG_SIZE1 >::type (MAX_SIG_SIZE+1((MAX_SIG_SIZE+1SEXP);
    rcpp_result_gen = Rcpp::wrap(INFORMATION((MAX_SIG_SIZE+1));
    return rcpp_result_gen;
END_RCPP
}
// timesTwo
NumericVector timesTwo(NumericVector x);
RcppExport SEXP _RadOnc_timesTwo(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(timesTwo(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_RadOnc_INFORMATION", (DL_FUNC) &_RadOnc_INFORMATION, 1},
    {"_RadOnc_timesTwo", (DL_FUNC) &_RadOnc_timesTwo, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_RadOnc(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
