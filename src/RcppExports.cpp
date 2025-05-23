// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// construct_ADchol_Rcpp
List construct_ADchol_Rcpp(Rcpp::S4 obj_spam, const List& P_list);
RcppExport SEXP _LMMsolver_construct_ADchol_Rcpp(SEXP obj_spamSEXP, SEXP P_listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type obj_spam(obj_spamSEXP);
    Rcpp::traits::input_parameter< const List& >::type P_list(P_listSEXP);
    rcpp_result_gen = Rcpp::wrap(construct_ADchol_Rcpp(obj_spam, P_list));
    return rcpp_result_gen;
END_RCPP
}
// dlogdet
NumericVector dlogdet(Rcpp::S4 obj, NumericVector theta, Nullable<NumericVector> b_);
RcppExport SEXP _LMMsolver_dlogdet(SEXP objSEXP, SEXP thetaSEXP, SEXP b_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type obj(objSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< Nullable<NumericVector> >::type b_(b_SEXP);
    rcpp_result_gen = Rcpp::wrap(dlogdet(obj, theta, b_));
    return rcpp_result_gen;
END_RCPP
}
// partialDerivCholesky
NumericVector partialDerivCholesky(Rcpp::S4 obj);
RcppExport SEXP _LMMsolver_partialDerivCholesky(SEXP objSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type obj(objSEXP);
    rcpp_result_gen = Rcpp::wrap(partialDerivCholesky(obj));
    return rcpp_result_gen;
END_RCPP
}
// diagXCinvXt
NumericVector diagXCinvXt(Rcpp::S4 obj, Rcpp::S4 transposeX);
RcppExport SEXP _LMMsolver_diagXCinvXt(SEXP objSEXP, SEXP transposeXSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type obj(objSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type transposeX(transposeXSEXP);
    rcpp_result_gen = Rcpp::wrap(diagXCinvXt(obj, transposeX));
    return rcpp_result_gen;
END_RCPP
}
// GetIntVector
IntegerVector GetIntVector(Rcpp::S4 obj, const String& slotName, int ArrayIndexing);
RcppExport SEXP _LMMsolver_GetIntVector(SEXP objSEXP, SEXP slotNameSEXP, SEXP ArrayIndexingSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type obj(objSEXP);
    Rcpp::traits::input_parameter< const String& >::type slotName(slotNameSEXP);
    Rcpp::traits::input_parameter< int >::type ArrayIndexing(ArrayIndexingSEXP);
    rcpp_result_gen = Rcpp::wrap(GetIntVector(obj, slotName, ArrayIndexing));
    return rcpp_result_gen;
END_RCPP
}
// RowKron
Rcpp::S4 RowKron(Rcpp::S4 sX1, Rcpp::S4 sX2);
RcppExport SEXP _LMMsolver_RowKron(SEXP sX1SEXP, SEXP sX2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type sX1(sX1SEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type sX2(sX2SEXP);
    rcpp_result_gen = Rcpp::wrap(RowKron(sX1, sX2));
    return rcpp_result_gen;
END_RCPP
}
// MatrixProduct
Rcpp::S4 MatrixProduct(Rcpp::S4 sA, Rcpp::S4 sB);
RcppExport SEXP _LMMsolver_MatrixProduct(SEXP sASEXP, SEXP sBSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type sA(sASEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type sB(sBSEXP);
    rcpp_result_gen = Rcpp::wrap(MatrixProduct(sA, sB));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_LMMsolver_construct_ADchol_Rcpp", (DL_FUNC) &_LMMsolver_construct_ADchol_Rcpp, 2},
    {"_LMMsolver_dlogdet", (DL_FUNC) &_LMMsolver_dlogdet, 3},
    {"_LMMsolver_partialDerivCholesky", (DL_FUNC) &_LMMsolver_partialDerivCholesky, 1},
    {"_LMMsolver_diagXCinvXt", (DL_FUNC) &_LMMsolver_diagXCinvXt, 2},
    {"_LMMsolver_GetIntVector", (DL_FUNC) &_LMMsolver_GetIntVector, 3},
    {"_LMMsolver_RowKron", (DL_FUNC) &_LMMsolver_RowKron, 2},
    {"_LMMsolver_MatrixProduct", (DL_FUNC) &_LMMsolver_MatrixProduct, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_LMMsolver(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
