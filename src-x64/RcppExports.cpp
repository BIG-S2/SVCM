// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// aspV2
List aspV2(NumericMatrix BetaIn, NumericMatrix EtaIn, NumericVector HsIn, NumericVector VnsIn, NumericMatrix VnsidIn, NumericMatrix VdsIn, NumericMatrix XIn, NumericVector ChiseqIn, NumericVector SeIn, double CnIn);
RcppExport SEXP SVCM_aspV2(SEXP BetaInSEXP, SEXP EtaInSEXP, SEXP HsInSEXP, SEXP VnsInSEXP, SEXP VnsidInSEXP, SEXP VdsInSEXP, SEXP XInSEXP, SEXP ChiseqInSEXP, SEXP SeInSEXP, SEXP CnInSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type BetaIn(BetaInSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type EtaIn(EtaInSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type HsIn(HsInSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type VnsIn(VnsInSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type VnsidIn(VnsidInSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type VdsIn(VdsInSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type XIn(XInSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ChiseqIn(ChiseqInSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type SeIn(SeInSEXP);
    Rcpp::traits::input_parameter< double >::type CnIn(CnInSEXP);
    rcpp_result_gen = Rcpp::wrap(aspV2(BetaIn, EtaIn, HsIn, VnsIn, VnsidIn, VdsIn, XIn, ChiseqIn, SeIn, CnIn));
    return rcpp_result_gen;
END_RCPP
}
// aspV2NYU
List aspV2NYU(NumericMatrix BetaIn, NumericMatrix EtaIn, NumericVector HsIn, NumericVector VnsIn, NumericMatrix VnsidIn, NumericMatrix VdsIn, NumericMatrix XIn, NumericVector ChiseqIn, NumericVector SeIn, double CnIn, NumericMatrix IXX);
RcppExport SEXP SVCM_aspV2NYU(SEXP BetaInSEXP, SEXP EtaInSEXP, SEXP HsInSEXP, SEXP VnsInSEXP, SEXP VnsidInSEXP, SEXP VdsInSEXP, SEXP XInSEXP, SEXP ChiseqInSEXP, SEXP SeInSEXP, SEXP CnInSEXP, SEXP IXXSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type BetaIn(BetaInSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type EtaIn(EtaInSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type HsIn(HsInSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type VnsIn(VnsInSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type VnsidIn(VnsidInSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type VdsIn(VdsInSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type XIn(XInSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ChiseqIn(ChiseqInSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type SeIn(SeInSEXP);
    Rcpp::traits::input_parameter< double >::type CnIn(CnInSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type IXX(IXXSEXP);
    rcpp_result_gen = Rcpp::wrap(aspV2NYU(BetaIn, EtaIn, HsIn, VnsIn, VnsidIn, VdsIn, XIn, ChiseqIn, SeIn, CnIn, IXX));
    return rcpp_result_gen;
END_RCPP
}
// mlrWresd
List mlrWresd(NumericMatrix X, NumericMatrix Y);
RcppExport SEXP SVCM_mlrWresd(SEXP XSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(mlrWresd(X, Y));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP SVCM_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// siv
List siv(NumericMatrix RImg, NumericMatrix XYZ, NumericVector H);
RcppExport SEXP SVCM_siv(SEXP RImgSEXP, SEXP XYZSEXP, SEXP HSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type RImg(RImgSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type XYZ(XYZSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type H(HSEXP);
    rcpp_result_gen = Rcpp::wrap(siv(RImg, XYZ, H));
    return rcpp_result_gen;
END_RCPP
}
