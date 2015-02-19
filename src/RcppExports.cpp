// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// MVN_test
NumericVector MVN_test(NumericMatrix X_, NumericMatrix mu_, NumericMatrix G_, bool isPrec);
RcppExport SEXP CppADutils_MVN_test(SEXP X_SEXP, SEXP mu_SEXP, SEXP G_SEXP, SEXP isPrecSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type X_(X_SEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type mu_(mu_SEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type G_(G_SEXP );
        Rcpp::traits::input_parameter< bool >::type isPrec(isPrecSEXP );
        NumericVector __result = MVN_test(X_, mu_, G_, isPrec);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// Sparse_MVN_test
NumericVector Sparse_MVN_test(NumericMatrix X_, NumericMatrix mu_, Rcpp::S4 S_, bool isPrec);
RcppExport SEXP CppADutils_Sparse_MVN_test(SEXP X_SEXP, SEXP mu_SEXP, SEXP S_SEXP, SEXP isPrecSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type X_(X_SEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type mu_(mu_SEXP );
        Rcpp::traits::input_parameter< Rcpp::S4 >::type S_(S_SEXP );
        Rcpp::traits::input_parameter< bool >::type isPrec(isPrecSEXP );
        NumericVector __result = Sparse_MVN_test(X_, mu_, S_, isPrec);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// invlogit
Rcpp::NumericVector invlogit(const NumericVector& X);
RcppExport SEXP CppADutils_invlogit(SEXP XSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const NumericVector& >::type X(XSEXP );
        Rcpp::NumericVector __result = invlogit(X);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// log1pmx
Rcpp::NumericVector log1pmx(const NumericVector& X);
RcppExport SEXP CppADutils_log1pmx(SEXP XSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const NumericVector& >::type X(XSEXP );
        Rcpp::NumericVector __result = log1pmx(X);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// loginvlogit
Rcpp::NumericVector loginvlogit(const NumericVector& X);
RcppExport SEXP CppADutils_loginvlogit(SEXP XSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const NumericVector& >::type X(XSEXP );
        Rcpp::NumericVector __result = loginvlogit(X);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// log1pexp
Rcpp::NumericVector log1pexp(const NumericVector& X);
RcppExport SEXP CppADutils_log1pexp(SEXP XSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const NumericVector& >::type X(XSEXP );
        Rcpp::NumericVector __result = log1pexp(X);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// lgamma1p
Rcpp::NumericVector lgamma1p(const NumericVector& X);
RcppExport SEXP CppADutils_lgamma1p(SEXP XSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const NumericVector& >::type X(XSEXP );
        Rcpp::NumericVector __result = lgamma1p(X);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// lgammaexp
Rcpp::NumericVector lgammaexp(const NumericVector& X);
RcppExport SEXP CppADutils_lgammaexp(SEXP XSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const NumericVector& >::type X(XSEXP );
        Rcpp::NumericVector __result = lgammaexp(X);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// CppADutils_tests
void CppADutils_tests(const NumericVector& X);
RcppExport SEXP CppADutils_CppADutils_tests(SEXP XSEXP) {
BEGIN_RCPP
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const NumericVector& >::type X(XSEXP );
        CppADutils_tests(X);
    }
    return R_NilValue;
END_RCPP
}
// cppad_lbeta1
Rcpp::List cppad_lbeta1(const NumericVector& X);
RcppExport SEXP CppADutils_cppad_lbeta1(SEXP XSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const NumericVector& >::type X(XSEXP );
        Rcpp::List __result = cppad_lbeta1(X);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// cppad_dnorm_log
Rcpp::List cppad_dnorm_log(const NumericVector& X);
RcppExport SEXP CppADutils_cppad_dnorm_log(SEXP XSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const NumericVector& >::type X(XSEXP );
        Rcpp::List __result = cppad_dnorm_log(X);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// cppad_dbeta_log
Rcpp::List cppad_dbeta_log(const NumericVector& X);
RcppExport SEXP CppADutils_cppad_dbeta_log(SEXP XSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const NumericVector& >::type X(XSEXP );
        Rcpp::List __result = cppad_dbeta_log(X);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// cppad_dlogitbeta_log
Rcpp::List cppad_dlogitbeta_log(const NumericVector& X);
RcppExport SEXP CppADutils_cppad_dlogitbeta_log(SEXP XSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const NumericVector& >::type X(XSEXP );
        Rcpp::List __result = cppad_dlogitbeta_log(X);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// cppad_lgamma1
Rcpp::List cppad_lgamma1(const NumericVector& X);
RcppExport SEXP CppADutils_cppad_lgamma1(SEXP XSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const NumericVector& >::type X(XSEXP );
        Rcpp::List __result = cppad_lgamma1(X);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// cppad_lgamma1p
Rcpp::List cppad_lgamma1p(const NumericVector& X);
RcppExport SEXP CppADutils_cppad_lgamma1p(SEXP XSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const NumericVector& >::type X(XSEXP );
        Rcpp::List __result = cppad_lgamma1p(X);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// cppad_lgammaexp
Rcpp::List cppad_lgammaexp(const NumericVector& X);
RcppExport SEXP CppADutils_cppad_lgammaexp(SEXP XSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const NumericVector& >::type X(XSEXP );
        Rcpp::List __result = cppad_lgammaexp(X);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// cppad_log1p
Rcpp::List cppad_log1p(const NumericVector& X);
RcppExport SEXP CppADutils_cppad_log1p(SEXP XSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const NumericVector& >::type X(XSEXP );
        Rcpp::List __result = cppad_log1p(X);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// cppad_expm1
Rcpp::List cppad_expm1(const NumericVector& X);
RcppExport SEXP CppADutils_cppad_expm1(SEXP XSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const NumericVector& >::type X(XSEXP );
        Rcpp::List __result = cppad_expm1(X);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// cppad_log1pexp
Rcpp::List cppad_log1pexp(const NumericVector& X);
RcppExport SEXP CppADutils_cppad_log1pexp(SEXP XSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const NumericVector& >::type X(XSEXP );
        Rcpp::List __result = cppad_log1pexp(X);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// cppad_log1pmx
Rcpp::List cppad_log1pmx(const NumericVector& X);
RcppExport SEXP CppADutils_cppad_log1pmx(SEXP XSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const NumericVector& >::type X(XSEXP );
        Rcpp::List __result = cppad_log1pmx(X);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// cppad_loginvlogit
Rcpp::List cppad_loginvlogit(const NumericVector& X);
RcppExport SEXP CppADutils_cppad_loginvlogit(SEXP XSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const NumericVector& >::type X(XSEXP );
        Rcpp::List __result = cppad_loginvlogit(X);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// cppad_invlogit
Rcpp::List cppad_invlogit(const NumericVector& X);
RcppExport SEXP CppADutils_cppad_invlogit(SEXP XSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const NumericVector& >::type X(XSEXP );
        Rcpp::List __result = cppad_invlogit(X);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// cppad_incgamma
Rcpp::List cppad_incgamma(const NumericVector& X);
RcppExport SEXP CppADutils_cppad_incgamma(SEXP XSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const NumericVector& >::type X(XSEXP );
        Rcpp::List __result = cppad_incgamma(X);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// cppad_incbeta
Rcpp::List cppad_incbeta(const NumericVector& X);
RcppExport SEXP CppADutils_cppad_incbeta(SEXP XSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const NumericVector& >::type X(XSEXP );
        Rcpp::List __result = cppad_incbeta(X);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// cppad_incbeta2
Rcpp::List cppad_incbeta2(const NumericVector& X);
RcppExport SEXP CppADutils_cppad_incbeta2(SEXP XSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const NumericVector& >::type X(XSEXP );
        Rcpp::List __result = cppad_incbeta2(X);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// cppad_lgammaLogExp
Rcpp::List cppad_lgammaLogExp(const NumericVector& X);
RcppExport SEXP CppADutils_cppad_lgammaLogExp(SEXP XSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const NumericVector& >::type X(XSEXP );
        Rcpp::List __result = cppad_lgammaLogExp(X);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// cppad_dt_log
Rcpp::List cppad_dt_log(const NumericVector& X);
RcppExport SEXP CppADutils_cppad_dt_log(SEXP XSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const NumericVector& >::type X(XSEXP );
        Rcpp::List __result = cppad_dt_log(X);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// cppad_dhalft_log
Rcpp::List cppad_dhalft_log(const NumericVector& X);
RcppExport SEXP CppADutils_cppad_dhalft_log(SEXP XSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const NumericVector& >::type X(XSEXP );
        Rcpp::List __result = cppad_dhalft_log(X);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// cppad_dhalft_log2
Rcpp::List cppad_dhalft_log2(const NumericVector& X);
RcppExport SEXP CppADutils_cppad_dhalft_log2(SEXP XSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const NumericVector& >::type X(XSEXP );
        Rcpp::List __result = cppad_dhalft_log2(X);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// cppad_dnormTrunc0_log
Rcpp::List cppad_dnormTrunc0_log(const NumericVector& X);
RcppExport SEXP CppADutils_cppad_dnormTrunc0_log(SEXP XSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const NumericVector& >::type X(XSEXP );
        Rcpp::List __result = cppad_dnormTrunc0_log(X);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// cppad_atan2a
Rcpp::List cppad_atan2a(const NumericVector& X);
RcppExport SEXP CppADutils_cppad_atan2a(SEXP XSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const NumericVector& >::type X(XSEXP );
        Rcpp::List __result = cppad_atan2a(X);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// cppad_atan2b
Rcpp::List cppad_atan2b(const NumericVector& X);
RcppExport SEXP CppADutils_cppad_atan2b(SEXP XSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const NumericVector& >::type X(XSEXP );
        Rcpp::List __result = cppad_atan2b(X);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// cppad_pnorm_log
Rcpp::List cppad_pnorm_log(const NumericVector& X);
RcppExport SEXP CppADutils_cppad_pnorm_log(SEXP XSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const NumericVector& >::type X(XSEXP );
        Rcpp::List __result = cppad_pnorm_log(X);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// Wish_test
double Wish_test(NumericMatrix X_, double nu_, NumericMatrix S_);
RcppExport SEXP CppADutils_Wish_test(SEXP X_SEXP, SEXP nu_SEXP, SEXP S_SEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type X_(X_SEXP );
        Rcpp::traits::input_parameter< double >::type nu_(nu_SEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type S_(S_SEXP );
        double __result = Wish_test(X_, nu_, S_);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// Inv_Wish_test
double Inv_Wish_test(NumericMatrix X_, double nu_, NumericMatrix S_);
RcppExport SEXP CppADutils_Inv_Wish_test(SEXP X_SEXP, SEXP nu_SEXP, SEXP S_SEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type X_(X_SEXP );
        Rcpp::traits::input_parameter< double >::type nu_(nu_SEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type S_(S_SEXP );
        double __result = Inv_Wish_test(X_, nu_, S_);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
