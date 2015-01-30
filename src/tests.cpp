/* Export code to run various tests of CppAD code
*/

#include <tests/test_functors.h>


//' @name CppADutils-tests
//' @title tests
//' @description tests
//' @param X numeric vector (not matrix at this time)
//[[Rcpp::export]]
void CppADutils_tests(const NumericVector& X) {}

//' @inheritParams CppADutils-tests
//' @rdname CppADutils-tests
//[[Rcpp::export]]
Rcpp::List cppad_lbeta1(const NumericVector& X) {
  run_test<lbeta1> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}

//' @inheritParams CppADutils-tests
//' @rdname CppADutils-tests
//[[Rcpp::export]]
Rcpp::List cppad_dnorm_log(const NumericVector& X) {
  run_test<dnorm_log_test> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}

//' @inheritParams CppADutils-tests
//' @rdname CppADutils-tests
//[[Rcpp::export]]
Rcpp::List cppad_dbeta_log(const NumericVector& X) {
  run_test<dbeta_log_test> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}

//' @inheritParams CppADutils-tests
//' @rdname CppADutils-tests
//[[Rcpp::export]]
Rcpp::List cppad_dlogitbeta_log(const NumericVector& X) {
  run_test<dlogitbeta_log_test> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}


//' @inheritParams CppADutils-tests
//' @rdname CppADutils-tests
//[[Rcpp::export]]
Rcpp::List cppad_lgamma1(const NumericVector& X) {
  run_test<lgamma1> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}

//' @inheritParams CppADutils-tests
//' @rdname CppADutils-tests
//[[Rcpp::export]]
Rcpp::List cppad_lgamma1p(const NumericVector& X) {
  run_test<lgamma1p_test> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}

//' @inheritParams CppADutils-tests
//' @rdname CppADutils-tests
//[[Rcpp::export]]
Rcpp::List cppad_lgammaexp(const NumericVector& X) {
  run_test<lgammaexp_test> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}

//' @inheritParams CppADutils-tests
//' @rdname CppADutils-tests
//[[Rcpp::export]]
Rcpp::List cppad_log1p(const NumericVector& X) {
  run_test<log1p_test> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}

//' @inheritParams CppADutils-tests
//' @rdname CppADutils-tests
//[[Rcpp::export]]
Rcpp::List cppad_expm1(const NumericVector& X) {
  run_test<expm1_test> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}

//' @inheritParams CppADutils-tests
//' @rdname CppADutils-tests
//[[Rcpp::export]]
Rcpp::List cppad_log1pexp(const NumericVector& X) {
  run_test<log1pexp_test> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}

//' @inheritParams CppADutils-tests
//' @rdname CppADutils-tests
//[[Rcpp::export]]
Rcpp::List cppad_log1pmx(const NumericVector& X) {
  run_test<log1pmx_test> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}

//' @inheritParams CppADutils-tests
//' @rdname CppADutils-tests
//[[Rcpp::export]]
Rcpp::List cppad_loginvlogit(const NumericVector& X) {
  run_test<loginvlogit_test> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}

//' @inheritParams CppADutils-tests
//' @rdname CppADutils-tests
//[[Rcpp::export]]
Rcpp::List cppad_invlogit(const NumericVector& X) {
  run_test<invlogit_test> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}

//' @inheritParams CppADutils-tests
//' @rdname CppADutils-tests
//[[Rcpp::export]]
Rcpp::List cppad_incgamma(const NumericVector& X) {
  run_test<incgamma_test> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}

//' @inheritParams CppADutils-tests
//' @rdname CppADutils-tests
//[[Rcpp::export]]
Rcpp::List cppad_incbeta(const NumericVector& X) {
  run_test<incbeta_test> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}



//' @inheritParams CppADutils-tests
//' @rdname CppADutils-tests
//[[Rcpp::export]]
Rcpp::List cppad_incbeta2(const NumericVector& X) {
  run_test<incbeta_test2> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}

//' @inheritParams CppADutils-tests
//' @rdname CppADutils-tests
//[[Rcpp::export]]
Rcpp::List cppad_lgammaLogExp(const NumericVector& X) {
  run_test<lgammaLogExp> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}

//' @inheritParams CppADutils-tests
//' @rdname CppADutils-tests
//[[Rcpp::export]]
Rcpp::List cppad_dt_log(const NumericVector& X) {
  run_test<dt_log_test> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}


//' @inheritParams CppADutils-tests
//' @rdname CppADutils-tests
//[[Rcpp::export]]
Rcpp::List cppad_dhalft_log(const NumericVector& X) {
  run_test<dhalft_log_test> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}


//' @inheritParams CppADutils-tests
//' @rdname CppADutils-tests
//[[Rcpp::export]]
Rcpp::List cppad_dhalft_log2(const NumericVector& X) {
  run_test<dhalft_log_test2> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}

//' @inheritParams CppADutils-tests
//' @rdname CppADutils-tests
//[[Rcpp::export]]
Rcpp::List cppad_dnormTrunc0_log(const NumericVector& X) {
  run_test<dnormTrunc0_log_test> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}

//' @inheritParams CppADutils-tests
//' @rdname CppADutils-tests
//[[Rcpp::export]]
Rcpp::List cppad_pnorm_log(const NumericVector& X) {
  run_test<pnorm_log_test> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}

