/* Export code to run various tests of CppAD code
*/

#include <tests/test_functors.h>
#include <tests/test_compound.h>

//[[Rcpp::export]]
Rcpp::List cppad_log1p_sum(const NumericVector& X) {
  run_test<log1p_sum> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}

/*
//[[Rcpp::export]]
Rcpp::List cppad_log1p_prod(const NumericVector& X) {
  run_test<log1p_prod> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}

//[[Rcpp::export]]
Rcpp::List cppad_lgamma1p_sum(const NumericVector& X) {
  run_test<lgamma1p_sum> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}

//[[Rcpp::export]]
Rcpp::List cppad_lgamma1p_prod(const NumericVector& X) {
  run_test<lgamma1p_prod> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}

//[[Rcpp::export]]
Rcpp::List cppad_log1pmx_sum(const NumericVector& X) {
  run_test<log1pmx_sum> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}

//[[Rcpp::export]]
Rcpp::List cppad_log1pmx_prod(const NumericVector& X) {
  run_test<log1pmx_prod> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}

//[[Rcpp::export]]
Rcpp::List cppad_expm1_sum(const NumericVector& X) {
  run_test<expm1_sum> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}

//[[Rcpp::export]]
Rcpp::List cppad_expm1_prod(const NumericVector& X) {
  run_test<expm1_prod> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}


//[[Rcpp::export]]
Rcpp::List cppad_lgamma_sum(const NumericVector& X) {
  run_test<lgamma_sum> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}

//[[Rcpp::export]]
Rcpp::List cppad_lgamma_prod(const NumericVector& X) {
  run_test<lgamma_prod> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}

//[[Rcpp::export]]
Rcpp::List cppad_log1pexp_prod(const NumericVector& X) {
  run_test<log1pexp_prod> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}

//[[Rcpp::export]]
Rcpp::List cppad_log1pexp_sum(const NumericVector& X) {
  run_test<log1pexp_sum> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}


//[[Rcpp::export]]
Rcpp::List cppad_invlogit_prod(const NumericVector& X) {
  run_test<invlogit_prod> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}

//[[Rcpp::export]]
Rcpp::List cppad_invlogit_sum(const NumericVector& X) {
  run_test<invlogit_sum> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}


//[[Rcpp::export]]
Rcpp::List cppad_loginvlogit_prod(const NumericVector& X) {
  run_test<loginvlogit_prod> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}

//[[Rcpp::export]]
Rcpp::List cppad_loginvlogit_sum(const NumericVector& X) {
  run_test<loginvlogit_sum> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}

//[[Rcpp::export]]
Rcpp::List cppad_lgammaexp_sum(const NumericVector& X) {
  run_test<lgammaexp_sum> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}


//[[Rcpp::export]]
Rcpp::List cppad_lgammaexp_prod(const NumericVector& X) {
  run_test<lgammaexp_prod> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}


//[[Rcpp::export]]
Rcpp::List cppad_comp1_sum(const NumericVector& X) {
  run_test<comp1_sum> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}

//[[Rcpp::export]]
Rcpp::List cppad_comp1_prod(const NumericVector& X) {
  run_test<comp1_prod> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}

//[[Rcpp::export]]
Rcpp::List cppad_comp1b(const NumericVector& X) {
  run_test<comp1b> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}

//[[Rcpp::export]]
Rcpp::List cppad_lbeta1(const NumericVector& X) {
  run_test<lbeta1> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}


//[[Rcpp::export]]
Rcpp::List cppad_lbeta_sum(const NumericVector& X) {
  run_test<lbeta_sum> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}

//[[Rcpp::export]]
Rcpp::List cppad_inc_beta_reg1(const NumericVector& X) {
  run_test<inc_beta_reg1> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}


//[[Rcpp::export]]
Rcpp::List cppad_inc_beta_reg_sum(const NumericVector& X) {
  run_test<inc_beta_reg_sum> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}


//[[Rcpp::export]]
Rcpp::List cppad_inc_gamma_reg1(const NumericVector& X) {
  run_test<inc_gamma_reg1> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}


//[[Rcpp::export]]
Rcpp::List cppad_inc_gamma_reg_sum(const NumericVector& X) {
  run_test<inc_gamma_reg_sum> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}



//[[Rcpp::export]]
Rcpp::List cppad_dnorm_log1(const NumericVector& X) {
  run_test<dnorm_log1> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}

//[[Rcpp::export]]
Rcpp::List cppad_dnorm_log_sum(const NumericVector& X) {
  run_test<dnorm_log_sum> test;
  test.record_tape(X);
  Rcpp::List res = test.cppad_results();
  return(res);
}
*/
