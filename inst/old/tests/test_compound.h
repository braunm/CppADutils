
#include <tests/cppad_test_define.h>
#include <cppad_atomic/lgamma1p_atomic.h>
#include <cppad_atomic/expm1_atomic.h>
#include <cppad_atomic/log1p_atomic.h>
#include <cppad_atomic/invlogit_atomic.h>

using CppAD::Atomic::lgamma1p;
using CppAD::Atomic::expm1;
using CppAD::Atomic::log1p;
using CppAD::Atomic::invlogit;

struct comp1_sum {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {    
    int nvars = Y.size();
    AScalar res = 0;
    for (int i=0; i < nvars-1; i++) {
      res += expm1(lgamma1p(Y(i)));
    }
    res += log1p(Y(0)*Y(nvars-1)*Y(nvars-1));
    return res;
  }
};

struct comp1_prod {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {    
    int nvars = Y.size();
    AScalar res = 1;
    for (int i=0; i < nvars-1; i++) {
      res *= expm1(lgamma1p(Y(i)));
    }
    res += log1p(Y(nvars-1)*Y(nvars-1));
    return res;
  }
};

struct comp1b {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {    
    int nvars = Y.size();
    AScalar res = 0;
    AScalar a = 0.4;
    for (int i=0; i < nvars-1; i++) {
      res += expm1(lgamma1p(Y(i)));
    }
    res += log1p(Y(0)*Y(nvars-1)*Y(nvars-1));
    res *= invlogit(a);
    return res;
  }
};

