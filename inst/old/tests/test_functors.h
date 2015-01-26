
#include <tests/cppad_test_define.h>
#include <cppad_atomic/log1p_atomic.h>
#include <cppad_atomic/lgamma_atomic.h>
#include <cppad_atomic/lgamma1p_atomic.h>
#include <cppad_atomic/log1pexp_atomic.h>
#include <cppad_atomic/invlogit_atomic.h>
#include <cppad_atomic/loginvlogit_atomic.h>
#include <cppad_atomic/expm1_atomic.h>
#include <cppad_atomic/log1pmx_atomic.h>
#include <cppad_atomic/lgammaexp_atomic.h>
#include <cppad_atomic/lbeta_atomic.h>
#include <cppad_atomic/inc_beta_atomic.h>
#include <cppad_atomic/inc_gamma_atomic.h>
#include <cppad_atomic/dnorm_log_atomic.h>

struct log1p_sum {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {    
    int nvars = Y.size();
    AScalar res = 0;
    for (int i=0; i<nvars; i++) {
      res += CppAD::Atomic::log1p(Y(i));
    }
    return res;
  }
};

struct log1pmx_prod {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {    
    int nvars = Y.size();
    AScalar res = 1;
    for (int i=0; i<nvars; i++) {
      res *= CppAD::Atomic::log1pmx(Y(i));
    }
    return res;
  }
};

struct log1pmx_sum {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {    
    int nvars = Y.size();
    AScalar res = 0;
    for (int i=0; i<nvars; i++) {
      res += CppAD::Atomic::log1pmx(Y(i));
    }
    return res;
  }
};

struct log1p_prod {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {    
    int nvars = Y.size();
    AScalar res = 1;
    for (int i=0; i<nvars; i++) {
      res *= CppAD::Atomic::log1p(Y(i));
    }
    return res;
  }
};


struct expm1_sum {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {    
    int nvars = Y.size();
    AScalar res = 0;
    for (int i=0; i<nvars; i++) {
      res += CppAD::Atomic::expm1(Y(i));
    }
    return res;
  }
};

struct expm1_prod {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {    
    int nvars = Y.size();
    AScalar res = 1;
    for (int i=0; i<nvars; i++) {
      res *= CppAD::Atomic::expm1(Y(i));
    }
    return res;
  }
};



struct lgammaexp_sum {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {    
    int nvars = Y.size();
    AScalar res = 0;
    for (int i=0; i<nvars; i++) {
      res += CppAD::Atomic::lgammaexp(Y(i));
    }
    return res;
  }
};


struct lgammaexp_prod {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {    
    int nvars = Y.size();
    AScalar res = 1;
    for (int i=0; i<nvars; i++) {
      res *= CppAD::Atomic::lgammaexp(Y(i));
    }
    return res;
  }
};


struct lgamma1p_sum {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {    
    int nvars = Y.size();
    AScalar res = 0;
    for (int i=0; i<nvars; i++) {
      res += CppAD::Atomic::lgamma1p(Y(i));
    }
    return res;
  }
};


struct lgamma1p_prod {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {    
    int nvars = Y.size();
    AScalar res = 1;
    for (int i=0; i<nvars; i++) {
      res *= CppAD::Atomic::lgamma1p(Y(i));
    }
    return res;
  }
};

struct log1pexp_sum {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {    
    int nvars = Y.size();
    AScalar res = 0;
    for (int i=0; i<nvars; i++) {
      res += CppAD::Atomic::log1pexp(Y(i));
    }
    return res;
  }
};


struct log1pexp_prod {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {    
    int nvars = Y.size();
    AScalar res = 1;
    for (int i=0; i<nvars; i++) {
      res *= CppAD::Atomic::log1pexp(Y(i));
    }
    return res;
  }
};


struct invlogit_sum {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {    
    int nvars = Y.size();
    AScalar res = 0;
    for (int i=0; i<nvars; i++) {
      res += CppAD::Atomic::invlogit(Y(i));
    }
    return res;
  }
};

struct invlogit_prod {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {    
    int nvars = Y.size();
    AScalar res = 1;
    for (int i=0; i<nvars; i++) {
      res *= CppAD::Atomic::invlogit(Y(i));
    }
    return res;
  }
};


struct loginvlogit_sum {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {    
    int nvars = Y.size();
    AScalar res = 0;
    for (int i=0; i<nvars; i++) {
      res += CppAD::Atomic::loginvlogit(Y(i));
    }
    return res;
  }
};

struct loginvlogit_prod {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {    
    int nvars = Y.size();
    AScalar res = 1;
    for (int i=0; i<nvars; i++) {
      res *= CppAD::Atomic::loginvlogit(Y(i));
    }
    return res;
  }
};

/* struct lbeta1 {   */
/*   template<typename TA>   */
/*   AScalar eval(const MatrixBase<TA>& Y) { */
/*     AScalar res = CppAD::Atomic::lbeta(Y(0), Y(1)); */
/*     return res; */
/*   } */
/* }; */

/* struct lbeta_sum {   */
/*   template<typename TA>   */
/*   AScalar eval(const MatrixBase<TA>& Y) { */
/*     size_t n = Y.size(); // must be even */
/*     AScalar res = 0; */
/*     for (size_t i=0; i<n/2; i++) { */
/*       res += CppAD::Atomic::lbeta(Y(2*i), Y(2*i+1)); */
/*     } */
/*     return res; */
/*   } */
/* }; */


struct inc_beta_reg1 {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {
    AScalar res = CppAD::Atomic::inc_beta(Y(0), Y(1), Y(2));
    return res;
  }
};

struct inc_beta_reg_sum {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {
    size_t n = Y.size(); // must be divisible by 3
    AScalar res = 0;
    for (size_t i=0; i < n/3; i++) {
      res += CppAD::Atomic::inc_beta(Y(3*i), Y(3*i+1), Y(3*i+2));
    }
    return res;
  }
};

struct inc_gamma_reg1 {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {
    AScalar res = CppAD::Atomic::inc_gamma(Y(0), Y(1));
    return res;
  }
};

struct inc_gamma_reg_sum {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {
    size_t n = Y.size(); // must be divisible by 2
    AScalar res = 0;
    for (size_t i=0; i < n/2; i++) {
      res += CppAD::Atomic::inc_gamma(Y(2*i), Y(2*i+1));
    }
    return res;
  }
};

struct dnorm_log1 {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {
    AScalar res = CppAD::Atomic::dnorm_log(Y(0), Y(1), Y(2));
    return res;
  }
};

struct dnorm_log_sum {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {
    size_t n = Y.size(); // must be divisible by 3
    AScalar res = 0;
    for (size_t i=0; i < n/3; i++) {
      res += CppAD::Atomic::dnorm_log(Y(3*i), Y(3*i+1), Y(3*i+2));
    }
    return res;
  }
};

