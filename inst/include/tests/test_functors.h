#include <tests/cppad_test_define.h>
#include <cppad_atomic/atomics.h>

struct lbeta1 {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {
    size_t n = Y.size()/2;
    AScalar res = 0;
    for (size_t i=0; i<n; i++) {
      res += lbeta(Y(2*i), Y(2*i+1));
    }
    return res;
  }
};

struct lgamma1 {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {
    size_t n = Y.size()/2;
    AScalar res = 0;
    for (size_t i=0; i<n; i++) {
      res += lgamma(Y(2*i)) * lgamma(Y(2*i+1));
    }
    return res;
  }
};

struct lgamma1p_test {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {
    size_t n = Y.size()/2;
    AScalar res = 0;
    for (size_t i=0; i<n; i++) {
      res += lgamma1p(Y(2*i)) * lgamma1p(Y(2*i+1));
    }
    return res;
  }
};

struct lgammaexp_test {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {
    size_t n = Y.size()/2;
    AScalar res = 0;
    for (size_t i=0; i<n; i++) {
      res += lgamma(exp(Y(2*i))) * lgamma(exp(Y(2*i+1)));
    }
    return res;
  }
};

struct lgammaLogExp {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {
    size_t n = Y.size()/2;
    AScalar res = 0;
    for (size_t i=0; i<n; i++) {
      AScalar a = log(Y(2*i));
      AScalar b = exp(Y(2*i+1));
      res += lgamma(a) * lgamma(b);
    }
    return res;
  }
};

struct log1p_test {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {
    size_t n = Y.size()/2;
    AScalar res = 0;
    for (size_t i=0; i<n; i++) {
      res += log1p(Y(2*i)) * log1p(Y(2*i+1));
    }
    return res;
  }
};

struct expm1_test {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {
    size_t n = Y.size()/2;
    AScalar res = 0;
    for (size_t i=0; i<n; i++) {
      res += expm1(Y(2*i)) * expm1(Y(2*i+1));
    }
    return res;
  }
};

struct log1pexp_test {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {
    size_t n = Y.size()/2;
    AScalar res = 0;
    for (size_t i=0; i<n; i++) {
      res += log1pexp(Y(2*i)) * lgamma(log1pexp(Y(2*i+1)));
    }
    return res;
  }
};

struct loginvlogit_test {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {
    size_t n = Y.size()/2;
    AScalar res = 0;
    for (size_t i=0; i<n; i++) {
      res += loginvlogit(Y(2*i)) * loginvlogit(Y(2*i+1));
    }
    return res;
  }
};

struct invlogit_test {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {
    size_t n = Y.size()/2;
    AScalar res = 0;
    for (size_t i=0; i<n; i++) {
      res += invlogit(Y(2*i)) * invlogit(Y(2*i+1));
    }
    return res;
  }
};


struct log1pmx_test {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {
    size_t n = Y.size()/2;
    AScalar res = 0;
    for (size_t i=0; i<n; i++) {
      res += log1pmx(Y(2*i)) * expm1(log1pmx(Y(2*i+1)));
    }
    return res;
  }
};

struct incgamma_test {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {
    size_t n = Y.size()/2;
    AScalar res = 0.0;
    for (size_t i=0; i<n; i++) {
      res += incgamma(Y(2*i), Y(2*i+1));
    }
    return res;
  }
};

struct incbeta_test {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {
    size_t n = Y.size()/3;
    AScalar res = 0.0;
    for (size_t i=0; i<n; i++) {
      res += incbeta(Y(3*i), Y(3*i+1), Y(3*i+2));
    }
    return (res*res);
  }
};


struct incbeta_test2 {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {
    size_t n = Y.size()/3;
    AScalar res1 = 0.0;
    AScalar res2 = 0.0;
    for (size_t i=0; i<n; i++) {
      res1 += incbeta(Y(3*i), Y(3*i+1), Y(3*i+2));
      res2 += dhalft_log(Y(3*i), Y(3*i+1), Y(3*i+2));
    }
    return (res1*res2);
  }
};

struct dnorm_log_test {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {
    size_t n = Y.size()/3;
    AScalar res = 0.0;
    for (size_t i=0; i<n; i++) {
      res += dnorm_log(Y(3*i), Y(3*i+1), Y(3*i+2));
    }
    return res;
  }
};

struct dbeta_log_test {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {
    size_t n = Y.size()/3;
    AScalar res = 0.0;
    for (size_t i=0; i<n; i++) {
      res += dbeta_log(Y(3*i), Y(3*i+1), Y(3*i+2));
    }
    return res;
  }
};

struct dlogitbeta_log_test {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {
    size_t n = Y.size()/3;
    AScalar res = 0.0;
    for (size_t i=0; i<n; i++) {
      res += dlogitbeta_log(Y(3*i), Y(3*i+1), Y(3*i+2));
    }
    return res;
  }
};

struct dnormTrunc0_log_test {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {
    size_t n = Y.size()/3;
    AScalar res = 0.0;
    for (size_t i=0; i<n; i++) {
      res += dnormTrunc0_log(Y(3*i), Y(3*i+1), Y(3*i+2));
    }
    return res;
  }
};

struct dt_log_test {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {
    size_t n = Y.size()/3;
    AScalar res = 0.0;
    for (size_t i=0; i<n; i++) {
      res += dt_log(Y(3*i), Y(3*i+1), Y(3*i+2));
    }
    return res;
  }
};

struct dhalft_log_test {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {
    size_t n = Y.size()/3;
    AScalar res = 0.0;
    for (size_t i=0; i<n; i++) {
      res += dhalft_log(Y(3*i), Y(3*i+1), Y(3*i+2));
    }
    return res;
  }
};


struct dhalft_log_test2 {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {
    size_t n = Y.size()/3;
    AScalar res1 = 0.0;
    AScalar res2 = 0.0;
    for (size_t i=0; i<n; i++) {
      res1 += dhalft_log(Y(3*i), Y(3*i+1), Y(3*i+2));
      res2 += dnorm_log(Y(3*i), Y(3*i+1), Y(3*i+2));
    }
    return res1*res2;
  }
};

struct pnorm_log_test {  
  template<typename TA>  
  AScalar eval(const MatrixBase<TA>& Y) {
    size_t n = Y.size()/3;
    AScalar res = 0.0;
    for (size_t i=0; i<n; i++) {
      res += pnorm_log(Y(3*i), Y(3*i+1), Y(3*i+2));
    }
    return res;
  }
};
