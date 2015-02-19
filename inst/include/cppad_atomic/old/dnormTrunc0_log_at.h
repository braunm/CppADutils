#ifndef __DNORMTRUNC0_LOG_AT
#define __DNORMTRUNC0_LOG_AT

#include <cppad_atomic/mb_atomic.h>

// NOTE:  This is a NORMALIZED incomplete beta function.
//  Equivalent to the cdf of a beta(z;a,b) distribution

using Eigen::MatrixBase;

class dnormTrunc0_log_cl {
       
 public:
  
  template<typename TX>
    void eval(const MatrixBase<TX>& x,
	      double& f) {

    const double z = x(0);
    const double m = x(1);
    const double s = x(2);
    
    f = R::dnorm(z, m, s, 1) - log(erfc(-M_SQRT1_2 * m/s)) + M_LN2;
  }

  template<typename TX, typename TD>
    void eval(const MatrixBase<TX>& x, 
	      double& f,
	      const MatrixBase<TD>& df_) {
    
    MatrixBase<TD>& df = const_cast<MatrixBase<TD>&>(df_);    
    
    const double z = x(0);
    const double m = x(1);
    const double s = x(2);
    const double R = erfc(-M_SQRT1_2 * m / s);

    f = R::dnorm(z, m, s, 1) - log(R) + M_LN2;
    const double expf0 = exp(R::dnorm(0,m,s,1) - log(R) + M_LN2);
    
    df(0) = (m-z) / (s*s);
    df(1) = -df(0) - expf0;
    df(2) = expf0*m/s + s*df(0)*df(0) - 1/s;
  }
  
  
  template<typename TX, typename TD, typename TH>
    void eval(const MatrixBase<TX>& x, 
	      double& f,
	      const MatrixBase<TD>& df_,
	      const MatrixBase<TH>& hess_) {
    
    MatrixBase<TD>& df = const_cast<MatrixBase<TD>&>(df_);
    MatrixBase<TH>& hess = const_cast<MatrixBase<TH>&>(hess_);

    const double z = x(0);
    const double m = x(1);
    const double s = x(2);
    const double R = erfc(-M_SQRT1_2 * m / s);

    f = R::dnorm(z, m, s, 1) - log(R) + M_LN2;
    const double expf0 = exp(R::dnorm(0,m,s,1) - log(R) + M_LN2);
    const double s2inv = 1/(s*s);
    const double m2 = m*m;
    
    df(0) = (m-z) * s2inv;
    df(1) = -df(0) - expf0;
    df(2) = (expf0*m - 1)/s + s*df(0)*df(0);

    hess(1,0) = s2inv;
    hess(0,0) = -s2inv;
    hess(0,1) = hess(1,0);
    hess(0,2) = -2.0 * df(0) / s;
    hess(2,0) = hess(0,2);
    hess(1,1) = m * expf0 * s2inv + expf0*expf0 - s2inv;
    hess(2,2) = expf0*m*s2inv*(expf0*m - 2 + m2*s2inv); 
    hess(2,2) += s2inv - 3*df(0)*df(0);
    hess(1,2) = expf0*(1- m2*s2inv - expf0*m)/s;
    hess(1,2) += 2*df(0)/s;
    hess(2,1) = hess(1,2); 
  }
};

#endif
