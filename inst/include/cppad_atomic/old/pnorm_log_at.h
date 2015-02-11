#ifndef __PNORM_LOG_AT
#define __PNORM_LOG_AT

#include <cppad_atomic/mb_atomic.h>

// NOTE:  This is a NORMALIZED incomplete beta function.
//  Equivalent to the cdf of a beta(z;a,b) distribution

using Eigen::MatrixBase;

class pnorm_log_cl {
       
 public:
  
  template<typename TX>
    void eval(const MatrixBase<TX>& x,
	      double& f) {

    const double z = x(0);
    const double m = x(1);
    const double s = x(2);
    const double R = erf((z-m)*M_SQRT1_2/s);
    
    f = log1p(R) - M_LN2;
  }

  template<typename TX, typename TD>
    void eval(const MatrixBase<TX>& x, 
	      double& f,
	      const MatrixBase<TD>& df_) {
    
    MatrixBase<TD>& df = const_cast<MatrixBase<TD>&>(df_);    
    
    const double z = x(0);
    const double m = x(1);
    const double s = x(2);
    const double R = erf((z-m)*M_SQRT1_2/s);
    const double d = R::dnorm(z, m, s, 1); // not logged

    f = log1p(R) - M_LN2;
    df(0) = exp(d-f);
    df(1) = -df(0);
    df(2) = df(0) * (m-z) / s;
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
    const double R = erf((z-m)*M_SQRT1_2/s);
    const double d = R::dnorm(z, m, s, 1); // not logged

    f = log1p(R) - M_LN2;
    df(0) = exp(d-f);
    df(1) = -df(0);
    df(2) = df(0) * (m-z) / s;
    const double s2inv = 1/(s*s);
    
    hess(0,0) = -df(0)*df(0) + df(2)/s;
    hess(1,1) = hess(0,0);
    hess(2,2) = df(2)*(z-m)*(z-m)*s2inv/s - df(2)*df(2);
    hess(2,2) -= 2*df(2)/s;
    hess(0,1) = df(0)*df(0) - df(2)/s;
    hess(1,0) = hess(0,1);
    hess(0,2) = -df(0)/s - df(2)*(z-m)*s2inv - df(0)*df(2);
    hess(2,0) = hess(0,2);
    hess(1,2) = df(0)/s + df(2)*(z-m)*s2inv + df(0)*df(2);
    hess(2,1) = hess(1,2);
  }
};

#endif
