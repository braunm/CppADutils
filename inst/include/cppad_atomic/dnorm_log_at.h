#ifndef __DNORM_LOG_AT
#define __DNORM_LOG_AT

#include <cppad_atomic/mb_atomic.h>

// NOTE:  This is a NORMALIZED incomplete beta function.
//  Equivalent to the cdf of a beta(z;a,b) distribution

using Eigen::MatrixBase;

class dnorm_log_cl {
       
 public:
  
  template<typename TX>
    void eval(const MatrixBase<TX>& x,
	      double& f) {
    
    f = R::dnorm(x(0), x(1), x(2), 1);   
  }

  template<typename TX, typename TD>
    void eval(const MatrixBase<TX>& x, 
	      double& f,
	      const MatrixBase<TD>& df_) {
    
    MatrixBase<TD>& df = const_cast<MatrixBase<TD>&>(df_);    
    
    const double y = x(0);
    const double m = x(1);
    const double s = x(2);
    
    f = R::dnorm(y, m, s, 1);
    df(0) = (m-y) / (s*s);
    df(1) = -df(0);
    df(2) = (y-m-s) * (y-m+s) / (s*s*s);  
  }
  
  
  template<typename TX, typename TD, typename TH>
    void eval(const MatrixBase<TX>& x, 
	      double& f,
	      const MatrixBase<TD>& df_,
	      const MatrixBase<TH>& hess_) {
    
    MatrixBase<TD>& df = const_cast<MatrixBase<TD>&>(df_);
    MatrixBase<TH>& hess = const_cast<MatrixBase<TH>&>(hess_);

    const double y = x(0);
    const double m = x(1);
    const double s = x(2);
    
    f = R::dnorm(y, m, s, 1);
    df(0) = (m-y) / (s*s);
    df(1) = -df(0);
    df(2) = (y-m-s) * (y-m+s) / (s*s*s);
    hess(1,1) = -1/(s*s);
    hess(0,0) = hess(1,1);
    hess(0,1) = -hess(1,1);
    hess(1,0) = -hess(1,1);
    hess(0,2) = 2*df(1)/s;
    hess(2,0) = hess(0,2);
    hess(1,2) = -hess(0,2);
    hess(2,1) = hess(1,2);
    hess(2,2) = hess(0,1) - 3*df(0)*df(0);
 
  }
};

#endif
