#ifndef __DLOGITBETA_LOG_AT
#define __DLOGITBETA_LOG_AT

#include <cppad_atomic/mb_atomic.h>

// NOTE:  This is a NORMALIZED incomplete beta function.
//  Equivalent to the cdf of a beta(z;a,b) distribution

using Eigen::MatrixBase;
using R::digamma;
using R::trigamma;

class dlogitbeta_log_cl {
       
 public:
  
  template<typename TX>
    void eval(const MatrixBase<TX>& x,
	      double& f) {

    const double p = exp(x(0) - R::log1pexp(x(0)));
    const double a = x(1);
    const double b = x(2);
    
    f = R::dbeta(p, a, b, 1) + log(p) + log(1-p);
  }

  template<typename TX, typename TD>
    void eval(const MatrixBase<TX>& x, 
	      double& f,
	      const MatrixBase<TD>& df_) {
    
    MatrixBase<TD>& df = const_cast<MatrixBase<TD>&>(df_);    

    const double p = exp(x(0) - R::log1pexp(x(0)));
    const double a = x(1);
    const double b = x(2);
    const double psi_ab = digamma(a+b);

    f = R::dbeta(p, a, b, 1) + log(p) + log(1-p);

    df(0) = (a+b)*(1-p)-b;
    df(1) = log(p) - digamma(a) + psi_ab;
    df(2) = log(1-p)- digamma(b) + psi_ab;
  }
  
  
  template<typename TX, typename TD, typename TH>
    void eval(const MatrixBase<TX>& x, 
	      double& f,
	      const MatrixBase<TD>& df_,
	      const MatrixBase<TH>& hess_) {
    
    MatrixBase<TD>& df = const_cast<MatrixBase<TD>&>(df_);
    MatrixBase<TH>& hess = const_cast<MatrixBase<TH>&>(hess_);

    const double p = exp(x(0) - R::log1pexp(x(0)));
    const double a = x(1);
    const double b = x(2);
    const double psi_ab = digamma(a+b);
    const double psi_1_ab = trigamma(a+b);
    
    f = R::dbeta(p, a, b, 1) + log(p) + log(1-p);
    df(0) = (a+b)*(1-p)-b;
    df(1) = log(p) - digamma(a) + psi_ab;
    df(2) = log(1-p)- digamma(b) + psi_ab;

    hess(0,0) = -(a+b) * p * (1-p);
    hess(1,1) = psi_1_ab - trigamma(a);
    hess(2,2) = psi_1_ab - trigamma(b);
    hess(0,1) = 1-p;
    hess(1,0) = hess(0,1);
    hess(0,2) = -p;
    hess(2,0) = hess(0,2);
    hess(1,2) = psi_1_ab;
    hess(2,1) = hess(1,2);
  }
};

#endif
