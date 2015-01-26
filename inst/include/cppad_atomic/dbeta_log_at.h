#ifndef __DBETA_LOG_AT
#define __DBETA_LOG_AT

#include <cppad_atomic/mb_atomic.h>

// NOTE:  This is a NORMALIZED incomplete beta function.
//  Equivalent to the cdf of a beta(z;a,b) distribution

using Eigen::MatrixBase;
using R::digamma;
using R::trigamma;

class dbeta_log_cl {
       
 public:
  
  template<typename TX>
    void eval(const MatrixBase<TX>& x,
	      double& f) {
    
    f = R::dbeta(x(0), x(1), x(2), 1);   
  }

  template<typename TX, typename TD>
    void eval(const MatrixBase<TX>& x, 
	      double& f,
	      const MatrixBase<TD>& df_) {
    
    MatrixBase<TD>& df = const_cast<MatrixBase<TD>&>(df_);    
    
    const double p = x(0);
    const double a = x(1);
    const double b = x(2);
    const double psi_ab = digamma(a+b);
    
    f = R::dbeta(p, a, b, 1);
    df(0) = (a-1)/p - (b-1)/(1-p);
    df(1) = psi_ab - digamma(a) + log(p);
    df(2) = psi_ab - digamma(b) + log(1-p);
  }
  
  
  template<typename TX, typename TD, typename TH>
    void eval(const MatrixBase<TX>& x, 
	      double& f,
	      const MatrixBase<TD>& df_,
	      const MatrixBase<TH>& hess_) {
    
    MatrixBase<TD>& df = const_cast<MatrixBase<TD>&>(df_);
    MatrixBase<TH>& hess = const_cast<MatrixBase<TH>&>(hess_);

    const double p = x(0);
    const double a = x(1);
    const double b = x(2);
    const double psi_a = digamma(a);
    const double psi_b = digamma(b);
    const double psi_ab = digamma(a+b);
    const double psi_1_ab = trigamma(a+b);
    const double p2 = p*p;
    const double pm12 = (1-p)*(1-p);
    
    f = R::dbeta(p, a, b, 1);
    df(0) = (a-1)/p - (b-1)/(1-p);
    df(1) = psi_ab - psi_a + log(p);
    df(2) = psi_ab - psi_b + log(1-p);

    hess(0,0) = (1-a)/p2 + (1-b)/pm12;
    hess(1,1) = psi_1_ab - trigamma(a);
    hess(2,2) = psi_1_ab - trigamma(b);
    hess(0,1) = 1/p;
    hess(1,0) = hess(0,1);
    hess(0,2) = 1/(p-1);
    hess(2,0) = hess(0,2);
    hess(1,2) = psi_1_ab;
    hess(2,1) = hess(1,2);
  }
};

#endif
