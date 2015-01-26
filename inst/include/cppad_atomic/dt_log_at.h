#ifndef __DT_LOG_AT
#define __DT_LOG_AT

#include <cppad_atomic/mb_atomic.h>

// NOTE:  This is a NORMALIZED incomplete beta function.
//  Equivalent to the cdf of a beta(z;a,b) distribution

using Eigen::MatrixBase;
using R::lgammafn;
using R::digamma;
using R::trigamma;

class dt_log_cl {

  
 public:
  
  template<typename TX>
    void eval(const MatrixBase<TX>& x,
	      double& f) {

    const double& z = x(0);
    const double& v = x(1);
    const double& s = x(2);

    f = lgammafn((v+1)/2) - lgammafn(v/2);
    f -= 0.5 * ((v+1) * log1p(z*z/(s*s*v)) + log(v));
    f -= log(s) + M_LN_SQRT_PI;
  }

  template<typename TX, typename TD>
    void eval(const MatrixBase<TX>& x, 
	      double& f,
	      const MatrixBase<TD>& df_) {
    
    MatrixBase<TD>& df = const_cast<MatrixBase<TD>&>(df_);    

    const double& z = x(0);
    const double& v = x(1);
    const double& s = x(2);
    const double zz = z*z;
    const double ssv = s*s*v;
    
    eval(x, f);
    df(0) = -z*(v+1) / (zz + ssv);
    
    df(1) = digamma((v+1)/2) - digamma(v/2);
    df(1) -= log1p(zz/ssv) + 1/v - zz*(v+1) / (v*(ssv+zz));
    df(1) *= 0.5;
    df(2) = (zz*v - ssv) / (s*(ssv + zz));
  }
  
  
  template<typename TX, typename TD, typename TH>
    void eval(const MatrixBase<TX>& x, 
	      double& f,
	      const MatrixBase<TD>& df_,
	      const MatrixBase<TH>& hess_) {
    
    MatrixBase<TD>& df = const_cast<MatrixBase<TD>&>(df_);
    MatrixBase<TH>& hess = const_cast<MatrixBase<TH>&>(hess_);

    
    const double& z = x(0);
    const double& v = x(1);
    const double& s = x(2);
    const double zz = z*z;
    const double ssv = s*s*v;
    const double s2vz22 = (ssv + zz) * (ssv + zz);
    const double s4vz4 = s*s*s*s*v + zz*zz;

    eval(x, f);
    df(0) = z*(v+1) / (zz + ssv);
    df(1) = df(0)/v - log1p(zz/(ssv)) - 1/v;
    df(1) += digamma((v+1)/2) - digamma(v/2);
    df(1) *= 0.5;
    df(2) = (zz*v - ssv) / (s*(ssv + zz));


    hess(0,0) = (v+1) * (zz - ssv) / ((zz + ssv) * (zz + ssv));
    hess(1,1) = 2*s4vz4 + s2vz22 * v * (trigamma((v+1)/2) - trigamma(v/2));
    hess(1,1) /= 4*v*s2vz22;
    hess(2,2) = (ssv*ssv - ssv*zz*(3*v+1) - v*zz*zz) / (s*s*s2vz22);
    hess(0,1) = z * (s-z) * (s+z) / s2vz22;
    hess(1,0) = hess(0,1);
    hess(0,2) = 2*s*v*z*(v+1) / s2vz22;
    hess(2,0) = hess(0,2);
    hess(1,2) = zz * (zz - s*s) / (s * s2vz22);
    hess(2,1) = hess(1,2);
 
  }
};

    

#endif
