#ifndef __INCGAMMA_AT
#define __INCGAMMA_AT

#include <cppad_atomic/mb_atomic.h>
using Eigen::MatrixBase;

// NOTE:  This is a NORMALIZED incomplete gamma function.
//  Equivalent to the cdf of a gamma(z;r,1) distribution
        
class incgamma_cl {

 private:

  // extern "C" int digami(double*, double*, double*, double*,
  //			double*, double*, double*, double*,
  //			double*, int*);

#include "digami.h"
    
 public:
  
  template<typename TX>
    void eval(const MatrixBase<TX>& x,
	      double& f) {
    f = R::pgamma(x(0), x(1), 1, 1, 0);	
  }
  
  template<typename TX, typename TD>
    void eval(const MatrixBase<TX>& x, 
	      double& f,
	      const MatrixBase<TD>& df_) {
    
    MatrixBase<TD>& df = const_cast<MatrixBase<TD>&>(df_);
    
    double z = x(0);
    double r = x(1);	
    double D[6]; // output	
    double grlog = R::lgammafn(r);
    double gr1log = log(r) + grlog;
    double psir = R::digamma(r);
    double psir1 = psir + 1/r;
    double psidr = R::trigamma(r);
    double psidr1 = psidr - 1/(r*r);
    int ifault;
    
    digami(D, &z, &r, &grlog, &gr1log, &psir,
	   &psir1, &psidr, &psidr1, &ifault);
    
    f = D[5];	
    df(0) = exp((r-1)*log(z)-z-grlog);
    df(1) = D[2];
  }
  
  template<typename TX, typename TD, typename TH>
    void eval(const MatrixBase<TX>& x, 
	      double& f,
	      const MatrixBase<TD>& df_,
	      const MatrixBase<TH>& hess_) {
    
    MatrixBase<TD>& df = const_cast<MatrixBase<TD>&>(df_);
    MatrixBase<TH>& hess = const_cast<MatrixBase<TH>&>(hess_);
    
    double z = x(0);
    double r = x(1);	
    double log_z = log(z);
    double D[6]; // output	
    double grlog = R::lgammafn(r);
    double gr1log = log(r) + grlog;
    double psir = R::digamma(r);
    double psir1 = psir + 1/r;
    double psidr = R::trigamma(r);
    double psidr1 = psidr - 1/(r*r);
    int ifault;

    digami(D, &z, &r, &grlog, &gr1log, &psir,
	   &psir1, &psidr, &psidr1, &ifault);
    
    f = D[5];	
    df(0) = exp((r-1)*log(z)-z-grlog);
    df(1) = D[2];		
    hess(0,0) = exp((r-2)*log_z - z-grlog)*(r-z-1);
    hess(1,1) = D[3];
    hess(0,1) = exp((r-1)*log_z-z-grlog)*(log_z-psir);
    hess(1,0) = hess(0,1);
  }
}; // end class
#endif
