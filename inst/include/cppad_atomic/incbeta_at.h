#ifndef __INCBETA_AT
#define __INCBETA_AT

#include <cppad_atomic/mb_atomic.h>
//#include "inc_beta_deriv.h"

// NOTE:  This is a NORMALIZED incomplete beta function.
//  Equivalent to the cdf of a beta(z;a,b) distribution

using Eigen::MatrixBase;

class incbeta_cl {

 private:

#include "inbeder.h"
  
  template <typename TG>
    int betader(double z,
		double a,
		double b,
		double& f,
		const MatrixBase<TG>& grad_) {
    
    MatrixBase<TG>& grad = const_cast<MatrixBase<TG>&>(grad_);
    
    double log_z = log(z);
    double log_1mz = log(1-z);
    
    double psi[7];
    psi[0] = R::lbeta(a, b);   // gsl_sf_lnbeta(a,b);
    psi[1] = R::digamma(a);    // gsl_sf_psi(a);
    psi[2] = R::trigamma(a);   // gsl_sf_psi_1(a);
    psi[3] = R::digamma(b);    // gsl_sf_psi(b);
    psi[4] = R::trigamma(b);   // gsl_sf_psi_1(b);
    psi[5] = R::digamma(a+b);  // gsl_sf_psi(a+b);
    psi[6] = R::trigamma(a+b); // gsl_sf_psi_1(a+b);
    f = R::pbeta(z, a, b, 1, 0);  // gsl_sf_beta_inc(a, b, z);
    
    double der[6]; // output vector
    
    int nappx;
    double errapx;
    int ifault;
    
    inbeder(&z, &a, &b, psi, der, &nappx, &errapx, &ifault);
    
    grad(0) = exp((a-1)*log_z+(b-1)*log_1mz-psi[0]);
    grad(1) = der[1];
    grad(2) = der[3];
    return (ifault);
  }
  
  
  template <typename TG, typename TH>
    int betader(double z,
		double a,
		double b,
		double& f,
		const MatrixBase<TG>& grad_,
		const MatrixBase<TH>& hess_)
  {
    
    MatrixBase<TG>& grad = const_cast<MatrixBase<TG>&>(grad_);
    MatrixBase<TH>& hess = const_cast<MatrixBase<TH>&>(hess_);
    
    double log_z = log(z);
    double log_1mz = log(1-z);
    
    double psi[7];
    psi[0] = R::lbeta(a, b);   // gsl_sf_lnbeta(a,b);
    psi[1] = R::digamma(a);    // gsl_sf_psi(a);
    psi[2] = R::trigamma(a);   // gsl_sf_psi_1(a);
    psi[3] = R::digamma(b);    // gsl_sf_psi(b);
    psi[4] = R::trigamma(b);   // gsl_sf_psi_1(b);
    psi[5] = R::digamma(a+b);  // gsl_sf_psi(a+b);
    psi[6] = R::trigamma(a+b); // gsl_sf_psi_1(a+b);
    f = R::pbeta(z, a, b, 1, 0);  // gsl_sf_beta_inc(a, b, z);
    
    double der[6]; // output vector
    
    int nappx;
    double errapx;
    int ifault;
    
    inbeder(&z, &a, &b, psi, der, &nappx, &errapx, &ifault);
    
    grad(0) = exp((a-1)*log_z+(b-1)*log_1mz-psi[0]);
    grad(1) = der[1];
    grad(2) = der[3];
    hess(0,0) = -exp((a-2)*log_z+(b-2)*log_1mz-psi[0])*(a*(z-1)+(b-2)*z+1);
    hess(0,1) = grad(0)*(psi[5]-psi[1]+log_z);
    hess(0,2) = grad(0)*(psi[5]-psi[3]+log_1mz);
    hess(1,0) = hess(0,1);
    hess(1,1) = der[2];
    hess(1,2) = der[5];
    hess(2,0) = hess(0,2);
    hess(2,1) = der[5];
    hess(2,2) = der[4];
    
    return (ifault);
  }
  
 public:
  
  template<typename TX>
    void eval(const MatrixBase<TX>& x,
	      double& f) {
    
    f = R::pbeta(x(0), x(1), x(2), 1, 0);   
  }

  template<typename TX, typename TD>
    void eval(const MatrixBase<TX>& x, 
	      double& f,
	      const MatrixBase<TD>& df_) {

    
    MatrixBase<TD>& df = const_cast<MatrixBase<TD>&>(df_);    
    betader(x(0), x(1), x(2), f, df);    
  }
  
  
  template<typename TX, typename TD, typename TH>
    void eval(const MatrixBase<TX>& x, 
	      double& f,
	      const MatrixBase<TD>& df_,
	      const MatrixBase<TH>& hess_) {
    
    MatrixBase<TD>& df = const_cast<MatrixBase<TD>&>(df_);
    MatrixBase<TH>& hess = const_cast<MatrixBase<TH>&>(hess_);
    
    betader(x(0), x(1), x(2), f, df, hess);    
  }
};

#endif
