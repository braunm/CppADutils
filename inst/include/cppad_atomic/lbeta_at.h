#ifndef __LBETA_AT
#define __LBETA_AT

#include <cppad_atomic/mb_atomic.h>
#include <RcppEigen.h>
#include <Eigen/Core>

using Eigen::MatrixBase;

class lbeta_cl {
 public:

  lbeta_cl() {
    //   Rcout << "Constructung lbeta_cl\n";
  }

  ~lbeta_cl() {
    //   Rcout << "Destructing lbeta_cl\n";
  }
  
  template<typename TX>    
    void eval(const MatrixBase<TX>& x,
	       double& f) {    	
    f = R::lbeta(x(0), x(1));	
  }
  
  template<typename TX, typename TD>
    void eval(const MatrixBase<TX>& x, 
	       double& f,
	       const MatrixBase<TD>& df_) {

    MatrixBase<TD>& df = const_cast<MatrixBase<TD>&>(df_);
    
    f = R::lbeta(x(0), x(1));			
    const double psi_a = R::digamma(x(0));
    const double psi_b = R::digamma(x(1));
    const double psi_ab = R::digamma(x.sum());
    df(0) = psi_a - psi_ab;
    df(1) = psi_b - psi_ab;
  }

  template<typename TX, typename TD, typename TH>
    void eval(const MatrixBase<TX>& x, 
	       double& f,
	       const MatrixBase<TD>& df_,
	       const MatrixBase<TH>& hess_)
    
  {
    
    MatrixBase<TD>& df = const_cast<MatrixBase<TD>&>(df_);
    MatrixBase<TH>& hess = const_cast<MatrixBase<TH>&>(hess_);
    
    f = R::lbeta(x(0), x(1));			
    
    const double psi_a = R::digamma(x(0));
    const double psi_b = R::digamma(x(1));
    const double psi_ab = R::digamma(x.sum());
    df(0) = psi_a - psi_ab;
    df(1) = psi_b - psi_ab;
    
    const double psi_1_a = R::trigamma(x(0));
    const double psi_1_b = R::trigamma(x(1));
    const double psi_1_ab = R::trigamma(x.sum());
    
    hess(0,0) = psi_1_a - psi_1_ab;
    hess(1,1) = psi_1_b - psi_1_ab;
    hess(0,1) = -psi_1_ab;
    hess(1,0) = -psi_1_ab;
  }
}; // end class

#endif
