#ifndef __INVLOGIT_AT
#define __INVLOGIT_AT

#include <cppad_atomic/mb_atomic.h>

using Eigen::MatrixBase;

class invlogit_cl {
       
    public:

  // eval function value
  template<typename TX>
    void eval(const MatrixBase<TX>& x,
	      double& f) {    	
    f = exp(x(0) -  R::log1pexp(x(0)));
  }

  // eval function value and gradient
  template<typename TX, typename TD>
    void eval(const MatrixBase<TX>& x, 
	      double& f,
	      const MatrixBase<TD>& df_) {
    
    MatrixBase<TD>& df = const_cast<MatrixBase<TD>&>(df_);

    f = exp(x(0) -  R::log1pexp(x(0)));
    df(0) = f * (1-f);
  }

  // eval function value, gradient and hessian
  template<typename TX, typename TD, typename TH>
    void eval(const MatrixBase<TX>& x, 
	       double& f,
	       const MatrixBase<TD>& df_,
	      const MatrixBase<TH>& hess_) {
   
    MatrixBase<TD>& df = const_cast<MatrixBase<TD>&>(df_);
    MatrixBase<TH>& hess = const_cast<MatrixBase<TH>&>(hess_);

    f = exp(x(0) -  R::log1pexp(x(0)));
    df(0) = f * (1-f);
    hess(0,0) = -df(0) * (1-f) * expm1(x(0));
  }
}; // end class
  

// Function to call AD version of atomic function
 

#endif
