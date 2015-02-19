#ifndef __LOG1PEXP_AT
#define __LOG1PEXP_AT

#include <cppad_atomic/mb_atomic.h>

using Eigen::MatrixBase;

class log1pexp_cl {
       
    public:

  // eval function value
  template<typename TX>
    void eval(const MatrixBase<TX>& x,
	      double& f) {    	
    f = log1pexp(x(0));
  }

  // eval function value and gradient
  template<typename TX, typename TD>
    void eval(const MatrixBase<TX>& x, 
	      double& f,
	      const MatrixBase<TD>& df_) {
    
    MatrixBase<TD>& df = const_cast<MatrixBase<TD>&>(df_);
    
    f = R::log1pexp(x(0));			
    df(0) = exp(x(0) - f);
  }

  // eval function value, gradient and hessian
  template<typename TX, typename TD, typename TH>
    void eval(const MatrixBase<TX>& x, 
	       double& f,
	       const MatrixBase<TD>& df_,
	       const MatrixBase<TH>& hess_)
    
  {
    
    MatrixBase<TD>& df = const_cast<MatrixBase<TD>&>(df_);
    MatrixBase<TH>& hess = const_cast<MatrixBase<TH>&>(hess_);
    
    f = log1pexp(x(0));
    df(0) = exp(x(0) - f);
    hess(0,0) = exp(x(0) - 2*f);
  }
}; // end class
  

#endif
