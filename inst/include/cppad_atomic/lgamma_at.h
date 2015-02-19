#ifndef __LGAMMA_AT
#define __LGAMMA_AT

#include <cppad_atomic/mb_atomic.h>

using Eigen::MatrixBase;

class lgamma_cl {
       
    public:

  // eval function value
  template<typename TX>
    void eval(const MatrixBase<TX>& x,
	      double& f) {    	
    f = R::lgammafn(x(0));
  }

  // eval function value and gradient
  template<typename TX, typename TD>
    void eval(const MatrixBase<TX>& x, 
	      double& f,
	      const MatrixBase<TD>& df_) {
    
    MatrixBase<TD>& df = const_cast<MatrixBase<TD>&>(df_);
    
    f = R::lgammafn(x(0));			
    df(0) = R::digamma(x(0));
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
    
    f = R::lgammafn(x(0));			      
    df(0) = R::digamma(x(0));    
    hess(0,0) = R::trigamma(x(0));
  }
}; // end class

#endif
