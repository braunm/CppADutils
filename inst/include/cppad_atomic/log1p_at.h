#ifndef __LOG1P_AT
#define __LOG1P_AT

#include <cppad_atomic/mb_atomic.h>

using Eigen::MatrixBase;

class log1p_cl {
       
    public:

  // eval function value
  template<typename TX>
    void eval(const MatrixBase<TX>& x,
	      double& f) {    	
    f = log1p(x(0));
  }

  // eval function value and gradient
  template<typename TX, typename TD>
    void eval(const MatrixBase<TX>& x, 
	      double& f,
	      const MatrixBase<TD>& df_) {
    
    MatrixBase<TD>& df = const_cast<MatrixBase<TD>&>(df_);
    
    f = log1p(x(0));			
    df(0) = 1/(1.0+x(0));
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
    
    f = log1p(x(0));
    df(0) = 1/(1+x(0));
    hess(0,0) = -df(0)/(1+x(0));
  }
}; // end class
  


    

#endif
