#ifndef __LOG1PMX_AT
#define __LOG1PMX_AT

#include <cppad_atomic/mb_atomic.h>

using Eigen::MatrixBase;

class log1pmx_cl {
       
    public:

  // eval function value
  template<typename TX>
    void eval(const MatrixBase<TX>& x,
	      double& f) {    	
    f = R::log1pmx(x(0));
  }

  // eval function value and gradient
  template<typename TX, typename TD>
    void eval(const MatrixBase<TX>& x, 
	      double& f,
	      const MatrixBase<TD>& df_) {
    
    MatrixBase<TD>& df = const_cast<MatrixBase<TD>&>(df_);
    
    f = R::log1pmx(x(0));			
    df(0) = 1.0 /(x(0) + 1.0) - 1.0;
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

    f = R::log1pmx(x(0));			
    df(0) = 1.0 / (x(0) + 1.0) - 1.;
    hess(0,0) = -1.0 / ((x(0) + 1.0) * (x(0) + 1.0));
  }
}; // end class
  

// Function to call AD version of atomic function


#endif
