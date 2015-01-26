#ifndef __LOGINVLOGIT_AT
#define __LOGINVLOGIT_AT

#include <cppad_atomic/mb_atomic.h>

using Eigen::MatrixBase;

class loginvlogit_cl {
       
    public:

  // eval function value
  template<typename TX>
    void eval(const MatrixBase<TX>& x,
	      double& f) {    	
    f = x(0) -  R::log1pexp(x(0));
  }

  // eval function value and gradient
  template<typename TX, typename TD>
    void eval(const MatrixBase<TX>& x, 
	      double& f,
	      const MatrixBase<TD>& df_) {
    
    MatrixBase<TD>& df = const_cast<MatrixBase<TD>&>(df_);

    const double tmp = R::log1pexp(x(0));
    f = x(0) - tmp;
    df(0) = exp(-tmp);
  }

  // eval function value, gradient and hessian
  template<typename TX, typename TD, typename TH>
    void eval(const MatrixBase<TX>& x, 
	       double& f,
	       const MatrixBase<TD>& df_,
	      const MatrixBase<TH>& hess_) {
   
    MatrixBase<TD>& df = const_cast<MatrixBase<TD>&>(df_);
    MatrixBase<TH>& hess = const_cast<MatrixBase<TH>&>(hess_);

    const double tmp = R::log1pexp(x(0));
    f = x(0) - tmp;
    df(0) = exp(-tmp);
    hess(0,0) = -df(0) * (1-df(0));
  }
}; // end class
  
#endif
