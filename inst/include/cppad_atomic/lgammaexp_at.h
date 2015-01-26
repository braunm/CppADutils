#ifndef __LGAMMAEXP_AT
#define __LGAMMAEXP_AT

#include <cppad_atomic/mb_atomic.h>

using Eigen::MatrixBase;

class lgammaexp_cl {
       
    public:

  // eval function value
  template<typename TX>
    void eval(const MatrixBase<TX>& x,
	      double& f) {
    f = R::lgamma1p(exp(x(0))) - x(0);
 
  }

  // eval function value and gradient
  template<typename TX, typename TD>
    void eval(const MatrixBase<TX>& x, 
	      double& f,
	      const MatrixBase<TD>& df_) {
    
    MatrixBase<TD>& df = const_cast<MatrixBase<TD>&>(df_);

    const double expx = exp(x(0));
    f = R::lgamma1p(expx) - x(0);
    df(0) = expx * R::digamma(1 + expx) - 1;
    
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

    const double expx = exp(x(0));
    f = R::lgamma1p(expx) - x(0);
    df(0) = expx * R::digamma(1 + expx) - 1;
    hess(0,0) = expx * (R::digamma(1 + expx) + expx*R::trigamma(1 + expx));

  }
}; // end class
  

// Function to call AD version of atomic function


#endif
