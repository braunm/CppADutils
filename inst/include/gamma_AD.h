#ifndef __AD_GAMMA
#define __AD_GAMMA

#include <utilfuncs.h>

using Eigen::MatrixBase;
using Eigen::Lower;
using Eigen::Dynamic;
using Eigen::Matrix;
using Eigen::Upper;
using Eigen::Lower;

typedef CppAD::AD<double> AScalar;
typedef Matrix<AScalar, Dynamic, Dynamic> MatrixXA;
typedef Matrix<AScalar, Dynamic, 1> VectorXA;
typedef Matrix<AScalar, 1, Dynamic> RowVectorXA;


template<typename TX, typename TR, typename TA, typename TQ>
  void gamma_logpdf(const MatrixBase<TX>& X,
		    const MatrixBase<TR>& R,
		    const MatrixBase<TA>& A,
		    const MatrixBase<TQ>& out_) {
  
  MatrixBase<TQ>& out = const_cast<MatrixBase<TQ>& >(out_);

  const int k = X.size();
  VectorXA lgammaR(k);
  VectorXA res(k);
  
  for (int i=0; i<k; i++) {
    out(i) = -lgamma(R(i));
  }

  out.array() += R.array()*A.array().log() + (R.array()-1) * X.array().log() - A.array()*X.array();
}

#endif

