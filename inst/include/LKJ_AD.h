#ifndef __AD_LKJ
#define __AD_LKJ

#include <cppad_atomics.h>

using Eigen::MatrixBase;
using Eigen::Lower;
using Eigen::Dynamic;
using Eigen::Matrix;


typedef CppAD::AD<double> AScalar;
typedef Matrix<AScalar, Dynamic, Dynamic> MatrixXA;
typedef Matrix<AScalar, Dynamic, 1> VectorXA;
typedef Matrix<AScalar, 1, Dynamic> RowVectorXA;


AScalar lkj_const(const AScalar& eta,
		  const int& K) {

  // log of constant in LKJ,  Eq. 16

  AScalar res = 0;
  AScalar s = 0;
  for (int i=1; i<K; i++) {
    AScalar diff = K-i;
    AScalar barg = eta + 0.5*(diff-1);
    res += diff*(2*lgamma(barg)-lgamma(2*barg));
    s += diff*(2*eta-2+diff);
  }
  res += s*M_LN2;
  return(res);
}


// Converts unconstrained vector to Cholesky of correlation matrix
template<typename TY, typename TW>
  AScalar lkj_unwrap(const MatrixBase<TY>& Y,
		     const MatrixBase<TW>& W_) {

  MatrixBase<TW>& W = const_cast<MatrixBase<TW>& >(W_);
  
  const int K = W.rows();
  const int Kch2 = R::choose(K,2);
  VectorXA Z(Kch2);
  AScalar logjac = 0;
  
  for (int i=0; i<Kch2; i++) {
    Z(i) = tanh(Y(i));
    //  logjac += log1p(-Z(i)*Z(i));
    logjac += -2*log(cosh(Y(i)));
  }
  
  W.setZero();
  W(0,0)=1;
  int idx=0;
  for (int i=1; i<K; ++i) {
    W(i,0) = Z(idx++);
    AScalar sum_sqs = W(i,0)*W(i,0);
    for (int j=1; j<i; ++j) {
      logjac += 0.5*log1p(-sum_sqs);
      W(i,j) = Z(idx++) * sqrt(1-sum_sqs);
      sum_sqs += W(i,j)*W(i,j);
    }
    W(i,i) = sqrt(1-sum_sqs);
  }
  
  return(logjac);

}



template<typename TY>
AScalar lkj_logpdf(const MatrixBase<TY>& Y,
		   const AScalar& eta,
		   const int& K) {

  // Log pdf of LKJ distribution, after passing in an unconstrained
  // vector.  All Jacobians are accounted for here.
  
  MatrixXA L = MatrixXA::Zero(K,K);
  AScalar logjac = lkj_unwrap(Y, L);
  VectorXA v(K-1);


  AScalar c = lkj_const(eta, K);
  for (int i=0; i<K-1; i++) {
    //   v(i) = (K-i+2*eta-4)*log(L(i+1,i+1));
    v(i) += (eta-1)*2*log(L(i+1,i+1));
  }
    
  AScalar res = c + v.sum() + logjac;
  return(res);

}





#endif

