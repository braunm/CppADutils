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
  
  AScalar t1 = 0;
  AScalar t2 = 0;
  for (int i=1; i<=(K-1); i++) {
    t1 += (2.0 * eta - 2.0 + K - i) * (K - i);
    AScalar tmp = eta + 0.5 * (K - i - 1.0);
    t2 += (K-i) * (2.0 * lgamma(tmp) - lgamma(2.0 * tmp));
  }
  AScalar c = t1 * M_LN2 + t2;
  return(c);
}


// Converts unconstrained vector to Cholesky of correlation matrix
template<typename TY, typename TW>
  void lkj_unwrap(const MatrixBase<TY>& Y,
		  const MatrixBase<TW>& W_) {

  MatrixBase<TW>& W = const_cast<MatrixBase<TW>& >(W_);
  
  const int K = W.rows();
  int idx = 0;
  MatrixXA Z = MatrixXA::Zero(K,K);
  
  W.setZero();
  W(0,0)=1;
  for (int i=1; i<K; i++) {    
    Z(i,0) = tanh(Y(idx++)); 
    W(i,0) = Z(i,0);
  }
  for (int j=1; j<K; j++) {
    W(j,j) = sqrt(1-W.block(j,0,1,j).cwiseAbs2().sum());
    for (int i=j+1; i<K; i++) {
      Z(i,j) = tanh(Y(idx++));
      W(i,j) = Z(i,j)*sqrt(1-W.block(i,0,1,j).cwiseAbs2().sum());
    }
  } 
}



template<typename TY>
AScalar lkj_logpdf(const MatrixBase<TY>& Y,
		   const AScalar& eta,
		   const int& K) {

  // Log pdf of LKJ distribution, after passing in an unconstrained
  // vector.  All Jacobians are accounted for here.
  
  int Kch2 = Y.size();
  MatrixXA L = MatrixXA::Zero(K,K);
  lkj_unwrap(Y, L);

  AScalar lkjChol = 0;

  for (int i=1; i<K; i++) {
    lkjChol += (K-i+2*eta-3)*L(i,i);
  }



  AScalar logjac1;

  for (int i=0; i<Kch2; i++) {
    logjac1 += -2*log(cosh(Y(i)));
  }

  AScalar logjac2 = 0;
  for (int j=0; j<K-1; j++) {
    for (int i=j+1; i<K; i++) {
      logjac2 += log(1-L.block(i,0,1,j).cwiseAbs2().sum());
    }
  }


  AScalar c = lkj_const(eta, K);
  AScalar a = (eta-1)*lkjChol;
  AScalar res = c + a + logjac1 + 0.5*logjac2;
  return(res);

}





#endif

