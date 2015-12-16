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

  // W is lower triangular, filled in by rows
  
  MatrixBase<TW>& W = const_cast<MatrixBase<TW>& >(W_);
  
  const int K = W.rows();
  MatrixXA Z(K,K);
  AScalar logjac = 0;

  Z.setZero();
  W.setZero();
  W(0,0)=1;
  int idx=0;
  for (int i=1; i<K; i++) {
    AScalar sum_sqs = 0;
    for (int j=0; j<i; j++) { 
      Z(i,j) = tanh(Y(idx++));
      logjac += log(1-pow(Z(i,j),2));
      W(i,j) = Z(i,j)*sqrt(1-sum_sqs);
      logjac += 0.5 * log(1-sum_sqs);
      sum_sqs += W(i,j)*W(i,j);
    }
    W(i,i) = sqrt(1-sum_sqs);
  }
  return(logjac);

}


template<typename TY>
AScalar lkj_chol_logpdf(const MatrixBase<TY>& L,
			const AScalar& eta) {


  // Log pdf of LKJ distribution, after passing in a Cholesky
  // matrix.

  int K = L.rows();
  VectorXA v(K-1);
  AScalar c = lkj_const(eta, K);
  VectorXA log_diags = L.diagonal().tail(K-1).array().log();
 
  for (int i=0; i<K-1; i++) {
    v(i) = (K-i+2*eta-4)*log_diags(i);
  }
    
  AScalar res = c + v.sum();
  return(res);

}

template<typename TY>
AScalar lkj_logpdf(const MatrixBase<TY>& Y,
		   const AScalar& eta,
		   const int& K) {

  // Log pdf of cholesky LKJ distribution, after passing in an unconstrained
  // vector.  All Jacobians are accounted for here.
  
  MatrixXA L = MatrixXA::Zero(K,K);
  AScalar logjac = lkj_unwrap(Y, L);
  AScalar lkj_chol = lkj_chol_logpdf(L, eta);
    
  AScalar res = lkj_chol + logjac;
  return(res);

}

#endif

