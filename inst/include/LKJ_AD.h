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


template<typename TY, typename TZ>
  void lkj_Y_to_Z(const MatrixBase<TY>& Y,
		  const MatrixBase<TZ>& Z_) {

  MatrixBase<TZ>& Z = const_cast<MatrixBase<TZ>& >(Z_);
  const int K = Z.rows();

  int idx=0;
  for (int j=0; j<K-1; j++) {    
    for (int i=j+1; i<K; i++) {
      Z(i,j) = tanh(Y(idx++));
    }
  }
}

template<typename TZ, typename TX>
  void lkj_Z_to_X(const MatrixBase<TZ>& Z,
		  const MatrixBase<TX>& X_) {

  MatrixBase<TX>& X = const_cast<MatrixBase<TX>& >(X_);
  const int K = Z.rows();


  Z.setZero();
  Z(0,0) = 1;
  for (int i=1; i<K; i++) {
    X(i,0) = Z(i,0);
  }
  for (int j=1; j<K; j++) {
    X(j,j) = sqrt(1-X.block(j,0,1,j).cwiseAbs2().sum());
    for (int i=j+1; i<K; i++) {
      X(i,j) = Z(i,j) * sqrt(1-X.block(i,0,1,j).cwiseAbs2().sum());
    }
  }
}



template<typename TZ>
AScalar lkj_logpdf_Z(const MatrixBase<TZ>& Z,
		     const AScalar& eta) {
  
  // W - lower triangular Z matrix of transforms
  // returns log pdf of LKJ prior
  // Includes Jacobian with respect to Z
  
  const int K = Z.rows();
  
  const AScalar c = lkj_const(eta, K); // log normalizing constant

  AScalar logdet = 0;
  AScalar logjac = 0;
  for (int j = 0; j<=K-1; j++) {
    for (int i=j+1; i<K; i++) {
      AScalar t = log1p(-pow(Z(i,j),2));
      logdet += t;
      logjac += 0.5 * (K-j) * t;
    }
  }
  AScalar res = c + (eta-1) * logdet + logjac;

  return(res);

}


template<typename TY>
AScalar lkj_logpdf(const MatrixBase<TY>& Y,
		   const AScalar& eta,
		   const int& K) {

  // Log pdf of LKJ distribution, after passing in an unconstrained
  // vector.  All Jacobians are accounted for here.
  
  MatrixXA Z = MatrixXA::Zero(K,K);
  lkj_Y_to_Z(Y, Z);
  
  AScalar logjac=0;
  for (int j=0; j<K-1; j++) {    
    for (int i=j+1; i<K; i++) {
      logjac += log1p(-pow(Z(i,j),2));
    }
  }

  AScalar logpdf_Z = lkj_logpdf_Z(Z, eta);
  AScalar res = logpdf_Z + logjac;
  return(res);
}




// Converts unconstrained vector to Cholesky of correlation matrix
template<typename TY, typename TW>
  void lkj_unwrap(const MatrixBase<TY>& Y,
		  const MatrixBase<TW>& W_) {

  MatrixBase<TW>& W = const_cast<MatrixBase<TW>& >(W_);
  
  const int K = W.rows();

  MatrixXA Z = MatrixXA::Zero(K,K);
  MatrixXA X = MatrixXA::Zero(K,K);
  lkj_Y_to_Z(Y, Z);
  lkj_Z_to_X(Z, X); // Cholesky of correlation matrix
  AScalar jac1=0, jac2=0;
  /* for (int j=1; j<K; j++) { */
  /*   jac1 += 0.5 * log1p(-X.block(j,0,1,j) */

  W.setZero();
  W(0,0)=1;
  W.bottomLeftCorner(K-1,1) = Z.bottomLeftCorner(K-1,1);

  for (int j=1; j<K; j++) {
    W(j,j) = (1-Z.block(j,0,1,j).array().square()).sqrt().prod();
    for (int i=j+1; i<K; i++) {
      W(i,j) = Z(i,j)*(1-Z.block(i,0,1,j).array().square()).sqrt().prod();
    }
  }


  /* Rcout << "\nW:\n"; */
  /* for (int i=0; i<K; i++) { */
  /*   for (int j=0; j<K; j++) { */
  /*     Rcout << W(i,j) << "\t"; */
  /*   } */
  /*   Rcout << "\n"; */
  /* } */
  


  
}


#endif

