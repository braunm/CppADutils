#ifndef __CPPADUTILS_LDLT__
#define __CPPADUTILS_LDLT__

#define EIGEN_MATRIX_PLUGIN <eigen_plugin.h>
#define EIGEN_ARRAY_PLUGIN <eigen_plugin.h>
#define EIGEN_SPARSEMATRIX_PLUGIN <eigen_sparse_plugin.h>

#include <RcppEigen.h>
#include <Eigen/Eigen>
#include "except.h"
#include <cppad/cppad.hpp>
#include <cppad/example/cppad_eigen.hpp>

using Eigen::Dynamic;
using Eigen::MatrixBase;

typedef CppAD::AD<double> AScalar;
typedef Eigen::Matrix<AScalar, Dynamic, Dynamic> MatrixXA;
typedef Eigen::Matrix<AScalar, Dynamic, 1> VectorXA;

template<typename TA, typename TL, typename TD>
void LDLT(const MatrixBase<TA>& A,
	  const MatrixBase<TL>& L_, 
	  const MatrixBase<TD>& D_) {
  // Assume A is lower for now

  MatrixBase<TL>& L = const_cast<MatrixBase<TL>&>(L_);
  MatrixBase<TD>& D = const_cast<MatrixBase<TD>&>(D_);
  
  size_t n = A.rows();
  assert(A.cols() == n);
  assert(L.rows() == n);
  assert(L.cols() == n);
  assert(D.size() == n); // D should be a vector

  L.setIdentity(n,n);
  D.setZero();
  
  for (size_t j=1; j<=n; j++) {
     D(j-1) = A(j-1,j-1);
    for (size_t s=1; s <= (j-1); s++) {
      D(j-1) -= D(s-1) * L(j-1, s-1) * L(j-1, s-1);
    }
    for (size_t i=j+1; i<=n; i++) {
      L(i-1, j-1) = A(i-1, j-1);
      for (size_t s=1; s<=j-1; s++) {
	L(i-1,j-1) -= D(s-1) * L(i-1, s-1) * L(j-1, s-1);
      }
      L(i-1, j-1) = L(i-1, j-1) / D(j-1);
    }
  }
}




#endif
