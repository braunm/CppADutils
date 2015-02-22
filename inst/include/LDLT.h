#ifndef __CPPADUTILS_LDLT__
#define __CPPADUTILS_LDLT__

#define EIGEN_MATRIX_PLUGIN <eigen_plugin.h>
#define EIGEN_ARRAY_PLUGIN <eigen_plugin.h>
#define EIGEN_SPARSEMATRIX_PLUGIN <eigen_sparse_plugin.h>

#include <RcppEigen.h>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include "except.h"
#include <cppad/cppad.hpp>
#include <cppad/example/cppad_eigen.hpp>

typedef Matrix<AScalar, Dynamic, Dynamic> MatrixXA;
typedef Matrix<AScalar, Dynamic, 1> VectorXA;

void LDLT(const Ref<MatrixXA>& A,
	  const Ref<MatrixXA>& L, 
	  const Ref<VectorXA>& D) {
  // Assume A is lower for now

  size_t n = A.rows();
  assert(A.cols() == n);
  assert(L.rows() == n);
  assert(L.cols() == n);
  assert(D.size() == n); // D should be a vector

  L.setZero();
  D.setZero();

  VectorXA V;

  for (size_t j=1; j<=n; j++) {
    AScalar cjj = A(j-1, j-1);
    for (size_t s=1; s <= (j-1); s++) {
      cjj -= D(s-1) * L(j-1, s-1) * L(j-1, s-1);
    }
    D(j-1) = cjj;
    for (size_t i=j+1; i<=n; i++) {
      AScalar cij = A(i-1, j-1);
      for (size_t s=1; s<=j+1; s++) {
	cij -= D(s-1) * L(i-1, s-1) * L(j-1, s-1);
      }
      L(i-1, j-1) = cij / D(j-1);
    }
  }
}
 

}
	    
  





#endif
