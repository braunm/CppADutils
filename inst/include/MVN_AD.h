#ifndef __AD_MVNORM
#define __AD_MVNORM

using Eigen::MatrixBase;
using Eigen::Lower;
using Eigen::Dynamic;
using Eigen::Matrix;
using Eigen::Upper;
using Eigen::Lower;
using Eigen::SparseMatrixBase;
using Eigen::PermutationMatrix;
using Eigen::PermutationBase;
using Eigen::SimplicialLLT;

typedef CppAD::AD<double> AScalar;
typedef Matrix<AScalar, Dynamic, Dynamic> MatrixXA;
typedef Matrix<AScalar, Dynamic, 1> VectorXA;
typedef Matrix<AScalar, 1, Dynamic> RowVectorXA;


template<typename TY, typename TM, typename TG, typename TQ>
void MVN_logpdf(const MatrixBase<TY>& Y,
		const MatrixBase<TM>& mu,
		const MatrixBase<TG>& chol_G,
		const MatrixBase<TQ>& out_,
		const bool isPrec = true) {
  
 
  
  MatrixBase<TQ>& out = const_cast<MatrixBase<TQ>& >(out_);

  /* 
     If isPrec is true, chol_G is Cholesky of the precision matrix.
     If false, it is the Cholesky of the covariance matrix. 
     Each draw is a column.
  */

  const int k = Y.rows();
  const int N = Y.cols();

  MatrixXA Z = Y;
  AScalar normConst;

  if (mu.rows() == k) {
    if (mu.cols() == N) {
      Z -= mu; 
    } else {
      if (mu.cols() == 1) {
	Z.colwise() -= mu.col(0); 
      } else {
	throw MyException("MVN: dimension mismatch",
			  __FILE__, __LINE__);
      }
    }
  }

  AScalar cholGlogDet = chol_G.diagonal().array().log().sum();

  if (isPrec) {
    normConst = cholGlogDet - k * M_LN_SQRT_2PI;
    Z = chol_G.transpose().template triangularView<Upper>() * Z;
  } else {
    normConst = - cholGlogDet - k * M_LN_SQRT_2PI;
    chol_G.template triangularView<Lower>().solveInPlace(Z);
  }

  VectorXA tmp = normConst - (0.5*Z.array()*Z.array()).colwise().sum();
  out = tmp;


    //out.array() = normConst - 0.5 * (Z.array()*Z.array()).colwise().sum();

}

template<typename TY, typename TM, typename TG, typename TQ, typename TP>
void MVN_logpdf(const MatrixBase<TY>& Y,
		const MatrixBase<TM>& mu,
		const SparseMatrixBase<TG>& L, // lower triangle of Cholesky
		const PermutationBase<TP>& P, // permutation matrix for Cholesky
		const MatrixBase<TQ>& out_,
		const bool isPrec = true) {
  

  
  MatrixBase<TQ>& out = const_cast<MatrixBase<TQ>& >(out_);

  /* if isPrec is true, chol_G is Cholesky of the precision matrix.
     If false, it is the Cholesky of the covariance matrix
  */

  const int k = Y.rows();
  const int N = Y.cols();

  MatrixXA Z = Y;
  AScalar normConst;

  if (mu.rows() == k) {
    if (mu.cols() == N) {
      Z -= mu; // Y and mu match
    } else {
      if (mu.cols() == 1) {
	Z.colwise() -= mu.col(0); //common mu for each column 
      } else {
	throw MyException("MVN: dimension mismatch",__FILE__, __LINE__);
      }
    }
  } 

  AScalar Ldet = 0.;
  for (int i=0; i<k; i++) {
    //    Ldet += log(L.template derived().coeff(i,i));
    Ldet += log(L.derived().coeff(i,i));
  }

  Z = P*Z;
  MatrixXA tmp(k,N);

  if (isPrec) {
    normConst = Ldet - k * M_LN_SQRT_2PI;
    Z = L.transpose().template triangularView<Upper>()  * Z;
  } else {
    normConst = -Ldet - k * M_LN_SQRT_2PI;
    Z = L.template triangularView<Lower>().solve(Z);
  }
  
  out.array() = normConst - 0.5 * (Z.array() * Z.array()).colwise().sum();


}






#endif

