#ifndef __MATNORM__AD
#define __MATNORM__AD


#include <iostream>
using Eigen::MatrixBase;
using Eigen::Lower;
using Eigen::Dynamic;
using Eigen::Matrix;
using Eigen::OnTheRight;
using Eigen::SparseMatrixBase;
using Eigen::PermutationBase;

template<typename TY, typename TM, typename TG, typename TA, typename TP>
typename TM::Scalar MatNorm_logpdf(const MatrixBase<TY>& Y,
				   const MatrixBase<TM>& M,
				   const SparseMatrixBase<TG>& LU, // k x k matrix
				   const PermutationBase<TP>& PU,
				   const SparseMatrixBase<TA>& LV, // N x N matrix
				   const PermutationBase<TP>& PV,
				   const bool isPrec
				   ) {
  
  typedef typename TM::Scalar Scalar; // TY, TM and TG must all have same scalar type
  typedef Matrix<Scalar, Dynamic, Dynamic> MatrixXS;  

  /* chol_U is Cholesky of row cov/prec
     chol_col_G is Cholesky of col cov/prec
  */

  const int k = Y.rows();
  const int N = Y.cols();

  assert(LV.rows() == N); // row cov has number of columns in Y
  assert(LU.cols() == k); // col cov has number of rows in Y 

  Scalar c = k * N * M_LN_SQRT_2PI;

  Scalar detU = 0.;
  for (int i=0; i<k; i++) {
    detU += log(LU.template derived().coeff(i,i));
  }

  Scalar detV = 0.;
  for (int i=0; i<N; i++) {
    detV += log(LV.template derived().coeff(i,i));
  }
  
  
  MatrixXS Z = PU*(Y-M)*PV.inverse();
  Scalar res;
 
  if (isPrec) {
 
    Z = (LU.template triangularView<Lower>().transpose() * Z).eval() * LV.template triangularView<Lower>();
    res = -0.5*Z.squaredNorm() - c + k*detV + N*detU; 
   
  } else {
    MatrixXS T = LU.template triangularView<Lower>().solve(Z);
    Z.transpose() = LV.template triangularView<Lower>().solve(T.transpose());
    res = -0.5*Z.squaredNorm() - c - k*detV - N*detU; 
  }
  return(res);

}





template<typename TY, typename TM, typename TG, typename TA>
typename TM::Scalar MatNorm_logpdf(const MatrixBase<TY>& Y,
				   const MatrixBase<TM>& M,
				   const MatrixBase<TG>& chol_U, // k x k matrix
				   const MatrixBase<TA>& chol_V, // N x N matrix
				   const bool isPrec
				   ) {
  
  typedef typename TM::Scalar Scalar; // TY, TM and TG must all have same scalar type
  typedef Matrix<Scalar, Dynamic, Dynamic> MatrixXS;  

  /* chol_U is Cholesky of row cov/prec
   chol_col_G is Cholesky of col cov/prec
  */

  const int k = Y.rows();
  const int N = Y.cols();

  assert(chol_V.rows() == N); // row cov has number of columns in Y
  assert(chol_U.cols() == k); // col cov has number of rows in Y 

  Scalar c = k * N * M_LN_SQRT_2PI;
  Scalar detV = chol_V.diagonal().array().log().sum();
  Scalar detU = chol_U.diagonal().array().log().sum();

  MatrixXS Z(k,N);
  Scalar res;
 
  if (isPrec) {
 
    Z = (chol_U.template triangularView<Lower>().transpose() * (Y-M)).eval() * chol_V.template triangularView<Lower>();
    res = -0.5*Z.squaredNorm() - c + k*detV + N*detU; 
   
  } else {
    Z = chol_U.template triangularView<Lower>().solve(Y-M);
    Z.transpose() = chol_V.template triangularView<Lower>().solve(Z.transpose());
    res = -0.5*Z.squaredNorm() - c - k*detV - N*detU; 
  }
  return(res);

}


#endif

