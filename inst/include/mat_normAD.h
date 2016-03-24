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



template<typename TM>
  typename TM::Scalar MatNorm_log_const(const MatrixBase<TM>& chol_U,
					const MatrixBase<TM>& chol_V,
					const bool isPrec) {

  typedef typename TM::Scalar Scalar;

  // log constant of matrix normal.  U and V must be of same type

  const int k = chol_U.rows();
  const int N = chol_V.rows();

  Scalar LD_U = chol_U.diagonal().array().log().sum();
  Scalar LD_V = chol_V.diagonal().array().log().sum();

  Scalar res;
  if (isPrec) {
    res = -N*k*M_LN_SQRT_2PI + N*LD_U + k*LD_V;
  } else {
    res = -N*k*M_LN_SQRT_2PI - N*LD_U - k*LD_V;
  }
  
  return(res);
}

template<typename TM>
  typename TM::Scalar MatNorm_log_const(const SparseMatrixBase<TM>& chol_U,
					const SparseMatrixBase<TM>& chol_V,
					const bool isPrec) {

  typedef typename TM::Scalar Scalar;


  // log constant of matrix normal.  U and V must be of same type

  const int k = chol_U.rows();
  const int N = chol_V.rows();

  Scalar LD_U = 0.;
  for (int i=0; i<k; i++) {
    LD_U += log(chol_U.derived().coeff(i,i));
  }

  Scalar LD_V = 0.;
  for (int i=0; i<N; i++) {
    //    LD_V += log(chol_V.template derived().coeff(i,i));
    LD_V += log(chol_V.derived().coeff(i,i));
  }
  

  Scalar res;
  if (isPrec) {
    res = -N*k*M_LN_SQRT_2PI + N*LD_U + k*LD_V;
  } else {
    res = -N*k*M_LN_SQRT_2PI - N*LD_U - k*LD_V;
  }
  
  return(res);
}



template<typename TY, typename TM, typename TG, typename TA>
typename TM::Scalar MatNorm_logpdf(const MatrixBase<TY>& Y, // k x N matrix
				   const MatrixBase<TM>& M, // k x N matrix
				   const MatrixBase<TG>& chol_U, // k x k matrix
				   const MatrixBase<TA>& chol_V, // N x N matrix
				   const bool isPrec
				   ) {
  
  typedef typename TM::Scalar Scalar; // TY, TM and TG must all have same scalar type
  typedef Matrix<Scalar, Dynamic, Dynamic> MatrixXS;  

  /* 
     Y is a k x N matrix
     chol_U is Cholesky of row cov/prec
     chol_col_G is Cholesky of col cov/prec
  */

  const int k = Y.rows();
  const int N = Y.cols();

  //  assert(chol_V.rows() == N); // row cov has number of columns in Y
  // assert(chol_U.cols() == k); // col cov has number of rows in Y 

  Scalar c = MatNorm_log_const(chol_U, chol_V, isPrec);

  MatrixXS A(k,N);
  MatrixXS Z(N,k);

  MatrixXS YM = Y-M;
 
  if (isPrec) {
    A = chol_V.template triangularView<Lower>().transpose() * YM.transpose();
    Z = A * chol_U.template triangularView<Lower>();
  } else {
    A = chol_U.template triangularView<Lower>().solve(YM);
    Z = chol_V.template triangularView<Lower>().solve(A.transpose());
  }

  Scalar res = c - 0.5*Z.squaredNorm();

  return(res);

}




template<typename TY, typename TG, typename TA, typename TP>
typename TY::Scalar MatNorm_logpdf(const MatrixBase<TY>& Y,
				   const MatrixBase<TY>& M,
				   const SparseMatrixBase<TG>& LU, // k x k matrix
				   const PermutationBase<TP>& PU,
				   const SparseMatrixBase<TA>& LV, // N x N matrix
				   const PermutationBase<TP>& PV,
				   const bool isPrec
				   ) {
  
  typedef typename TY::Scalar Scalar; // TY  and TG must have same scalar type
  typedef Matrix<Scalar, Dynamic, Dynamic> MatrixXS;  

  // chol_U is Cholesky of row cov/prec
    // chol_col_G is Cholesky of col cov/prec

  Scalar c = MatNorm_log_const(LU, LV, isPrec);
  
  MatrixXS Z = PU*(Y-M)*PV.inverse();
 
  if (isPrec) {

    MatrixXS T = LU.template triangularView<Lower>().transpose() * Z;
    Z = T * LV.template triangularView<Lower>();
   
  } else {
    MatrixXS T = LU.template triangularView<Lower>().solve(Z);
    Z.transpose() = LV.template triangularView<Lower>().solve(T.transpose()); 
  }

  Scalar res = c - 0.5 * Z.squaredNorm();
  
  return(res);

}







#endif

