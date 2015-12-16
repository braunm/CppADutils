
#define EIGEN_MATRIX_PLUGIN <eigen_plugin.h>
#define EIGEN_SPARSEMATRIXBASE_PLUGIN <eigen_sparse_plugin.h>

#include <RcppEigen.h>
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <Eigen/Cholesky>
#include <cppad/cppad.hpp>
#include <except.h>
#include <mat_normAD.h>


using Rcpp::NumericVector;
using Rcpp::NumericMatrix;
using Rcpp::as;
using Rcpp::List;
using Rcpp::Named;
using Rcpp::wrap;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::MatrixBase;
using Rcpp::Rcout;
using Eigen::Dynamic;
using Rcpp::S4;
using Eigen::SimplicialLLT;


typedef CppAD::AD<double> AScalar;
typedef Eigen::Matrix<AScalar, Dynamic, Dynamic> MatrixXA;
typedef Eigen::Matrix<AScalar, Dynamic, 1> VectorXA;
typedef Eigen::SparseMatrix<AScalar> SparseMatrixXA;

//' @title MatNorm_test
//' @param X_ numeric matrix (k x N)
//' @param M_ numeric matrix (k x N)
//' @param chol_U_ k x k lower triangular matrix of lower triangle of chol(U)
//' @param chol_V_ N x N lower triangular matrix of lower triangle of chol(V)
//' @param isPrec TRUE if U and V are precision matrices.  FALSE is U and V are covariance matrices
//' @return log pdf of matrix normal distribution
//[[Rcpp::export]]
double MatNorm_test(NumericMatrix X_, NumericMatrix M_,
		    NumericMatrix chol_U_,
		    NumericMatrix chol_V_,
		    bool isPrec){
  
  size_t k = chol_U_.rows();
  size_t N = chol_V_.rows();

  MatrixXA X = MatrixXd::Map(X_.begin(), k, N).cast<AScalar>();
  MatrixXA M = MatrixXd::Map(M_.begin(), k, N).cast<AScalar>();
  MatrixXA chol_U = MatrixXd::Map(chol_U_.begin(), k, k).cast<AScalar>();
  MatrixXA chol_V = MatrixXd::Map(chol_V_.begin(), N, N).cast<AScalar>();
 
  AScalar res = MatNorm_logpdf(X, M, chol_U, chol_V, isPrec);  
  return(Value(res));
}

//' @rdname MatNorm_test
//' @param X_ numeric matrix (k x N)
//' @param M_ numeric matrix (k x N)
//' @param U_ k x k dgCMatrix (full, not Cholesky)
//' @param V_ N x N dgCMatrix (full, not Cholesky)
//' @param isPrec TRUE if U and V are precision matrices.  FALSE is U and V are covariance matrices
//' @return log pdf of matrix normal distribution
//[[Rcpp::export]]
double MatNorm_sparse_test(NumericMatrix X_, NumericMatrix M_,
			   S4 U_, S4 V_, bool isPrec) {

  typedef Eigen::MappedSparseMatrix<double> MSparseXd;
  typedef Eigen::PermutationMatrix<Dynamic, Dynamic, int> Perm;
 
  
  size_t k = X_.rows();
  size_t N = X_.cols();

  MatrixXA X = MatrixXd::Map(X_.begin(), k, N).cast<AScalar>();
  MatrixXA M = MatrixXd::Map(M_.begin(), k, N).cast<AScalar>();

  SparseMatrixXA U = as<MSparseXd>(U_).cast<AScalar>();
  SimplicialLLT<SparseMatrixXA> ChU(U);
  Perm PU = ChU.permutationP();
  SparseMatrixXA LU = ChU.matrixL();


  SparseMatrixXA V = as<MSparseXd>(V_).cast<AScalar>();
  SimplicialLLT<SparseMatrixXA> ChV(V);
  Perm PV = ChV.permutationP();
  SparseMatrixXA LV = ChV.matrixL();

  AScalar res = MatNorm_logpdf(X, M, LU, PU, LV, PV, isPrec);
  return(Value(res));
}

 
