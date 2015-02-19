
#define EIGEN_MATRIX_PLUGIN <eigen_plugin.h>
#define EIGEN_SPARSEMATRIXBASE_PLUGIN <eigen_sparse_plugin.h>

#include <RcppEigen.h>
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <Eigen/Cholesky>
#include <cppad/cppad.hpp>
#include <except.h>
#include <wish_AD.h>


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

typedef CppAD::AD<double> AScalar;
typedef Eigen::Matrix<AScalar, Dynamic, Dynamic> MatrixXA;
typedef Eigen::Matrix<AScalar, Dynamic, 1> VectorXA;
typedef Eigen::SparseMatrix<AScalar> SparseMatrixXA;

//[[Rcpp::export]]
double Wish_test(NumericMatrix X_, double nu_,
		       NumericMatrix S_){
	
  size_t k = X_.rows();

  MatrixXA X = MatrixXd::Map(X_.begin(), k, k).cast<AScalar>();
  AScalar nu = nu_;
  MatrixXA S = MatrixXd::Map(S_.begin(), k, k).cast<AScalar>();

  Eigen::LLT<MatrixXA> chol_X(X);
  MatrixXA chol_XL = chol_X.matrixL();
  Eigen::LLT<MatrixXA> chol_S(S);
  MatrixXA chol_SL = chol_S.matrixL();

  AScalar out = Wishart_logpdf(chol_XL, nu, chol_SL);
  double res = Value(out);
  
  return(res);
}


//[[Rcpp::export]]
double Inv_Wish_test(NumericMatrix X_, double nu_,
		       NumericMatrix S_){
	
  size_t k = X_.rows();

  MatrixXA X = MatrixXd::Map(X_.begin(), k, k).cast<AScalar>();
  AScalar nu = nu_;
  MatrixXA S = MatrixXd::Map(S_.begin(), k, k).cast<AScalar>();

  Eigen::LLT<MatrixXA> chol_X(X);
  MatrixXA chol_XL = chol_X.matrixL();
  Eigen::LLT<MatrixXA> chol_S(S);
  MatrixXA chol_SL = chol_S.matrixL();

  AScalar out = Inv_Wishart_logpdf(chol_XL, nu, chol_SL);
  double res = Value(out);
  
  return(res);
}

