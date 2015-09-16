
#define EIGEN_MATRIX_PLUGIN <eigen_plugin.h>
#define EIGEN_SPARSEMATRIXBASE_PLUGIN <eigen_sparse_plugin.h>

#include <RcppEigen.h>
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <Eigen/Cholesky>
#include <cppad/cppad.hpp>
#include <except.h>
#include <LKJ_AD.h>


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

//' @title LKJ pdf
//' @param Y_ unconstrained vector
//' @param eta_ parameter >0
//' @param K dimension integer
//' @return log pdf
//' @export
//[[Rcpp::export]]
double LKJ(NumericVector Y_, double eta_, int K) {

  size_t v = Y_.size();
  VectorXA Y = VectorXd::Map(Y_.begin(),v).cast<AScalar>();
  AScalar eta = eta_;
  AScalar res = lkj_logpdf(Y, eta, K);
  return(Value(res));
}


//' @title Unwrap LKJ
//' @param Y input numeric vector
//' @param K integer dimension
//' @return matrix
//' @export
//[[Rcpp::export]]
List LKJ_unwrap(NumericVector Y_, int K) {


  VectorXA Y = MatrixXd::Map(Y_.begin(), K, K).cast<AScalar>();
  MatrixXA W(K,K);
  NumericMatrix M(K,K);
  AScalar log_jac = lkj_unwrap(Y, W);
  for (int j=0; j<K; j++) {
    for (int i=0; i<K; i++) {
      M(i,j) = Value(W(i,j));
    }
  }

  List res = List::create(
			  Named("chol_corr")=wrap(M),
			  Named("log_jac")=wrap(Value(log_jac))
			  );

  return(wrap(res));
}

  


