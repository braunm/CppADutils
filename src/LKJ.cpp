
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


//[[Rcpp::export]]
double LKJ(NumericVector Y_, double eta_, int K) {

  size_t v = Y_.size();
  VectorXA Y = VectorXd::Map(Y_.begin(),v).cast<AScalar>();
  AScalar eta = eta_;
  AScalar res = lkj_logpdf(Y, eta, K);
  return(Value(res));
}

//[[Rcpp::export]]
NumericMatrix LKJ_unwrap(NumericVector Y_, int K) {


  VectorXA Y = MatrixXd::Map(Y_.begin(), K, K).cast<AScalar>();
  MatrixXA W(K,K);
  NumericMatrix res(K,K);
  lkj_unwrap(Y, W);
  for (int j=0; j<K; j++) {
    for (int i=0; i<K; i++) {
      res(i,j) = Value(W(i,j));
    }
  }

  return(wrap(res));
}

  


