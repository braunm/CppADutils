
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


//' @title LKJ test
//' @param Y_ unconstrained vector
//' @param eta_ parameter >0
//' @param K dimension integer
//' @return List
//' @export
//[[Rcpp::export]]
List LKJ_test(NumericVector Y_, double eta_, int K) {

  size_t v = Y_.size();
  VectorXA Y = VectorXd::Map(Y_.begin(),v).cast<AScalar>();
  AScalar eta = eta_;
  MatrixXA L = MatrixXA::Zero(K,K);
  AScalar logjac = lkj_unwrap(Y, L);
  AScalar lkj_chol = lkj_chol_logpdf(L, eta);
  AScalar logpdf = lkj_chol + logjac;
  NumericMatrix M(K,K);
  for (int j=0; j<K; j++) {
    for (int i=0; i<K; i++) {
      if (i<j) {
	M(i,j) = 0;
      } else {
	M(i,j) = Value(L(i,j));
      }
    }
  }

  List res = List::create(
			  Named("L")=wrap(M),
			  Named("logpdf")=wrap(Value(logpdf)),
			  Named("logjac")=wrap(Value(logjac))
			  );
  
  return(wrap(res));
}


//' @title LKJ constant term in pdf
//' @param eta_ parameter >0
//' @param K dimension integer
//' @return log pdf
//' @export
//[[Rcpp::export]]
double LKJ_const(double eta, int K) {
  AScalar c = lkj_const(eta, K);
  return(Value(c));
}


