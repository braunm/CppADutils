
#define EIGEN_MATRIX_PLUGIN <eigen_plugin.h>
#define EIGEN_SPARSEMATRIXBASE_PLUGIN <eigen_sparse_plugin.h>

// #define CPPAD_USE_CPLUSPLUS_2011 1

#include <RcppEigen.h>
#include <Eigen/Eigen>
#include <cppad/cppad.hpp>
#include <except.h>
#include <gamma_AD.h>


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


//' @title gamma_logpdf test
//' @param X_ vector
//' @param R_ vector
//' @param A_ vector
//' @return Numeric vector
//' @export
//[[Rcpp::export]]
NumericVector dgamma_test(NumericVector X_, NumericVector R_,
		       NumericVector A_){
	

  size_t k = X_.size();

  VectorXA X = VectorXd::Map(X_.begin(), k).cast<AScalar>();
  VectorXA R = VectorXd::Map(R_.begin(), k).cast<AScalar>();
  VectorXA A = VectorXd::Map(A_.begin(), k).cast<AScalar>();

  VectorXA out(k);

  gamma_logpdf(X, R, A, out);

  NumericVector res(k);
  for (size_t i=0; i<k; i++)
    res(i) = Value(out(i));
  
  return(res);
}




