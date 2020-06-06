
#define EIGEN_MATRIX_PLUGIN <eigen_plugin.h>
#define EIGEN_SPARSEMATRIXBASE_PLUGIN <eigen_sparse_plugin.h>

 #define CPPAD_USE_CPLUSPLUS_2011 1

#include <RcppEigen.h>
#include <Eigen/Eigen>
#include <cppad/cppad.hpp>
#include <except.h>
#include <distributions.h>


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


//' @title CP_erf test
//' @param X_ vector
//' @return Numeric vector
//' @export
//[[Rcpp::export]]
NumericVector CP_erf(NumericVector X_) {
       
  size_t k = X_.size();

  VectorXA X = VectorXd::Map(X_.begin(), k).cast<AScalar>();

  VectorXA out(k);

  // MB_erf(X, out);

  NumericVector res(k);
  for (size_t i=0; i<k; i++)
    //    res(i) = Value(out(i));
    res(i)=Value(CppAD::erf(X(i)));
  return(res);
}


// TO DO:  Not convinced the CppAD erf is accurate
//  Might need to create a complete MB_erf atomic

//' @title MB_erf test
//' @param X_ vector
//' @return Numeric vector
//' @export
//[[Rcpp::export]]
NumericVector MB_erf(NumericVector X_) {
       
  size_t k = X_.size();

  VectorXA X = VectorXd::Map(X_.begin(), k).cast<AScalar>();

  VectorXA out(k);

  // MB_erf(X, out);

  NumericVector res(k);
  for (size_t i=0; i<k; i++)
    //    res(i) = Value(out(i));
    res(i)=Value(MB_erf(X(i)));
  return(res);
}




