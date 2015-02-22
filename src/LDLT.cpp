
#define EIGEN_MATRIX_PLUGIN <eigen_plugin.h>
#define EIGEN_SPARSEMATRIXBASE_PLUGIN <eigen_sparse_plugin.h>

#include <RcppEigen.h>
#include <Eigen/Eigen>
#include <except.h>
#include <LDLT_cppad.h>



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
List LDLT_test(const NumericMatrix& A_) {

  MatrixXA A = as<MatrixXd>(A_).cast<AScalar>();
  int n = A.rows();
  MatrixXA L_(n,n);
  VectorXA D_(n);

  LDLT(A, L_, D_);

  NumericMatrix L(n,n);
  NumericVector D(n);

  std::fill(L.begin(), L.end(), 0.0);

  for (int j=0; j<n; j++) {
    D(j) = Value(D_(j));
    for (int i=j; i<n; i++) {
      L(i,j) = Value(L_(i,j));
    }
  }

  List res = List::create(Named("L")=wrap(L),
			  Named("D")=wrap(D));

  return(res);
}
