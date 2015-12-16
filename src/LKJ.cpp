
#define EIGEN_MATRIX_PLUGIN <eigen_plugin.h>
#define EIGEN_SPARSEMATRIXBASE_PLUGIN <eigen_sparse_plugin.h>

#include <RcppEigen.h>
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <Eigen/Cholesky>
#include <cppad/cppad.hpp>
#include <except.h>
#include <LKJ_AD.h>
#include <stan/math/prim/mat/prob/lkj_corr_log.hpp>
#include <stan/math/prim/mat/prob/lkj_corr_cholesky_log.hpp>



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




//' @title LKJ chol stan
//' @param L Cholesky
//' @param eta parameter
//' @return lkj pdf
//' @export
//[[Rcpp::export]]
double LKJ_chol_stan(NumericMatrix L_, double eta) {

  int K = L_.rows();
  MatrixXd L = MatrixXd::Map(L_.begin(), K, K);
  double res = stan::math::lkj_corr_cholesky_log(L,eta);
  return(res);
}


//' @title LKJ stan
//' @param L correlation matrix
//' @param eta parameter
//' @return lkj pdf
//' @export
//[[Rcpp::export]]
double LKJ_stan(NumericMatrix S_, double eta) {

  int K = S_.rows();
  MatrixXd S = MatrixXd::Map(S_.begin(), K, K);
  double res = stan::math::lkj_corr_log(S,eta);
  return(res);
}

//' @title LKJ stan constant
//' @param eta parameter
//' @param K integer
//' @return lkj constant
//' @export
//[[Rcpp::export]]
double LKJ_stan_const(double eta, int K) {
  double res = stan::math::do_lkj_constant(eta, K);
  return(res);
}

//' @export
//[[Rcpp::export]]
List corr_matrix(NumericVector X_,
		 int K) {
  
  int Kch2 = K*(K-1)/2;
  VectorXd X = VectorXd::Map(X_.begin(),Kch2);
  MatrixXd Z(K,K);
  double lp=0;
  Z = stan::math::corr_matrix_constrain(X, K, lp);

  List res = List::create(
			  Named("S")=wrap(Z),
			  Named("lp")=wrap(lp)
			  );
  
  return(wrap(res));
}

//' @export
//[[Rcpp::export]]
List chol_corr_matrix(NumericVector X_,
		 int K) {
  
  int Kch2 = K*(K-1)/2;
  VectorXd X = VectorXd::Map(X_.begin(),Kch2);
  MatrixXd Z(K,K);
  double lp=0;
  Z = stan::math::cholesky_corr_constrain(X, K, lp);

  List res = List::create(
			  Named("S")=wrap(Z),
			  Named("lp")=wrap(lp)
			  );
  
  return(wrap(res));

  
}

