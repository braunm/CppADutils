
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


//' @title LKJ pdf
//' @param Y_ unconstrained vector
//' @param eta_ parameter >0
//' @param K dimension integer
//' @return log pdf
//' @export
//[[Rcpp::export]]
double LKJ_chol(NumericMatrix L_, double eta_) {

  size_t K = L_.rows();
  MatrixXA L = MatrixXd::Map(L_.begin(),K,K).cast<AScalar>();
  AScalar eta = eta_;
  AScalar res = lkj_chol_logpdf(L, eta);
  return(Value(res));
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




//' @title Unwrap LKJ
//' @param Y input numeric vector
//' @param K integer dimension
//' @return matrix
//' @export
//[[Rcpp::export]]
List LKJ_unwrap(NumericVector Y_, int K) {

  const int Kch2 = R::choose(K,2);
  VectorXA Y = VectorXd::Map(Y_.begin(), Kch2).cast<AScalar>();
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

//' @ LKJ chol stan
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


//' @ LKJ stan
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

//' @ LKJ stan constant
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

