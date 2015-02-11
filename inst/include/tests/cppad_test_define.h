#ifndef __CPPADUTILS_TEST_DEFINE
#define __CPPADUTILS_TEST_DEFINE

#define EIGEN_MATRIX_PLUGIN <eigen_plugin.h>
#define EIGEN_SPARSEMATRIX_PLUGIN <eigen_sparse_plugin.h>

#include <RcppEigen.h>
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <cppad/cppad.hpp>
#include <cppad/example/cppad_eigen.hpp>


using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::VectorXi;
using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::SparseMatrix;
using Eigen::Lower;
using Rcpp::List;
using Rcpp::NumericVector;
using Rcpp::Named;
using Rcpp::wrap;
using Rcpp::Rcout;

typedef CppAD::AD<double> AScalar;
typedef Matrix<AScalar, Dynamic, Dynamic> MatrixXA;
typedef Matrix<AScalar, Dynamic, 1> VectorXA;
typedef Eigen::Triplet <double> TT;
typedef SparseMatrix<double> SparseMatrixXd;
typedef std::vector< std::set<size_t> > sp_set;


// functor class that points to another function 

template<typename TM>
class run_test {

  CppAD::ADFun<double> tape;
  
  CppAD::sparse_hessian_work hess_info;

  int nvars;
  VectorXA f; // to hold result
  VectorXA P;

  std::vector<TT> tripletList;
  VectorXd X;

  std::shared_ptr<TM> func;


 public:
  
  run_test() {
    // constructor
    f.resize(1);
    func = std::make_shared<TM>();
  }

  void record_tape(const NumericVector& P_) {
    
    nvars = P_.size();
    X.resize(nvars);
    X = VectorXd::Map(P_.begin(), nvars);
    VectorXA P = X.cast<AScalar>();
    
    CppAD::Independent(P);
    f(0) = func->eval(P);
    tape.Dependent(P, f);
    tape.optimize();
  }
  
  Rcpp::List cppad_results() {


    VectorXd w(1);
    w(0) = 1.0;

    
 
      
    VectorXd val(1);
    val = tape.Forward(0, X);
    VectorXd grad = tape.Jacobian(X);
    MatrixXd hess_dense = MatrixXd::Map(tape.Hessian(X, size_t(0)).data(), nvars, nvars);
    
    /* Rcout << "grad:\n" << grad << "\n\n"; */
    /* Rcout << "hess_dense:\n" << hess_dense << "\n\n"; */

    MatrixXd hess_sparse = MatrixXd::Map(tape.SparseHessian(X, w).data(),
					 nvars, nvars);
   

    
    //    Rcout << "hess_sparse from driver:\n" << hess_sparse << "\n\n";
    
    // make identity matrices
    
    sp_set identSet(nvars); 
    for (size_t ii=0; ii<nvars; ii++) {
      identSet[ii].insert(ii);
    }
    
    // forward Jacobian and Reverse Hessian sparsity
    val = tape.Forward(0, X);
    grad = tape.Jacobian(X);
    sp_set sp = tape.ForSparseJac(nvars, identSet);
    sp_set zset(1);
    zset[0].insert(0);
    sp_set hs = tape.RevSparseHes(nvars, zset);  
    

    int nnz_all = 0;
    for (int ii=0; ii<nvars; ii++) {
      nnz_all += hs[ii].size();
    }
    int nnz_LT = (nnz_all + nvars)/2;
    
    VectorXi iRow(nnz_LT);
    VectorXi jCol(nnz_LT);
    size_t idx = 0;
    for (size_t ii=0; ii<nvars; ii++) {      
      for (auto jj = hs[ii].begin(); jj != hs[ii].end(); ++jj) {
      	if (*jj <= ii) {
      	  iRow(idx) = ii;
      	  jCol(idx) = *jj;
      	  idx++;
      	}
      }      
    }

    VectorXd hess_vals(nnz_LT);
    hess_info.clear();
    size_t n_sweep = tape.SparseHessian(X, w, hs, iRow, jCol, hess_vals, hess_info);
    
    std::vector<TT> trips;
    trips.resize(nnz_LT);
    for (int k=0; k<hess_vals.size(); k++) {
      TT tmp = TT(iRow(k), jCol(k), hess_vals(k));
      trips.push_back(tmp);
    }
    
    SparseMatrixXd H_LT;
    H_LT.resize(nvars, nvars);
    H_LT.setFromTriplets(trips.begin(), trips.end());
    H_LT.makeCompressed();
    
    /* CppAD::sparse_hessian_work info_new; */
    /* std::vector<TT> tripsNew; */
    /* tripsNew = tape.SparseHessian(X, w, hs, info_new); */
    /* SparseMatrixXd H_new; */
    /* H_new.resize(nvars, nvars); */
    /* H_new.setFromTriplets(tripsNew.begin(), tripsNew.end()); */
    /* H_new.makeCompressed(); */
    
    /* Rcpp::NumericVector err(1); */
    /* int vb = INT_MAX; */
    /* err[vb] = 0; */


    List res = List::create(Named("val") = wrap(val),
			  Named("grad") = wrap(grad),
			  Named("hess.dense") = wrap(hess_dense),
			  Named("hess.sp") = wrap(hess_sparse),
			  Named("hess.spLT") = wrap(H_LT),
			  //			  Named("new.driver") = wrap(H_new)
			  Named("n_sweep") = wrap(n_sweep)			  
			  );
  return(res);
  }
};

#endif
