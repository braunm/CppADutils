#ifndef __mb_base
#define __mb_base

#define EIGEN_MATRIX_PLUGIN <eigen_plugin.h>
#define EIGEN_ARRAY_PLUGIN <eigen_plugin.h>
#define EIGEN_SPARSEMATRIX_PLUGIN <eigen_sparse_plugin.h>

#include <RcppEigen.h>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include "except.h"
#include <cppad/cppad.hpp>
#include <cppad/example/cppad_eigen.hpp>

using Rcpp::NumericVector;
using Rcpp::NumericMatrix;
using Rcpp::List;
using Eigen::Map;
using Eigen::VectorXd;
using Eigen::MatrixXd;


typedef CppAD::AD<double> AScalar;
typedef Eigen::Matrix<AScalar, Eigen::Dynamic, Eigen::Dynamic> MatrixXA;
typedef Eigen::Matrix<AScalar, Eigen::Dynamic, 1> VectorXA;
typedef Eigen::Triplet<double> TT;
typedef Eigen::SparseMatrix<double> SparseMatrixXd;
typedef std::vector< std::set<size_t> > sp_type;

template<typename TM>
class MB_Base {

 public :
  
  MB_Base(const List&);
  ~MB_Base();
  std::shared_ptr<TM> model;


  void hessian_init(const NumericVector&);
  
  double get_f(const NumericVector&);
  NumericVector get_df(const NumericVector&);
  List get_fdf(const NumericVector&);
  List get_fdfh(const NumericVector&);
  NumericMatrix get_hessian(const NumericVector&);
  Rcpp::S4 get_hessian_sparse(const NumericVector&);

  List get_hessian_test(const NumericVector&);

  
  List get_tape_stats();
  void record_tape(const NumericVector&);
  void record_to_abort(const NumericVector&, const size_t&);

  CppAD::ADFun<double> tape;
  bool tape_ready;

  inline size_t nvars() {return(nvars_);}

#ifdef MB_BASE_PLUGIN
#include MB_BASE_PLUGIN
#endif

 private:

  size_t nvars_;
  Eigen::VectorXi iRow; // row indices of nonzero elements
  Eigen::VectorXi jCol; // col indices of nonzero elements
  Eigen::VectorXi jpntr; // for each col, pointer to first element
  size_t nnz;
  CppAD::sparse_hessian_work hess_info;
  sp_type sp;
  bool hessian_initialized;
  std::vector<TT> trips;
  
};

template<typename TM>
MB_Base<TM>::MB_Base(const List& D) {
  model = std::make_shared<TM>(D);
  tape_ready = false;
  hessian_initialized = false;
}

template<typename TM>
MB_Base<TM>::~MB_Base() {
  model.reset();
}

template<typename TM>
void MB_Base<TM>::record_tape(const NumericVector& P_) {
  
  nvars_ = P_.size();
  VectorXA f(1); // to hold result 
  VectorXd Pd = VectorXd::Map(P_.begin(), nvars_);
  VectorXA P = Pd.cast<AScalar>();
  CppAD::Independent(P);
  f(0) = model->eval_f(P);
  tape.Dependent(P, f);
  //#ifdef NDEBUG
  tape.optimize();
  //#endif
  tape.check_for_nan(false);
  tape_ready = true;
}

template<typename TM>
void MB_Base<TM>::record_to_abort(const NumericVector& P_, const size_t& op) {
  
  nvars_ = P_.size();
  VectorXA f(1); // to hold result 
  VectorXd Pd = VectorXd::Map(P_.begin(), nvars_);
  VectorXA P = Pd.cast<AScalar>();
  CppAD::Independent(P, op);
  f(0) = model->eval_f(P);
  tape.Dependent(P, f);
  tape_ready = false;
}


template<typename TM>
void MB_Base<TM>::hessian_init(const NumericVector& P_)
{
			
// Initialization when sparsity pattern is unknown
  
  VectorXd P = VectorXd::Map(P_.begin(), nvars_); 
  hess_info.clear();
 

  sp_type identSet(nvars_);
  for (size_t ii=0; ii<nvars_; ii++) {
    identSet[ii].insert(ii);
  }
 
  // Rcout << "Computing sparsity pattern\n";
  tape.ForSparseJac(nvars_, identSet);
  sp_type s(1);
  s[0].insert(0);
  sp = tape.RevSparseHes(nvars_, s);

  // Rcout << "Computing nnz for one triangle
  size_t nnz_all = 0;

  for (auto ii=sp.begin(); ii != sp.end(); ++ii) {
    nnz_all += (*ii).size();
  }
  
  nnz = (nnz_all + nvars_)/2;
  
  // extracting row and column indices for lower triangle of Hessian
  iRow.resize(nnz);
  jCol.resize(nnz);
  size_t idx = 0;
  for (size_t ii=0; ii<sp.size(); ii++) {
    for (auto jj=sp[ii].begin(); jj != sp[ii].end(); ++jj) {
      if (*jj <= ii) {
    	iRow(idx) = ii;
    	jCol(idx) = *jj;
    	idx++;
      }
    }
  }

  VectorXd hess_vals(nnz);
  VectorXd w(1);
  w(0) = 1;
  tape.SparseHessian(P, w, sp, iRow, jCol, hess_vals, hess_info);

  hessian_initialized = true;

}

template<typename TM>
double MB_Base<TM>::get_f(const NumericVector& P_) {
  
  if (!tape_ready)
    throw MyException("tape not ready",__FILE__,__LINE__);

  VectorXd P = VectorXd::Map(P_.begin(), nvars_);
  VectorXd f(1);
  f = tape.Forward(0, P);
  
  //#ifndef NDEBUG
  size_t bad_tape = tape.compare_change_op_index();
  if (bad_tape != 0) {
    Rcout << "At index " << bad_tape << ".  Aborting now:\n";
    record_to_abort(P_, bad_tape);
    throw MyException("Tape operations incorrect.  Must retape\n", __FILE__, __LINE__);
  }
  //#endif
  
  return(f(0)); 
}

template<typename TM>
NumericVector MB_Base<TM>::get_df(const NumericVector& P_) {

  if (!tape_ready)
    throw MyException("tape not ready",__FILE__,__LINE__);

  VectorXd P = VectorXd::Map(P_.begin(), nvars_); 
  VectorXd df = tape.Jacobian(P);
  return(Rcpp::wrap(df));
}

template<typename TM>
List MB_Base<TM>::get_fdf(const NumericVector& P_) {

  if (!tape_ready)
    throw MyException("tape not ready",__FILE__,__LINE__);

  VectorXd P = VectorXd::Map(P_.begin(),nvars_); 
  VectorXd f(1);
  f = tape.Forward(0,P);

#ifndef NDEBUG
  size_t bad_tape = tape.compare_change_op_index();
  if (bad_tape != 0) {
    Rcout << "At index " << bad_tape << ".  Aborting now:\n";
    record_to_abort(P_, bad_tape);
    throw MyException("Tape operations incorrect.  Must retape\n", __FILE__, __LINE__);
  }
#endif
  
  VectorXd df = tape.Jacobian(P);

  //  assert(!df.hasNaN());
  List res = List::create(Rcpp::Named("val") = f,
			  Rcpp::Named("grad") = Rcpp::wrap(df)
			  );
  return(res);
}

template<typename TM>
List MB_Base<TM>::get_fdfh(const NumericVector& P_) {

  if (!tape_ready)
    throw MyException("tape not ready",__FILE__,__LINE__);

  // for DENSE Hessian

  VectorXd P = VectorXd::Map(P_.begin(),nvars_); 
  VectorXd f(1);
  f = tape.Forward(0,P);
  VectorXd df = tape.Jacobian(P);

  MatrixXd hess = MatrixXd::Map(tape.Hessian(P,size_t(0)).data(), nvars_, nvars_);
 
  List res = List::create(Rcpp::Named("val") = f,
			  Rcpp::Named("grad") = Rcpp::wrap(df),
			  Rcpp::Named("hessian") = Rcpp::wrap(hess)
			  );
  return(res);
}


template<typename TM>
NumericMatrix MB_Base<TM>::get_hessian(const NumericVector& P_) {

  // for DENSE Hessian
  if (!tape_ready)
    throw MyException("tape not ready",__FILE__,__LINE__);

  VectorXd P = VectorXd::Map(P_.begin(), nvars_);
  MatrixXd hess = MatrixXd::Map(tape.Hessian(P, size_t(0)).data(), nvars_, nvars_);
  return(Rcpp::wrap(hess));
}


template<typename TM>
List MB_Base<TM>::get_hessian_test(const NumericVector& P_) {

  // for DENSE Hessian
  if (!tape_ready)
    throw MyException("tape not ready",__FILE__,__LINE__);

  VectorXd P = VectorXd::Map(P_.begin(), nvars_);
 

  VectorXd w(1);
  w(0) = 1;
  
  VectorXd hessvec = tape.Hessian(P, w);
 

  //  MatrixXd hess = MatrixXd::Map(tape.Hessian(P, size_t(0)).data(), nvars_, nvars_);
  MatrixXd hess = MatrixXd::Map(hessvec.data(), nvars_, nvars_);


  List res = List::create(Rcpp::Named("pars") = Rcpp::wrap(P),
			  Rcpp::Named("hessvec") = Rcpp::wrap(hessvec),
			  Rcpp::Named("hess") = Rcpp::wrap(hess),
			  Rcpp::Named("nvars") = nvars_);

  return(res);
  
}




template<typename TM>
Rcpp::S4 MB_Base<TM>::get_hessian_sparse(const NumericVector& P_)
{
			
  VectorXd P = VectorXd::Map(P_.begin(), nvars_); 
 
  if (hessian_initialized) {
    // Rcout << "Computing hessian values\n"; 
    VectorXd hess_vals(nnz);
    VectorXd w(1);
    w(0) = 1;
    tape.SparseHessian(P,w,sp, iRow, jCol, hess_vals, hess_info);
    
    SparseMatrixXd H(nvars_, nvars_);
    std::vector<TT> trips;
    trips.resize(nnz);
    for (size_t k=0; k<nnz; k++) {
      TT tmp = TT(iRow(k), jCol(k), hess_vals(k));
      trips.push_back(tmp);
    }
    H.setFromTriplets(trips.begin(), trips.end());
    H.makeCompressed();
    
    SparseMatrixXd res(nvars_, nvars_);
    res = H.selfadjointView<Eigen::Lower>();
    return(Rcpp::wrap(res));

  } else {
    throw MyException("Hessian not initialized",__FILE__, __LINE__); 
  }
}

template<typename TM>
List MB_Base<TM>::get_tape_stats() {

  if (!tape_ready)
    throw MyException("tape not ready",__FILE__,__LINE__);

  List res = List::create(Rcpp::Named("nvars")=Rcpp::wrap(tape.Domain()),
			  Rcpp::Named("neq") = Rcpp::wrap(tape.Range()),
			  Rcpp::Named("size.var") = Rcpp::wrap(tape.size_var()),
			  Rcpp::Named("size.par") = Rcpp::wrap(tape.size_par()),
			  Rcpp::Named("size.op") = Rcpp::wrap(tape.size_op()),
			  Rcpp::Named("size.op.arg") = Rcpp::wrap(tape.size_op_arg()),
			  Rcpp::Named("size.op.seq") = Rcpp::wrap(tape.size_op_seq())
			  );
  
  return(res);
}

/*
template<typename TM>
void MB_Base<TM>::hessian_init(const NumericVector& P_,
			       const IntegerVector& hess_iRow,
			       const IntegerVector& hess_jCol)
{
  
  VectorXd P = VectorXd::Map(P_.begin(), nvars_); 
  hess_info.clear();

  nnz = hess_iRow.size();
  iRow.resize(nnz);
  jCol.resize(nnz);

  // Construct sparsity pattern
  sp.resize(nvars_);
  for (size_t k=0; k<nnz; k++) {
    iRow(k) = hess_iRow(k) - 1;
    jCol(k) = hess_jCol(k) - 1;
    sp[iRow[k]].insert(jCol[k]);
    sp[jCol[k]].insert(iRow[k]);
  }

  VectorXd w(1);
  w(0) = 1.0;
  VectorXd hess_vals(nnz);
  tape.SparseHessian(P, w, sp, iRow, jCol, hess_vals, hess_info);

  hessian_initialized = true;
 
}
*/


#endif



