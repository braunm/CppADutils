#ifndef __mb_base
#define __mb_base

#define EIGEN_MATRIX_PLUGIN <cppad/example/eigen_plugin.hpp>
#define EIGEN_ARRAY_PLUGIN <cppad/example/eigen_plugin.hpp>
#define EIGEN_SPARSEMATRIXBASE_PLUGIN <cppad/example/eigen_plugin.hpp>

#include <RcppEigen.h>
#include <cppad/cppad.hpp>
#include <cppad/example/cppad_eigen.hpp>


using Rcpp::NumericVector;
using Rcpp::NumericMatrix;
using Rcpp::IntegerVector;
using Rcpp::List;
using Eigen::Map;
using Eigen::Matrix;
using Eigen::Dynamic;
using Rcpp::as;
using Rcpp::wrap;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::VectorXi;
using Eigen::MatrixXi;
using Eigen::SparseMatrix;

template<typename TM>
RCPP_EXPOSED_CLASS(MB_Base<TM>)
class MB_Base {
 public :
  
  typedef int Index;
  typedef double Scalar;
  typedef CppAD::AD<Scalar> AScalar;
  
  typedef Matrix<Scalar, Dynamic, Dynamic> MatrixXS;
  typedef Matrix<Scalar, Dynamic, 1> VectorXS;

  typedef Matrix<AScalar, Dynamic, Dynamic> MatrixXA;
  typedef Matrix<AScalar, Dynamic, 1> VectorXA;

  int nvars;

#ifdef COMPILE_SPARSE_HESSIAN

  typedef std::vector<std::set<size_t> > sp_type;

  VectorXi iRow; // row indices of nonzero elements
  VectorXi jCol; // col indices of nonzero elements
  VectorXi jpntr; // for each col, pointer to first element

  int nnz;

  CppAD::sparse_hessian_work hess_info;
  sp_type sp;

  void hessian_init_nopattern(const NumericVector&);
  void hessian_init(const NumericVector&,
		    const IntegerVector&,
		    const IntegerVector&);

  bool hessian_initialized;

#endif

  MB_Base(const List&);

  std::shared_ptr<TM> model;
 
  NumericVector get_f(const NumericVector&);
  NumericVector get_f_direct(const NumericVector&);
  double get_LL(const NumericVector&);
  NumericVector get_df(const NumericVector&);
  List get_fdf(const NumericVector&);
  List get_fdfh(const NumericVector&);
  NumericMatrix get_hessian(const NumericVector&);
  Rcpp::S4 get_hessian_sparse(const NumericVector&);

  List get_tape_stats();
  void record_tape(const NumericVector&);

  CppAD::ADFun<Scalar> tape;
  bool tape_ready;


#ifdef MB_BASE_PLUGIN
#include MB_BASE_PLUGIN
#endif

};

template<typename TM>
MB_Base<TM>::MB_Base(const List& D) {

  model = std::make_shared<TM>();
  tape_ready = false;
#ifdef COMPILE_SPARSE_HESSIAN
  hessian_initialized = false;
#endif


}


template<typename TM>								       
void MB_Base<TM>::record_tape(const NumericVector& P_) {
  
  nvars = P_.size();
  VectorXA f(1); // to hold result 
  VectorXS Pd = VectorXS::Map(P_.begin(),nvars);
  VectorXA P = Pd.cast<AScalar>();
  std::cout << "Tape recording.  nvars = " << nvars << "\n";
  CppAD::Independent(P);
  f(0) = model->eval_f(P);
  tape.Dependent(P, f);
  Rcpp::Rcout << "Finished recording\n";
#ifdef NDEBUG
  Rcpp::Rcout << "Optimizing tape\n";
  tape.optimize();
  Rcpp::Rcout << "Tape optimized\n";
#endif
  tape_ready = true;
}

template<typename TM>
NumericVector MB_Base<TM>::get_f(const NumericVector& P_) {
  
  if (!tape_ready)
    throw MyException("tape not ready",__FILE__,__LINE__);

  VectorXS P = VectorXS::Map(P_.begin(),nvars);
  VectorXS f(1);
  f = tape.Forward(0, P);
  
#ifndef NDEBUG
  size_t bad_tape = tape.CompareChange();
  if (bad_tape != 0) {
    Rcpp::Rcout << "Warning:  " << bad_tape << " tape operations incorrect.  Must retape\n";
  }
#endif
  
  return(wrap(f)); 
}

template<typename TM>
NumericVector MB_Base<TM>::get_f_direct(const NumericVector& P_) {

 if (!tape_ready)
    throw MyException("tape not ready",__FILE__,__LINE__);


  VectorXS Pd = VectorXS::Map(P_.begin(),nvars);
  VectorXA P = Pd.cast<AScalar>();
  AScalar fr = model->eval_f(P); // inherited from base class
  VectorXS f(1);
  f(0) = Value(fr);
  return(wrap(f));
}

template<typename TM>
double MB_Base<TM>::get_LL(const NumericVector& P_) {

 if (!tape_ready)
    throw MyException("tape not ready",__FILE__,__LINE__);


  Map<VectorXS> Pd(P_.begin(),nvars);
  VectorXA P = Pd.cast<AScalar>();
  AScalar fr = model->eval_LL(P); // inherited from base class
  Scalar res = Value(fr);
  return(res);

}


template<typename TM>
NumericVector MB_Base<TM>::get_df(const NumericVector& P_) {

  if (!tape_ready)
    throw MyException("tape not ready",__FILE__,__LINE__);


  // need to copy to be sure P is a SimpleVector for CppAD
  VectorXS P = VectorXS::Map(P_.begin(),nvars); 
  VectorXS df = tape.Jacobian(P);
  return(wrap(df));
}

template<typename TM>
List MB_Base<TM>::get_fdf(const NumericVector& P_) {

  if (!tape_ready)
    throw MyException("tape not ready",__FILE__,__LINE__);

  VectorXS P = VectorXS::Map(P_.begin(),nvars); 
  VectorXS f(1);
  f = tape.Forward(0,P);
  VectorXS df = tape.Jacobian(P);
  
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

  VectorXS P = VectorXS::Map(P_.begin(),nvars); 
  VectorXS f(1);
  f = tape.Forward(0,P);
  VectorXS df = tape.Jacobian(P);

  MatrixXS hess = tape.Hessian(P,size_t(0));
  hess.conservativeResize(nvars, nvars);
 
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

  VectorXS P = VectorXS::Map(P_.begin(), nvars);
  MatrixXS hess = tape.Hessian(P, size_t(0));
  hess.conservativeResize(nvars, nvars);
  return(Rcpp::wrap(hess));

}

#ifdef COMPILE_SPARSE_HESSIAN

template<typename TM>
void MB_Base<TM>::hessian_init_nopattern(const NumericVector& P_)
{
			
// Initialization when sparsity pattern is unknown


  VectorXS P = VectorXS::Map(P_.begin(), nvars); 
  hess_info.clear();
 
  //  cout << "Identity matrix\n";

  std::vector<std::set<size_t> > identSet(nvars);

  for (size_t ii=0; ii<nvars; ii++) {
    identSet[ii].insert(ii);
  }
 
  // cout << "Computing sparsity pattern\n";
  tape.ForSparseJac(nvars, identSet);
  sp_type s(1);
  s[0].insert(0);
  sp = tape.RevSparseHes(nvars, s);

  // cout << "Computing nnz for one triangle
  int nnz_all = 0;
  for (int ii=0; ii<nvars; ii++) {
    nnz_all += sp[ii].size();
  }
  nnz = (nnz_all + nvars)/2;
  
  // extracting row and column indices for lower triangle of Hessian
  iRow.resize(nnz);
  jCol.resize(nnz);
  int idx = 0;
  for (size_t ii=0; ii<sp.size(); ii++) {
    for (std::set<size_t>::iterator jj=sp[ii].begin(); jj != sp[ii].end(); ++jj) {
      if (*jj <= ii) {
	iRow(idx) = ii;
	jCol(idx) = *jj;
	idx++;
      }
    }
  }
  

  // cout << "Computing hessian values\n"; 
  VectorXd hess_vals(nnz);
  VectorXd w(1);
  w(0)=1;
  size_t n_sweep = tape.SparseHessian(P,w,sp, iRow, jCol, hess_vals, hess_info);

  hessian_initialized = true;

}


template<typename TM>
void MB_Base<TM>::hessian_init(const NumericVector& P_,
			       const IntegerVector& hess_iRow,
			       const IntegerVector& hess_jCol)
{
  
  
  VectorXS P = VectorXS::Map(P_.begin(), nvars); 
  hess_info.clear();

  nnz = hess_iRow.size();
  iRow.resize(nnz);
  jCol.resize(nnz);

  // Construct sparsity pattern
  sp.resize(nvars);
  for (int k=0; k<nnz; k++) {
    iRow(k) = hess_iRow(k)-1;
    jCol(k) = hess_jCol(k)-1;
    sp[iRow[k]].insert(jCol[k]);
    sp[jCol[k]].insert(iRow[k]);
  }

  VectorXd w(1);
  w(0)=1;
  VectorXd hess_vals(nnz);
  size_t n_sweep = tape.SparseHessian(P,w,sp, iRow, jCol, hess_vals, hess_info);

  hessian_initialized = true;
 
}


template<typename TM>
Rcpp::S4 MB_Base<TM>::get_hessian_sparse(const NumericVector& P_)
{
			
  VectorXS P = VectorXS::Map(P_.begin(), nvars); 
 
  if (hessian_initialized) {

    // cout << "Computing hessian values\n"; 
    VectorXd hess_vals(nnz);
    VectorXd w(1);
    w(0)=1;
    size_t n_sweep = tape.SparseHessian(P,w,sp, iRow, jCol, hess_vals, hess_info);
    
    SparseMatrix<double> H(nvars, nvars);
    
    // typedef Eigen::Triplet<double> TT;
    // std::vector<TT> tripletList;
    // tripletList.resize(nnz);
    // for (int k=0; k<nnz; k++) {
    //   TT tmp(iRow(k), jCol(k), hess_vals(k));
    //   tripletList.push_back(tmp);
    // }
    // H.setFromTriplets(tripletList.begin(), tripletList.end());
    
    H.reserve(nnz);
    for (int k=0; k<nnz; k++) {
      H.insert(iRow(k), jCol(k)) = hess_vals(k);
    }
    H.makeCompressed();
    
    Eigen::SparseMatrix<double> res(nvars, nvars);
    res = H.selfadjointView<Eigen::Lower>();
    return(Rcpp::wrap(res));

  } else {
    throw MyException("Hessian not initialized",__FILE__, __LINE__); 
  }
}

#endif






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

#endif




