
#define EIGEN_MATRIX_PLUGIN "/Users/braunm/ownCloud/CppADutils/inst/include/eigen_plugin.h"
#define EIGEN_SPARSEMATRIX_PLUGIN "/Users/braunm/ownCloud/CppADutils/inst/include/eigen_sparse_plugin.h" 

// undefine this to record function directly
#define MB_USE_ATOMIC

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <cppad/cppad.hpp>

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;

using CppAD::vector;

typedef CppAD::AD<double> AScalar;
typedef Eigen::Matrix<AScalar, Eigen::Dynamic, 1> VectorXA;
typedef Eigen::Matrix<AScalar, Eigen::Dynamic, Eigen::Dynamic> MatrixXA;
typedef Eigen::Triplet <double> TT;
typedef Eigen::SparseMatrix<double> SparseMatrixXd;

inline void my_union( 
		     std::set<size_t>&         result  , 
		     const std::set<size_t>&   left    , 
		     const std::set<size_t>&   right   ) 
{  
  std::set<size_t> temp; 
  std::set_union( 
		 left.begin()              , 
		 left.end()                , 
		 right.begin()             , 
		 right.end()               , 
		 std::inserter(temp, temp.begin()) 
		  ); 
  result.swap(temp); 
}

class atomic_cross : public CppAD::atomic_base<double> {
public:
  atomic_cross(const std::string& name) :
    CppAD::atomic_base<double>(name) {}
  
private:
  
  // ----------------------------------------------------------------------
  // forward mode routine called by CppAD
  virtual bool forward(
		       size_t                   q , // lowest order Taylor coeff
		       size_t                    p , // highest order Taylor coeff
		       const vector<bool>&      vx , // no idea
		       vector<bool>&            vy , // no idea
		       const vector<double>&     tx , // input vector of x
		       vector<double>&           ty   // output vector of y
		       )
  {
       
    size_t p1 = p+1;
    size_t n = tx.size() / p1;
    size_t m = ty.size() / p1;

    double f;
    
    const Map<const MatrixXd> x = MatrixXd::Map(&(tx[0]),p1,n);
    Map<VectorXd> y = VectorXd::Map(&(ty[0]),p1);

    const double a = x(0);
    const double b = x(1);
    
    bool ok = (p <= 2) && (q <= p);
    
    if( vx.size() > 0 )
      {
	vy[0] = vx[0] || vx[1];
      }
    
    if ((q <= 0) && (p == 0)) {
      y(0) = a*a*b*b;      
    }
    
    if ((q <= 1) && (p == 1)) {
      y(0) = f;      
      VectorXd df(n);
      df(0) = 2*a*b*b;
      df(1) = 2*a*a*b;
      
      y(1) = df.dot(x.row(1));    
    }
    
    if ((q <= 2) && (p == 2)) {
      VectorXd df(n);
      df(0) = 2*a*b*b;
      df(1) = 2*a*a*b;	  
      MatrixXd hess(n,n);
      hess(0,0) = 2*b*b;
      hess(1,1) = 2*a*a;
      hess(0,1) = 4*a*b;
      hess(1,0) = 4*a*b;
      
      y(0) = f;      
      y(1) = df.dot(x.row(1)); 	
      y(2) = x.row(1) * hess * x.row(1).transpose();    
      y(2) *= 0.5;
      y(2) += df.dot(x.row(2));
    }
    
    return ok;
  }
  
  // ----------------------------------------------------------------------
  // reverse mode routine called by CppAD
  virtual bool reverse(
		       size_t                   p ,
		       const vector<double>&     tx ,
		       const vector<double>&     ty ,
		       vector<double>&           px_ ,
		       const vector<double>&     py_
		       )
  {
    
    size_t n = tx.size() / (p+1);
    size_t m = ty.size() / (p+1);
    
    const Map<const MatrixXd> x = MatrixXd::Map(&(tx[0]),p+1,n);
    const Map<const VectorXd> py = VectorXd::Map(&(py_[0]),p+1);
    Map<MatrixXd> px = MatrixXd::Map(&(px_[0]),p+1,n);

    const double a = x(0,0);
    const double b = x(0,1);
    
    double f;
    VectorXd df(n);     
    MatrixXd dy(n,p+1);
    
    bool ok = (p <= 2);
    
    if (p == 0) {
      dy(0,0) = 2*a*b*b;
      dy(1,0) = 2*a*a*b;
    }
    if (p >= 1) {
      dy(0,0) = 2*a*b*b;
      dy(1,0) = 2*a*a*b;

      MatrixXd hess(n,n);
      hess(0,0) = 2*b*b;
      hess(1,1) = 2*a*a;
      hess(0,1) = 4*a*b;
      hess(1,0) = 4*a*b; 
      dy.col(1) = hess * x.row(1).transpose();
    }
	
    px.setZero();   
    
    for (int j=0; j<n; j++){
      for (int k=0; k<=p; k++) {
	for (int i=k; i<=p; i++) {
	  px(k,j) += py(i) * dy(j,i-k);
	}
      }
    }
    
    return ok;
  }
  
  // ----------------------------------------------------------------------
  // forward Jacobian sparsity routine called by CppAD
  virtual bool for_sparse_jac(
			      size_t                               q ,
			      const vector<std::set<size_t> >&     r , // r.size() >= n
			      vector<std::set<size_t> >&           s  // s.size() >= m
			      ) {
    size_t n = r.size();
    if (n > 1) {
      for (size_t i=0; i<n; i++) {
	my_union(s[0], r[0], r[i]);
      }
    } else {
      s[0] = r[0];
    }
      
    return true; 
  }
  
  
  // ----------------------------------------------------------------------
  // reverse Jacobian sparsity routine called by CppAD
  virtual bool rev_sparse_jac(
			      size_t  q ,
			      const vector< std::set<size_t> >& rt , 
			      vector< std::set<size_t> >& st 
			      ) {
    // m == 1
    size_t n = st.size();
    for (size_t i=0; i<n; i++) {
      st[i] = rt[0];
    }
    
    return true; 
  }

     // ----------------------------------------------------------------------
      // reverse Hessian sparsity routine called by CppAD
      virtual bool rev_sparse_hes(
				  const vector<bool>&             vx, 
				  const vector<bool>&             s, 
				  vector<bool>&                   t, 
				  size_t                          p ,             
				  const vector<std::set<size_t> >&  r ,
				  const vector<std::set<size_t> >&  u ,
				  vector<std::set<size_t> >&  v
				  ) 
      {
	size_t n = vx.size();

	for (size_t i=0; i<n; i++) { 
	  t[i] = s[0];
	  v[i] = u[0];
	  if( s[0] ) {
	    for (size_t j = 0; j < n; j++) {
	      my_union(v[i], v[i], r[j] );
	    }
	  }
	}	

	return true;	
      }
}; // end class


AScalar cross(const AScalar& a, const AScalar& b) {
  
  static atomic_cross func("atomic_cross");
  
  VectorXA x(2);
  VectorXA y(1);
  x(0) = a;
  x(1) = b;
  func(x,y);
  return(y[0]);
}

int main() {
  
  typedef vector< std::set<size_t> > sp_set;
  
  CppAD::ADFun<double> tape_set;
  CppAD::sparse_hessian_work hess_info_set;
  
  int nvars = 2;
  VectorXA f(1); // to hold result

  const double a = 2;
  const double b = 3;
  
  VectorXd X(nvars);
  X(0) = a;
  X(1) = b;
  VectorXA P = X.cast<AScalar>();        
  
  // correct values
  
  double trueval = a*a*b*b;
  VectorXd truegrad(2);
  truegrad(0) = 2*a*b*b;
  truegrad(1) = 2*a*a*b;
  
  MatrixXd truehess(2,2);
  truehess(0,0) = 2*b*b;
  truehess(1,1) = 2*a*a;
  truehess(1,0) = 4*a*b;
  truehess(0,1) = truehess(1,0);
  
  std::cout << "True value:\n" << trueval << "\n\n";
  std::cout << "True gradient:\n" << truegrad << "\n\n";
  std::cout << "True Hessian:\n" << truehess << "\n\n";
  
  // True Hessian sparsity structure
  vector<std::set<size_t> > pset;
  pset.resize(nvars);
  for (size_t i=0; i<nvars; i++) {
    for (size_t j=0; j<nvars; j++) {
      pset[i].insert(j);
    }
  }  
  
  // record tape for testing set structure
  CppAD::Independent(P);
#ifdef MB_USE_ATOMIC
  f(0) = cross(P(0), P(1));
#else
  f(0) = P(0)*P(0)*P(1)*P(1);
#endif
  
  tape_set.Dependent(P, f);
  tape_set.optimize();   
  
 
  VectorXd val_set(1);
  val_set = tape_set.Forward(0, X);
  VectorXd grad_set = tape_set.Jacobian(X);
  std::cout << "Gradient:\n" << grad_set << "\n\n";
  
  MatrixXd hess_dense_set = tape_set.Hessian(X, size_t(0));
  hess_dense_set.conservativeResize(nvars, nvars);

  std::cout << "Dense Hessian:\n" << hess_dense_set << "\n\n";
  
  VectorXd w(1);
  w(0) = 1;
  // make identity matrices
  sp_set identSet(nvars); 
  for (size_t ii=0; ii<nvars; ii++) {
    identSet[ii].insert(ii);
  }
    
  // sparse, with set structure
  MatrixXd hess_sparse_set = tape_set.SparseHessian(X, w, pset);
  hess_sparse_set.conservativeResize(nvars, nvars);

  bool test_val = trueval == val_set(0);
  bool test_grad =  truegrad.cwiseEqual(grad_set).all();
  bool test_hess_dense =  truehess.cwiseEqual(hess_dense_set).all();
  bool test_hess_sparse_driver = truehess.cwiseEqual(hess_sparse_set).all(); 



  
  std::cout << "Tests:  pass = 1, fail = 0\n";
  std::cout << "Function value: " << test_val << "\n";
  std::cout << "Gradient: " << test_grad << "\n";
  std::cout << "Dense Hessian: " << test_hess_dense << "\n";
  std::cout << "Sparse Hessian, computed from driver: " << test_hess_sparse_driver << "\n\n";

  std::cout << "All of the tests below fail when MB_USE_ATOMIC is defined,\n";
  std::cout << "  but pass when it is not.\n\n";

  /////////////////
  // THIS IS WHERE THE PROBLEMS START TO SHOW UP
  
  // forward Jacobian structure
   sp_set sp = tape_set.ForSparseJac(nvars, identSet);
  
  // std::cout << "Jacobian structure, set:\n";
  // for (auto jj=sp[0].begin(); jj != sp[0].end(); ++jj) {
  //   std::cout << *jj << " ";
  // }
  // std::cout << "\n";
  // std::cout << "Should be 0 1 (not sure how to compare sets)\n\n";
   
  // reverse Hessian structure 
  sp_set zset(1);
  zset[0].insert(0);
  sp_set hs_set = tape_set.RevSparseHes(nvars, zset);
  
  std::cout << "Hessian structure, set:\n";
  for (size_t ii=0; ii<nvars; ii++) {
    for (auto jj=hs_set[ii].begin(); jj != hs_set[ii].end(); ++jj) {
      std::cout << *jj << " ";
    }
    std::cout << "\n";
  }
  std::cout << "Should be:  (not sure how to compare sets)\n" << "0  1\n0  1\n\n";

  //
  //////////////////

  
  /* THE REST OF THE CODE IS ALL OF THE STEPS NEEDED
     TO RETURN THE LOWER TRIANGLE OF A SPARSE HESSIAN
  */
  
  // std::cout << "Computing nnz for one triangle
  int nnz_all_set = 0;
  for (int ii=0; ii<nvars; ii++) {
    nnz_all_set += hs_set[ii].size();
  }
  int nnz_LT_set = (nnz_all_set + nvars)/2;
  bool test_nnz_LT = nnz_LT_set == 3;
  
  // Extracting row/col indices for lower triangle
  VectorXi iRow_set(nnz_LT_set);
  VectorXi jCol_set(nnz_LT_set);
  size_t idx = 0;
  for (size_t ii=0; ii<nvars; ii++) {
    for (std::set<size_t>::iterator jj=hs_set[ii].begin(); jj != hs_set[ii].end(); ++jj) {
      if (*jj <= ii) {
	iRow_set(idx) = ii;
	jCol_set(idx) = *jj;
	idx++;
      }
    }
  }
  
  // Get LT elements of sparse Hessian
  VectorXd hess_vals_set(nnz_LT_set);
  size_t n_sweep_set = tape_set.SparseHessian(X, w, hs_set, iRow_set, jCol_set, hess_vals_set, hess_info_set);
  
  // Make Eigen sparse matrix of LT of Hessian
  std::vector<TT> trips_set;
  trips_set.resize(nnz_LT_set);
  for (int k=0; k<nnz_LT_set; k++) {
    TT tmp = TT(iRow_set(k), jCol_set(k), hess_vals_set(k));
    trips_set.push_back(tmp);
  }
  
  SparseMatrixXd H_LT_set;
  H_LT_set.resize(nvars, nvars);
  H_LT_set.setFromTriplets(trips_set.begin(), trips_set.end());
  H_LT_set.makeCompressed();
  
  // Compare LT of sparse Hessian to LT of true Hessian
  
  MatrixXd truehessLT = truehess.triangularView<Eigen::Lower>();
  MatrixXd sparsehessLT = H_LT_set.triangularView<Eigen::Lower>();
  bool test_sparse_LT = truehessLT.cwiseEqual(sparsehessLT).all();
  
  std::cout << "NNZ in sparse lower triangle (should be 3): " << nnz_LT_set << "\n";
  std::cout << "Sparse Hessian (Lower, using computed pattern)\n";
  std::cout <<  "pass = 1, fail = 0:  " << test_sparse_LT << "\n\n";
  std::cout << sparsehessLT << "\n\n"; 
  return 0;
}
