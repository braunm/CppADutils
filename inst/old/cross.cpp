
#define EIGEN_MATRIX_PLUGIN "/Users/braunm/ownCloud/CppADutils/inst/include/eigen_plugin.h"
#define EIGEN_SPARSEMATRIX_PLUGIN "/Users/braunm/ownCloud/CppADutils/inst/include/eigen_sparse_plugin.h" 


// undefine this to record function directly
#define MB_USE_ATOMIC

#define CPPAD_FORWARD1SWEEP_TRACE 1
#define CPPAD_FORWARD0SWEEP_TRACE 1
#define CPPAD_FOR_JAC_SWEEP_TRACE 1
#define CPPAD_REV_HES_SWEEP_TRACE 1


#include <Eigen/Core>
#include <Eigen/Sparse>
#include <cppad/cppad.hpp>



using Eigen::Dynamic;
using Eigen::MatrixBase;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;


using CppAD::vector;

typedef CppAD::AD<double> AScalar;
typedef Eigen::Matrix<AScalar, Dynamic, 1> VectorXA;
typedef Eigen::Matrix<AScalar, Dynamic, Dynamic> MatrixXA;
typedef Eigen::Matrix<bool, Dynamic, 1> VectorXb;
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
    
    assert(n==2);
    assert(m==1);
    
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
    assert (m==1);
    assert(n==2);
    
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
			      size_t                                p ,
			      const vector<bool>&     r , // r.size() >= n
			      vector<bool>&           s  // s.size() >= m
			      ) {
    size_t n = r.size() / p;
    size_t m = s.size() / p;
    assert( n == 2 );
    assert( m == 1 );
    
    // sparsity for S(x) = f'(x) * R 
    // where f'(x) = 2 * [ x_0, x_1 ]
    for(size_t j = 0; j < p; j++) {
      s[j] = true;
    }
    //    for(size_t i = 0; i < n; i++)
    //	s[j] |= r[i * p + j];
    // }
    return true; 
    
    
  }
  
  // ----------------------------------------------------------------------
  // forward Jacobian sparsity routine called by CppAD
  virtual bool for_sparse_jac(
			      size_t                               q ,
			      const vector<std::set<size_t> >&     r , // r.size() >= n
			      vector<std::set<size_t> >&           s  // s.size() >= m
			      ) {
  
    my_union(s[0], s[0], r[1]);
    
    return true; 
  }
  
  // ----------------------------------------------------------------------
  // reverse Jacobian sparsity routine called by CppAD
  virtual bool rev_sparse_jac(
			      size_t  q ,
			      const vector<bool>&  rt , 
			      vector<bool>& st 
			      ) {
    
    size_t n = st.size() / q;
    for(size_t j = 0; j < q; j++)
      for(size_t i = 0; i < n; i++)
	st[i * q + j] = rt[j];
    
    return true; 
  }
  
  // ----------------------------------------------------------------------
  // reverse Jacobian sparsity routine called by CppAD
  virtual bool rev_sparse_jac(
			      size_t  p ,
			      const vector< std::set<size_t> >& rt , 
			      vector< std::set<size_t> >& st 
			      ) {
    st[0] = rt[0];
    st[1] = rt[0];
    
    return true; 
  }
  
  
  // ----------------------------------------------------------------------
  // reverse Hessian sparsity routine called by CppAD
  virtual bool rev_sparse_hes(
			      const vector<bool>&             vx, 
			      const vector<bool>&             s, 
				  vector<bool>&                   t, 
				  size_t                          p ,             
				  const vector<bool>&                  r ,
				  const vector<bool>&                  u ,
				  vector<bool>&                         v
				  ) 
      {

	size_t n = s.size();
	
	t[0] = s[0];
	t[1] = s[0];
	
	for(size_t j = 0; j < p; j++) {
	  for(size_t i = 0; i < n; i++) {
	    v[ i * p + j] = u[j];
	  }
	}
	
	if( s[0] ) {
	  for(size_t j = 0; j < p; j++) {
	    for(size_t i = 0; i < n; i++) {
	      v[ i * p + j] |= r[ i * p + j];
	    }
	  }
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
	size_t m = s.size();

	t[0] = s[0];
	t[1] = s[0];

	v[0] = u[0];
	v[1] = u[0];
	
	if( s[0] ) {
	  for(size_t i = 0; i < n; i++) {
	      my_union(v[i], v[i], r[i] );
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
  
  typedef VectorXb sp_bool;
  typedef vector< std::set<size_t> > sp_set;
  
  
  CppAD::ADFun<double> tape_bool;
  CppAD::ADFun<double> tape_set;
  CppAD::sparse_hessian_work hess_info_bool;
  CppAD::sparse_hessian_work hess_info_set;
  
  int nvars = 2;
  VectorXA f(1); // to hold result

  const double a = 2;
  const double b = 3;
  
  VectorXd X(nvars);
  X(0) = a;
  X(1) = b;
  
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
  
  VectorXA P = X.cast<AScalar>();      
  
  // record tape for testing boolean structure
  CppAD::Independent(P);

#ifdef MB_USE_ATOMIC
  f(0) = cross(P(0), P(1));
#else
  f(0) = P(0)*P(0)*P(1)*P(1);
#endif
  
  tape_bool.Dependent(P, f);
  tape_bool.optimize();   
  
  
  // record tape for testing set structure
  CppAD::Independent(P);
#ifdef MB_USE_ATOMIC
  f(0) = cross(P(0), P(1));
#else
  f(0) = P(0)*P(0)*P(1)*P(1);
#endif
  
  tape_set.Dependent(P, f);
  tape_set.optimize();   
  
  
  
  VectorXd val_bool(1);
  val_bool = tape_bool.Forward(0, X);
  VectorXd grad_bool = tape_bool.Jacobian(X);
  MatrixXd hess_dense_bool = tape_bool.Hessian(X, size_t(0));
  hess_dense_bool.conservativeResize(nvars, nvars);

  VectorXd val_set(1);
  val_set = tape_set.Forward(0, X);
  VectorXd grad_set = tape_set.Jacobian(X);
  MatrixXd hess_dense_set = tape_set.Hessian(X, size_t(0));
  hess_dense_set.conservativeResize(nvars, nvars);


  std::cout << "True Value:  " << trueval << "\n";
  std::cout << "Value, tape_bool:  " << val_bool << "\n";
  std::cout << "Value, tape_set:  " << val_set << "\n\n";

  std::cout << "True grad:\n" << truegrad << "\n\n";
  std::cout << "grad, tape_bool:\n" << grad_bool << "\n\n";
  std::cout << "grad, tape_set:\n" << grad_set << "\n\n";
  
  
  
  
  VectorXd w(1);
  w(0) = 1;
  
  VectorXb pbool(nvars*nvars);
  pbool.setConstant(true);
  
  vector<std::set<size_t> > pset;
  pset.resize(nvars);
  for (size_t i=0; i<nvars; i++) {
    for (size_t j=0; j<nvars; j++) {
      pset[i].insert(j);
    }
  }
  
  
  // sparse, with boolean structure
  MatrixXd hess_sparse_bool = tape_bool.SparseHessian(X, w, pbool);
  hess_sparse_bool.conservativeResize(nvars, nvars);
  
  // sparse, with set structure
  MatrixXd hess_sparse_set = tape_set.SparseHessian(X, w, pset);
  hess_sparse_set.conservativeResize(nvars, nvars);


  

  std::cout << "True Hessian:\n" << truehess << "\n\n";
  std::cout << "Dense Hessian, tape_bool:\n" << hess_dense_bool << "\n\n";
  std::cout << "Dense Hessian, tape_set:\n" << hess_dense_set << "\n\n";
  std::cout << "SparseHessian, boolean structure:\n" << hess_sparse_bool << "\n\n";
  std::cout << "SparseHessian, set structure:\n" << hess_sparse_set << "\n\n";


  
  // make identity matrices
  
  sp_bool identBool(nvars*nvars);    
  identBool.setConstant(false);
  for (size_t ii=0; ii<nvars; ii++) {
    identBool(nvars*ii+ii) = true;
  }
  
  sp_set identSet(nvars); 
  for (size_t ii=0; ii<nvars; ii++) {
    identSet[ii].insert(ii);
  }
  
  // forward Jacobian structure
  sp_bool sb = tape_bool.ForSparseJac(nvars, identBool);
  sp_set sp = tape_set.ForSparseJac(nvars, identSet);

  std::cout << "Jacobian structure, bool:\n" << sb << "\n\n";

  std::cout << "Jacobian structure, set:\n";
    for (auto jj=sp[0].begin(); jj != sp[0].end(); ++jj) {
      std::cout << *jj << " ";
    }
    std::cout << "\n";
  std::cout << "\n\n";

  sp_bool zbool(1);
  zbool(0) = true;
  
  sp_bool hs_bool = tape_bool.RevSparseHes(nvars, zbool);

  std::cout << "Hessian structure, bool:\n" << hs_bool << "\n\n";
  
  int nnz_all_bool = hs_bool.count();
  int nnz_LT_bool = (nnz_all_bool + nvars)/2;

  std::cout << "nnz_bool = " << nnz_all_bool << "\n\n";
  std::cout << "nnz_LT_bool = " << nnz_LT_bool << "\n\n";
    
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
  std::cout << "\n\n";
  
  
  // std::cout << "Computing nnz for one triangle
  int nnz_all_set = 0;
  for (int ii=0; ii<nvars; ii++) {
    nnz_all_set += hs_set[ii].size();
  }
  int nnz_LT_set = (nnz_all_set + nvars)/2;

  std::cout << "nnz_set = " << nnz_all_set << "\n\n";
  
    
  VectorXi iRow_bool(nnz_LT_bool);
  VectorXi jCol_bool(nnz_LT_bool);
  int idx = 0;
  for (size_t jj=0; jj < nvars; jj++) {
    for (size_t ii=jj; ii < nvars; ii++) {
      if (hs_bool(nvars*jj+ii)) {
	iRow_bool(idx) = ii;
	jCol_bool(idx) = jj;
	idx++;
      }
    }
  }

  std::cout << "Computing hessian values\n"; 
  VectorXd hess_vals_bool(nnz_LT_bool);
  size_t n_sweep_bool = tape_bool.SparseHessian(X, w, hs_bool, iRow_bool, jCol_bool, hess_vals_bool, hess_info_bool);

  for (size_t k=0; k<nnz_LT_bool; k++) {
    std::cout << "i = " << iRow_bool(k) << "  j = " << jCol_bool(k) << "  val = " << hess_vals_bool(k) << "\n";
  }
  std::cout << "\n";
  
  std::vector<TT> trips_bool;
  trips_bool.resize(nnz_LT_bool);
  for (int k=0; k<nnz_LT_bool; k++) {
    TT tmp = TT(iRow_bool(k), jCol_bool(k), hess_vals_bool(k));
    trips_bool.push_back(tmp);
  }
  
  SparseMatrixXd H_LT_bool;
  H_LT_bool.resize(nvars, nvars);
  H_LT_bool.setFromTriplets(trips_bool.begin(), trips_bool.end());
  H_LT_bool.makeCompressed();
  
  std::cout << "Sparse hessian bool (Lower only):\n" << H_LT_bool << "\n\n";

  
  VectorXi iRow_set(nnz_LT_set);
  VectorXi jCol_set(nnz_LT_set);
  idx = 0;
  for (size_t ii=0; ii<nvars; ii++) {
    for (std::set<size_t>::iterator jj=hs_set[ii].begin(); jj != hs_set[ii].end(); ++jj) {
      if (*jj <= ii) {
	iRow_set(idx) = ii;
	jCol_set(idx) = *jj;
	idx++;
      }
    }
  }
  
  VectorXd hess_vals_set(nnz_LT_set);
  size_t n_sweep_set = tape_set.SparseHessian(X, w, hs_set, iRow_set, jCol_set, hess_vals_set, hess_info_set);

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
  
  std::cout << "Sparse hessian set (Lower only):\n" << H_LT_set << "\n\n";
  
  return 0;
}

    



