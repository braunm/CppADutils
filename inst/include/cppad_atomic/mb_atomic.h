#ifndef __MB_ATOMIC
#define __MB_ATOMIC

#define EIGEN_MATRIX_PLUGIN <eigen_plugin.h>
#define EIGEN_SPARSEMATRIX_PLUGIN <eigen_sparse_plugin.h>

#include <RcppEigen.h>
#include <Eigen/Core>
#include <cppad/cppad.hpp>

using Eigen::Dynamic;
using Eigen::MatrixBase;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;

using CppAD::vector;

typedef CppAD::AD<double> AScalar;
typedef Eigen::Matrix<AScalar, Dynamic, 1> VectorXA;
typedef Eigen::Matrix<AScalar, Dynamic, Dynamic> MatrixXA;

template<typename TF>
class mb_atomic : public CppAD::atomic_base<double> {
       
 public:
  std::shared_ptr<TF> func;
      
 mb_atomic(const std::string& name) :
  CppAD::atomic_base<double>(name) {
    func = std::make_shared<TF>();
  }

  ~mb_atomic() {}
      
 private:

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
  
      
  // ----------------------------------------------------------------------
  // forward mode routine called by CppAD
  virtual bool forward(
		       size_t                p , // lowest order Taylor coeff
		       size_t                q , // highest order Taylor coeff
		       const vector<bool>&   vx , // no idea
		       vector<bool>&         vy , // no idea
		       const vector<double>& tx , // input vector of x
		       vector<double>&       ty   // output vector of y
		       )
  {
	
    size_t n = tx.size() / (p+1);
	
    const Map<const MatrixXd> x = MatrixXd::Map(&(tx[0]),p+1,n);
    Map<VectorXd> y = VectorXd::Map(&(ty[0]),p+1);
	
    bool ok = (p <= 2) && (q <= p);
	
    if( vx.size() > 0 ) {
      vy[0] = vx[0];
      for (size_t j=1; j<n; j++) {
	vy[0] = vy[0] || vx[j];
      }
    }
    
    if ((q <= 0) && (p == 0)) {
      double f;
      func->eval(x.row(0),f);
      y(0) = f;      
    }
	
    if ((q <= 1) && (p == 1)) {
      double f;
      VectorXd df(n);
      func->eval(x.row(0),f,df);
      y(0) = f;      
      y(1) = df.dot(x.row(1));    
    }
	
    if ((q <= 2) && (p == 2)) {
      double f;
      VectorXd df(n);
      MatrixXd hess(n,n);
      func->eval(x.row(0), f, df, hess);
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
		       size_t                 q ,
		       const vector<double>&  tx ,
		       const vector<double>&  ty ,
		       vector<double>&        px ,
		       const vector<double>&  py
		       )
  {		
    size_t n = tx.size() / (q+1);
    
    const Map<const MatrixXd> x = MatrixXd::Map(&(tx[0]), q+1, n);
	
    double f;
    VectorXd df(n);     
    MatrixXd dy(n, q+1);

    bool ok = (q <= 2);
	
    if (q == 0) {
      func->eval(x.row(0), f, df);
      dy.col(0) = df;
    }
    
    if (q >= 1) {
      MatrixXd hess(n,n);
      func->eval(x.row(0), f, df, hess);
      dy.col(0) = df;
      dy.col(1) = hess * x.row(1).transpose();
    }
	
    for (size_t j=0; j<n; j++){
      for (size_t L=0; L <= q; L++) {
	px[j*(q+1)+L] = 0;
	for (size_t k=L; k<=q; k++) {
	  px[j*(q+1)+L] += py[k]*dy(j, k-L);
	}
      }
    }	
    return ok;
  }
      
   
  // ----------------------------------------------------------------------
  // forward Jacobian sparsity routine called by CppAD
  virtual bool for_sparse_jac(
			      size_t                             q ,
			      const vector<std::set<size_t> >&   r , 
			      vector<std::set<size_t> >&         s
			      ) {
    size_t n = r.size();
    s[0] = r[0];
    for (size_t i=0; i<n; i++) {
      my_union(s[0], s[0], r[i]);
    }  
    
    return true; 
  }
      
  
  // ----------------------------------------------------------------------
  // reverse Jacobian sparsity routine called by CppAD
  virtual bool rev_sparse_jac(
			      size_t  p ,
			      const vector< std::set<size_t> >& rt , 
			      vector< std::set<size_t> >& st 
			      ) {

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
#endif



  
