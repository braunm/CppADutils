#ifndef __LOG1PMX_ATOMIC
#define __LOG1PMX_ATOMIC

#define EIGEN_MATRIX_PLUGIN <eigen_plugin.h>
#define EIGEN_SPARSEMATRIX_PLUGIN <eigen_sparse_plugin.h>

#include <RcppEigen.h>
#include <Eigen/Core>
#include <cppad/cppad.hpp>
#include "union.h"

using Eigen::Dynamic;
using Eigen::MatrixBase;
using Eigen::Map;
using CppAD::vector;

typedef CppAD::AD<double> AScalar;
typedef Eigen::Matrix<AScalar, Dynamic, 1> VectorXA;
typedef Eigen::Matrix<AScalar, Dynamic, Dynamic> MatrixXA;

namespace CppAD {
  namespace Atomic{
    
    class atomic_log1pmx : public CppAD::atomic_base<double> {
       
    public:
    atomic_log1pmx(const std::string& name) :
      CppAD::atomic_base<double>(name) {}
      
    private:
      
      template<typename TX>
	bool log1pmx_deriv(const MatrixBase<TX>& x,
			 double& f){
 
	f = R::log1pmx(x(0));
	return true;	
      }
      
      template<typename TX, typename TD>
	bool log1pmx_deriv(const MatrixBase<TX>& x, 
			 double& f,
			 const MatrixBase<TD>& df_)	
      {
	MatrixBase<TD>& df = const_cast<MatrixBase<TD>&>(df_);       
	
	f = R::log1pmx(x(0));
	df(0) = 1.0 / (1.0 + x(0)) - 1.0;
	return true;
      }

      template<typename TX, typename TD, typename TH>
	bool log1pmx_deriv(const MatrixBase<TX>& x, 
			 double& f,
			 const MatrixBase<TD>& df_,
			 const MatrixBase<TH>& hess_)
	
      {	
	MatrixBase<TD>& df = const_cast<MatrixBase<TD>&>(df_);
	MatrixBase<TH>& hess = const_cast<MatrixBase<TH>&>(hess_);
	
	f = R::log1pmx(x(0));
	df(0) = 1.0 / (1.0 + x(0)) - 1.0;
	hess(0,0) = -1.0/(x(0)*x(0) + 2*x(0) + 1.0);
	
	return true;
      }
      
      
      
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
	
	using Eigen::MatrixXd;
	using Eigen::VectorXd;
	
	size_t p1 = p+1;
	size_t n = tx.size() / p1;
	//	size_t m = ty.size() / p1;
	
	double f;
	
	const Map<const MatrixXd> x = MatrixXd::Map(&(tx[0]),p1,n);
	Map<VectorXd> y = VectorXd::Map(&(ty[0]),p1);
	
	bool ok = (p <= 2) && (q <= p);
	
	if( vx.size() > 0 )
	  vy[0] = vx[0];
	
	
	if ((q <= 0) && (p == 0)) {
	  log1pmx_deriv(x.row(0),f);
	  y(0) = f;      
	}
	
	if ((q <= 1) && (p == 1)) {
	  VectorXd df(n);
	  log1pmx_deriv(x.row(0),f,df);
	  y(0) = f;      
	  y(1) = df.dot(x.row(1));    
	}
	
	if ((q <= 2) && (p == 2)) {
	  VectorXd df(n);
	  MatrixXd hess(n,n);
	  log1pmx_deriv(x.row(0),f,df,hess);
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
	
	using Eigen::MatrixXd;
	using Eigen::VectorXd;
	
	size_t n = tx.size() / (p+1);
	//	size_t m = ty.size() / (p+1);
	//	assert (m==1);
	
	const Map<const MatrixXd> x = MatrixXd::Map(&(tx[0]),p+1,n);
	const Map<const VectorXd> py = VectorXd::Map(&(py_[0]),p+1);
	Map<MatrixXd> px = MatrixXd::Map(&(px_[0]),p+1,n);
	
	double f;
	VectorXd df(n);     
	MatrixXd dy(n,p+1);

	bool ok = (p <= 2);
	
	if (p == 0) {
	  log1pmx_deriv(x.row(0), f, df);
	  dy.col(0) = df;
	}
	if (p >= 1) {
	  MatrixXd hess(n,n);
	  log1pmx_deriv(x.row(0), f, df, hess);
	  dy.col(0) = df;
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
				  size_t                                q ,
				  const vector<bool>&     r , // r.size() >= n
				  vector<bool>&           s  // s.size() >= m
				  ) {
	s[0] = r[0];
	
	return true;
	
      }

      // ----------------------------------------------------------------------
      // forward Jacobian sparsity routine called by CppAD
      virtual bool for_sparse_jac(
				  size_t                               q ,
				  const vector<std::set<size_t> >&     r , // r.size() >= n
				  vector<std::set<size_t> >&           s  // s.size() >= m
				  ) {
	s[0] = r[0];
	
	return true;
	
      }
      
      // ----------------------------------------------------------------------
      // reverse Jacobian sparsity routine called by CppAD
      virtual bool rev_sparse_jac(
				  size_t  q ,
				  const vector<bool>&  r , 
				  vector<bool>& s 
				  ) {
	
	for (size_t i=0; i<q; i++)
	  s[i] = r[i];
	
	return true;
      }
      
      
      // ----------------------------------------------------------------------
      // reverse Jacobian sparsity routine called by CppAD
      virtual bool rev_sparse_jac(
				  size_t  q ,
				  const vector< std::set<size_t> >& rt , 
				  vector< std::set<size_t> >& st 
				  ) {
	
	st[0] = rt[0];
	
	return true;
      }
      

      
      // ----------------------------------------------------------------------
      // reverse Hessian sparsity routine called by CppAD
      virtual bool rev_sparse_hes(
				  const vector<bool>&             vx, 
				  const vector<bool>&             s, 
				  vector<bool>&                   t, 
				  size_t                          q ,             
				  const vector<bool>&                  r ,
				  const vector<bool>&                  u ,
				  vector<bool>&                         v
				  ) 
      {
	
	t[0] = s[0];
	
	size_t j;
	for(j = 0; j < q; j++)
	  v[j] = u[j];
	
	if( s[0] ) {
	  for(j = 0; j < q; j++) {
	    v[j] |= r[j];
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
				  size_t                          q ,             
				  const vector<std::set<size_t> >&  r ,
				  const vector<std::set<size_t> >&  u ,
				  vector<std::set<size_t> >&  v
				  ) 
      {
	
	t[0] = s[0];
	v[0] = u[0];

	if (s[0]) 
	  my_union(v[0], v[0], r[0]);
		
	return true; 
	
      }

    }; // end class

 
    AScalar log1pmx(const AScalar& a) {
      
      static atomic_log1pmx func("atomic_log1pmx");
      
      VectorXA x(1);
      VectorXA y(1);
      
      x(0) = a;
      func(x,y);
      return(y[0]);
    }
    
    
  } // end namespace atomic
} // end namespace CppAD

#endif
