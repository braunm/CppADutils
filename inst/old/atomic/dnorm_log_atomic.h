#ifndef __DNORM_LOG_ATOMIC
#define __DNORM_LOG_ATOMIC

#define EIGEN_MATRIX_PLUGIN <eigen_plugin.h>
#define EIGEN_SPARSEMATRIX_PLUGIN <eigen_sparse_plugin.h>

#include <RcppEigen.h>
#include <Eigen/Core>
#include <cppad/cppad.hpp>
#include "union.h"

using Eigen::Dynamic;
using Eigen::MatrixBase;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;

using CppAD::vector;

typedef CppAD::AD<double> AScalar;
typedef Eigen::Matrix<AScalar, Dynamic, 1> VectorXA;
typedef Eigen::Matrix<AScalar, Dynamic, Dynamic> MatrixXA;

namespace CppAD { 
  namespace Atomic{
        
    class atomic_dnorm_log : public CppAD::atomic_base<double> {
       
    public:
      atomic_dnorm_log(const std::string& name) :
	CppAD::atomic_base<double>(name) {}
      
    private:
      
      template<typename TX>
	bool dnorm_log_deriv(const MatrixBase<TX>& x,
			 double& f)
      
      {	
	f = R::dnorm(x(0), x(1), x(2), 1);  // log p	
	return true;
      }

      template<typename TX, typename TD>
	bool dnorm_log_deriv(const MatrixBase<TX>& x, 
			     double& f,
			     const MatrixBase<TD>& df_)	
      {
	
	MatrixBase<TD>& df = const_cast<MatrixBase<TD>&>(df_);

	const double y = x(0);
	const double m = x(1);
	const double s = x(2);
	
	f = R::dnorm(y, m, s, 1);  // log p
	df(0) = (m-y) / (s*s);
	df(1) = -df(0);
	df(2) = (y-m-s) * (y-m+s) / (s*s*s);
	
	
	return true;	
	
      }

      template<typename TX, typename TD, typename TH>
	bool dnorm_log_deriv(const MatrixBase<TX>& x, 
			     double& f,
			     const MatrixBase<TD>& df_,
			     const MatrixBase<TH>& hess_)
	
      {
	
	MatrixBase<TD>& df = const_cast<MatrixBase<TD>&>(df_);
	MatrixBase<TH>& hess = const_cast<MatrixBase<TH>&>(hess_);
	
	const double y = x(0);
	const double m = x(1);
	const double s = x(2);
	
	f = R::dnorm(y, m, s, 1);  // log p
	df(0) = (m-y) / (s*s);
	df(1) = -df(0);
	df(2) = (y-m-s) * (y-m+s) / (s*s*s);
	hess(1,1) = -1/(s*s);
	hess(0,0) = hess(1,1);
	hess(0,1) = -hess(1,1);
	hess(1,0) = -hess(1,1);
	hess(0,2) = 2*df(1)/s;
	hess(2,0) = hess(0,2);
	hess(1,2) = -hess(0,2);
	hess(2,1) = hess(1,2);
	hess(2,2) = hess(0,1) - 3*df(0)*df(0);
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
	
	
	size_t p1 = p+1;
	size_t n = tx.size() / p1;

	double f;
	
	const Map<const MatrixXd> x = MatrixXd::Map(&(tx[0]),p1,n);
	Map<VectorXd> y = VectorXd::Map(&(ty[0]),p1);
	
	bool ok = (p <= 2) && (q <= p);
	
	if( vx.size() > 0 )
	  {
	    vy[0] = vx[0] || vx[1] || vx[2];
	  }
	
	if ((q <= 0) && (p == 0)) {
	  dnorm_log_deriv(x.row(0),f);
	  y(0) = f; 
	  ok = true;     
	}
	
	if ((q <= 1) && ((p == 1) || (p == 2))) {
	  VectorXd df(n); // make these Eigen static allocation?
	  MatrixXd hess(n,n);
	  dnorm_log_deriv(x.row(0), f, df, hess);
	  y(0) = f;      
	  y(1) = df.dot(x.row(1));    
	  if (p==2) {
	    y(2) = x.row(1) * hess * x.row(1).transpose();    
	    y(2) *= 0.5;
	    y(2) += df.dot(x.row(2));
	  }
	  ok = true;
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
	
	const Map<const MatrixXd> x = MatrixXd::Map(&(tx[0]),p+1,n);
	const Map<const VectorXd> py = VectorXd::Map(&(py_[0]),p+1);
	Map<MatrixXd> px = MatrixXd::Map(&(px_[0]),p+1,n);
	
	double f;
	VectorXd df(n);     
	MatrixXd dy(n,p+1);
	MatrixXd hess(n,n);

	bool ok = (p <= 2);
	
	dnorm_log_deriv(x.row(0), f, df, hess);
	dy.col(0) = df;
	if (p >= 1) {
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

	size_t n = r.size() / q;

	for (size_t i=0; i<n; i++) {
	  my_union(s[0], s[0], r[i]);
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
				  size_t                          q ,             
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
    
    AScalar dnorm_log(const AScalar& z, const AScalar& a, const AScalar& b) {
      
      static atomic_dnorm_log func("atomic_dnorm_log");
      
      VectorXA x(3);
      VectorXA y(1);
      
      x(0) = z;
      x(1) = a;
      x(2) = b;
      func(x,y);
      return(y[0]);
    }
    
    
  } // end namespace atomic
} // end namespace CppAD

#endif
