#ifndef __INC_GAMMA_ATOMIC
#define __INC_GAMMA_ATOMIC

#define EIGEN_MATRIX_PLUGIN <eigen_plugin.h>
#define EIGEN_SPARSEMATRIX_PLUGIN <eigen_sparse_plugin.h>

#include <RcppEigen.h>
#include <Eigen/Core>
#include <cppad/cppad.hpp>
#include "union.h"
#include "inc_gamma_deriv.h"

using Eigen::Dynamic;
using Eigen::MatrixBase;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;

using CppAD::vector;

typedef CppAD::AD<double> AScalar;
typedef Eigen::Matrix<AScalar, Dynamic, 1> VectorXA;
typedef Eigen::Matrix<AScalar, Dynamic, Dynamic> MatrixXA;


// NOTE:  This is a NORMALIZED incomplete gamma function.
//  Equivalent to the cdf of a gamma(z;r,1) distribution

namespace CppAD {  
  namespace Atomic{
        
    class atomic_inc_gamma : public CppAD::atomic_base<double> {
       
    public:
      atomic_inc_gamma(const std::string& name) :
	CppAD::atomic_base<double>(name) {}
      
    private:
      
      template<typename TX>
      bool inc_gamma_deriv(const MatrixBase<TX>& x,
			double& f)
	
      {	
	f = R::pgamma(x(0), x(1), 1, 1, 0);	
	return true;
      }
            
      template<typename TX, typename TD>
      bool inc_gamma_deriv(const MatrixBase<TX>& x, 
			double& f,
			const MatrixBase<TD>& df_)
	
      {
	MatrixBase<TD>& df = const_cast<MatrixBase<TD>&>(df_);

	double z = x(0);
	double r = x(1);	
	double D[6]; // output	
	double grlog = R::lgammafn(r);
	double gr1log = log(r) + grlog;
	double psir = R::digamma(r);
	double psir1 = psir + 1/r;
	double psidr = R::trigamma(r);
	double psidr1 = psidr - 1/(r*r);
	int ifault;
	
	digami(D, &z, &r, &grlog, &gr1log, &psir,
	       &psir1, &psidr, &psidr1, &ifault);
	
	f = D[5];	
	df(0) = exp((r-1)*log(z)-z-grlog);
	df(1) = D[2];
	
	return ifault;		
      }

      template<typename TX, typename TD, typename TH>
      bool inc_gamma_deriv(const MatrixBase<TX>& x, 
			double& f,
			const MatrixBase<TD>& df_,
			const MatrixBase<TH>& hess_)
	
      {
	MatrixBase<TD>& df = const_cast<MatrixBase<TD>&>(df_);
	MatrixBase<TH>& hess = const_cast<MatrixBase<TH>&>(hess_);
	
	double z = x(0);
	double r = x(1);	
	double log_z = log(z);
	double D[6]; // output	
	double grlog = R::lgammafn(r);
	double gr1log = log(r) + grlog;
	double psir = R::digamma(r);
	double psir1 = psir + 1/r;
	double psidr = R::trigamma(r);
	double psidr1 = psidr - 1/(r*r);
	int ifault;
	
	digami(D, &z, &r, &grlog, &gr1log, &psir,
	       &psir1, &psidr, &psidr1, &ifault);
	
	f = D[5];	
	df(0) = exp((r-1)*log(z)-z-grlog);
	df(1) = D[2];		
	hess(0,0) = exp((r-2)*log_z - z-grlog)*(r-z-1);
	hess(1,1) = D[3];
	hess(0,1) = exp((r-1)*log_z-z-grlog)*(log_z-psir);
	hess(1,0) = hess(0,1);
	return ifault;	
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

	double f;
	
	const Map<const MatrixXd> x = MatrixXd::Map(&(tx[0]),p1,n);
	Map<VectorXd> y = VectorXd::Map(&(ty[0]),p1);
	
	bool ok = (p <= 2) && (q <= p);
	
	if( vx.size() > 0 )
	  {
	    vy[0] = vx[0] || vx[1];
	  }
	
	if ((q <= 0) && (p == 0)) {
	  inc_gamma_deriv(x.row(0),f);
	  y(0) = f;      
	}
	
	if ((q <= 1) && (p == 1)) {
	  VectorXd df(n);
	  inc_gamma_deriv(x.row(0),f,df);
	  y(0) = f;      
	  y(1) = df.dot(x.row(1));    
	}
	
	if ((q <= 2) && (p == 2)) {
	  VectorXd df(n);
	  MatrixXd hess(n,n);
	  inc_gamma_deriv(x.row(0),f,df,hess);
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
	
	const Map<const MatrixXd> x = MatrixXd::Map(&(tx[0]),p+1,n);
	const Map<const VectorXd> py = VectorXd::Map(&(py_[0]),p+1);
	Map<MatrixXd> px = MatrixXd::Map(&(px_[0]),p+1,n);
	
	double f;
	VectorXd df(n);     
	MatrixXd dy(n,p+1);

	bool ok = (p <= 2);
	
	if (p == 0) {
	  inc_gamma_deriv(x.row(0), f, df);
	  dy.col(0) = df;
	}
	if (p >= 1) {
	  MatrixXd hess(n,n);
	  inc_gamma_deriv(x.row(0), f, df, hess);
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
    
    AScalar inc_gamma(const AScalar& z, const AScalar& r) {
      
      static atomic_inc_gamma func("atomic_inc_gamma");
      
      VectorXA x(2);
      VectorXA y(1);    
      x(0)=z;
      x(1)=r;
      func(x,y);
      return(y[0]);
    }    
  } // end namespace atomic
} // end namespace CppAD

#endif
