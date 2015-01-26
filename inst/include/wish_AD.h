#ifndef __AD_WISHART
#define __AD_WISHART

#include <utilfuncs.h>
#include <cppad_atomic/atomics.h>


using Eigen::MatrixBase;
using Eigen::Lower;
using Eigen::Dynamic;
using Eigen::Matrix;

template<typename TG, typename TM>
typename TG::Scalar Wishart_logpdf(const MatrixBase<TG>& chol_G,
				   const typename TG::Scalar& nu,
				   const MatrixBase<TM>& chol_S) {
  
  // S is the standard Wishart scale parameter, such that E(G) = nu*S


  typedef typename TG::Scalar Scalar;
  typedef Matrix<Scalar, Dynamic, Dynamic> MatrixXS;
  
  const int k = chol_G.rows();
  
  Scalar mvg=0;
  mvg = log_MVgamma(nu/2., k);

  Scalar c = nu * k * M_LN2 / 2. + mvg + nu*chol_S.diagonal().array().log().sum();
  
  Scalar d = (nu-k-1)*chol_G.diagonal().array().log().sum();
  
  
  MatrixXS Z = chol_S.template triangularView<Lower>().solve(chol_G);
  
  
  Scalar res = d - c - 0.5 * Z.squaredNorm();
  
  return(res);
}



template<typename TG, typename TM>
typename TG::Scalar Inv_Wishart_logpdf(const MatrixBase<TG>& chol_G,
				       const typename TG::Scalar& nu,
				       const MatrixBase<TM>& chol_S) {
  
  // S is the standard inverse Wishart scale parameter, such that E(G) = S/(nu-k-1)


  typedef typename TG::Scalar Scalar;
  typedef Matrix<Scalar, Dynamic, Dynamic> MatrixXS;
  
  const int k = chol_G.rows();

  Scalar mvg=0;
  mvg = log_MVgamma(nu/2., k);
  
  Scalar c = -nu * k * M_LN2/2. - mvg;

  Scalar d = -(nu+k+1)*chol_G.diagonal().array().log().sum() + nu*chol_S.diagonal().array().log().sum();

  MatrixXS Z = chol_G.template triangularView<Lower>().solve(chol_S);


  Scalar res = c + d - Z.squaredNorm()/2;

  return(res);
}



#endif

