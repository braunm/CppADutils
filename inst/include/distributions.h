#ifndef __CPPADUTILS_DISTRIBUTIONS__
#define __CPPADUTILS_DISTRIBUTIONS__

#include <utilfuncs.h>

//#define CPPAD_USE_CPLUSPLUS_2011 1

inline AScalar dnorm_log(const AScalar&, const AScalar&, const AScalar&);
inline AScalar dnormTrunc0_log(const AScalar&, const AScalar&, const AScalar&);
inline AScalar pnorm_log(const AScalar&, const AScalar&, const AScalar&);
inline AScalar pnorm(const AScalar&, const AScalar&, const AScalar&);
inline AScalar dt_log(const AScalar&, const AScalar&, const AScalar&);
inline AScalar dhalft_log(const AScalar&, const AScalar&, const AScalar&);
inline AScalar MB_erf_pos(const AScalar&);
inline AScalar MB_erf(const AScalar&);



AScalar dnorm_log(const AScalar& x,
		  const AScalar& m,
		  const AScalar& s) {

  AScalar res = -M_LN_SQRT_2PI -log(s) -(x-m)*(x-m)/(2*s*s);
  return(res);
}

AScalar pnorm_log(const AScalar& x,
		  const AScalar& m,
		  const AScalar& s) {

  AScalar z = (x-m)/s;
  AScalar res = log1p(erf(z*M_SQRT1_2)) - M_LN2;

  return(res);
}

AScalar pnorm(const AScalar& x,
	      const AScalar& m,
	      const AScalar& s) {
  
  AScalar z = (x-m)/s;
  AScalar res = 0.5*(1+erf(z*M_SQRT1_2));
  
  return(res);
}

AScalar dnormTrunc0_log(const AScalar& x,
		  const AScalar& m,
		  const AScalar& s) {
  
  AScalar R = erf(-M_SQRT1_2 * m/s);
  AScalar res = dnorm_log(x, m, s) - log(1.0-R) + M_LN2;

  return(res);
}

AScalar dt_log(const AScalar& z, const AScalar& v, const AScalar& s) {

  AScalar res = lgamma(0.5*(v+1)) - lgamma(0.5*v);
  res -= 0.5*((v+1)*log1p(z*z/(s*s*v)) + log(v));
  res -= log(s) + M_LN_SQRT_PI;
  return(res);	     
}


AScalar dhalft_log(const AScalar& z, const AScalar& v, const AScalar& s) {
  AScalar res = M_LN2 + dt_log(z, v, s);
  return(res);
}

AScalar MB_erf_pos(const AScalar& x) {

  const AScalar p = 0.3275911;
  const AScalar a1 = 0.254829592;
  const AScalar a2 = -0.284496736;
  const AScalar a3 = 1.421413741;
  const AScalar a4 = -1.453152027;
  const AScalar a5 = 1.061405429;
  
  const AScalar t = 1/(1+p*x);

  AScalar res = t*(a1 + t*(a2 + t*(a3 + t*(a4 + t*a5))));
  res = 1 - res * exp(-x*x);
  return(res);
}


AScalar MB_erf(const AScalar& x) {

  AScalar r1 = MB_erf_pos(x);
  AScalar r2 = -MB_erf_pos(-x);
  AScalar res = CppAD::CondExpGe(x, AScalar(0), r1, r2);

  return(res);
}

#endif
