#ifndef __CPPADUTILS_UTILFUNCS
#define __CPPADUTILS_UTILFUNCS

#include <cppad_atomics.h>

inline AScalar log_MVgamma(const AScalar&, const int&);
inline bool my_finite(const AScalar&);



AScalar log_MVgamma(const AScalar& v,
		    const int& k) {  
  AScalar res = 0.5 * k * (k-1) * M_LN_SQRT_PI;
  for (int j=1; j<=k; j++) {
    res += lgamma(v + 0.5*(1-j));
  }
  return res;
}

bool my_finite(const AScalar& x) {
  return( (abs(x) <= __DBL_MAX__ ) && ( x == x ) );
}

#endif
