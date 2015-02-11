
inline AScalar dnorm_log(const AScalar&, const AScalar&, const AScalar&);
inline AScalar dnormTrunc0_log(const AScalar&, const AScalar&, const AScalar&);
inline AScalar pnorm_log(const AScalar&, const AScalar&, const AScalar&);
inline AScalar dt_log(const AScalar&, const AScalar&, const AScalar&);
inline AScalar dhalft_log(const AScalar&, const AScalar&, const AScalar&);


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
  AScalar res = log(1+erf(z*M_SQRT1_2)) - M_LN2;

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
