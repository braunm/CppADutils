#ifndef __CPPADUTILS_ATOMICS
#define __CPPADUTILS_ATOMICS

#include <cppad_atomic/expm1_at.h>
//#include <cppad_atomic/dt_log_at.h>
#include <cppad_atomic/incgamma_at.h>
#include <cppad_atomic/incbeta_at.h>
#include <cppad_atomic/invlogit_at.h>
#include <cppad_atomic/lbeta_at.h>
#include <cppad_atomic/lgamma_at.h>
#include <cppad_atomic/lgamma1p_at.h>
#include <cppad_atomic/lgammaexp_at.h>
#include <cppad_atomic/log1p_at.h>
#include <cppad_atomic/log1pexp_at.h>
#include <cppad_atomic/log1pmx_at.h>
#include <cppad_atomic/loginvlogit_at.h>
#include <cppad_atomic/dbeta_log_at.h>
#include <cppad_atomic/dlogitbeta_log_at.h>

inline AScalar expm1(const AScalar&);
inline AScalar dbeta_log(const AScalar&, const AScalar&, const AScalar&);
inline AScalar dlogitbeta_log(const AScalar&, const AScalar&, const AScalar&);
//inline AScalar dhalft_log(const AScalar&, const AScalar&, const AScalar&);
//inline AScalar dt_log(const AScalar&, const AScalar&, const AScalar&);
inline AScalar incbeta(const AScalar&, const AScalar&, const AScalar&);
inline AScalar incgamma(const AScalar&, const AScalar&);
inline AScalar invlogit(const AScalar&);
inline AScalar lbeta(const AScalar&, const AScalar&);
inline AScalar lgamma(const AScalar&);
inline AScalar lgamma1p(const AScalar&);
inline AScalar lgammaexp(const AScalar&);
inline AScalar log1p(const AScalar&);
inline AScalar log1pexp(const AScalar&);
inline AScalar log1pmx(const AScalar&);
inline AScalar loginvlogit(const AScalar&);


#include <cppad_atomic/distributions.h>

AScalar lgamma(const AScalar& a) {
  
  static mb_atomic<lgamma_cl> lgamma_func("atomic_lgamma");
  VectorXA y(1);  
  VectorXA x(1);
  x << a;
  lgamma_func(x,y);
  return(y[0]);
}


AScalar dbeta_log(const AScalar& p, const AScalar& a, const AScalar& b) {
      
  static mb_atomic<dbeta_log_cl> dbeta_log_func("atomic_dbeta_log");
  VectorXA y(1);     
  VectorXA x(3);
  x << p, a, b; 
  dbeta_log_func(x,y);
  return(y[0]);
}

AScalar dlogitbeta_log(const AScalar& z, const AScalar& a, const AScalar& b) {
      
  static mb_atomic<dlogitbeta_log_cl> dlogitbeta_log_func("atomic_dlogitbeta_log");
  VectorXA y(1);     
  VectorXA x(3);
  x << z, a, b;
  dlogitbeta_log_func(x,y);
  return(y[0]);
}

/* AScalar dt_log(const AScalar& z, const AScalar& v, const AScalar& s) { */
      
/*   static mb_atomic<dt_log_cl> dt_log_func("atomic_dt_log"); */
/*   VectorXA y(1);      */
/*   VectorXA x(3); */
/*   x << z, v, s;  */
/*   dt_log_func(x,y); */
/*   dt_log_func.clear(); */
/*   return(y[0]); */
/* } */

/* AScalar dhalft_log(const AScalar& z, const AScalar& v, const AScalar& s) { */
      
/*   AScalar res = M_LN2 + dt_log(z, v, s); */
/*   return(res); */
/* } */

AScalar expm1(const AScalar& a) {
  
  static mb_atomic<expm1_cl> expm1_func("atomic_expm1");
  VectorXA y(1);
  VectorXA x(1);
  x << a;
  expm1_func(x,y);
  return(y[0]);
}

AScalar incbeta(const AScalar& z, const AScalar& a, const AScalar& b) {
      
  static mb_atomic<incbeta_cl> incbeta_func("atomic_incbeta");
  VectorXA y(1);     
  VectorXA x(3);
  x << z, a, b; 
  incbeta_func(x,y);
  incbeta_func.clear();
  return(y[0]);
}

AScalar incgamma(const AScalar& z, const AScalar& r) {
  
  static mb_atomic<incgamma_cl> incgamma_func("atomic_incgamma");
  VectorXA y(1);    
  VectorXA x(2);
  x << z, r;
  incgamma_func(x,y);
  incgamma_func.clear();
  return(y[0]);
}

AScalar invlogit(const AScalar& a) {
  
  static mb_atomic<invlogit_cl> invlogit_func("atomic_invlogit");
  VectorXA y(1);
  VectorXA x(1);
  x << a;
  invlogit_func(x,y);
  return(y[0]);
}

AScalar lbeta (const AScalar& a, const AScalar& b) {
  
  static mb_atomic<lbeta_cl> lbeta_func("atomic_lbeta"); 
  VectorXA y(1);
  VectorXA x(2);
  x << a, b;
  lbeta_func(x,y);
  lbeta_func.clear();
  return(y[0]);
}

AScalar lgamma1p(const AScalar& a) {
  static mb_atomic<lgamma1p_cl> lgamma1p_func("atomic_lgamma1p");
  VectorXA x(1);
  VectorXA y(1);
  x << a;
  lgamma1p_func(x,y);
  return(y[0]);
}

AScalar lgammaexp(const AScalar& a) {
  
  static mb_atomic<lgammaexp_cl> lgammaexp_func("atomic_lgammaexp");
  VectorXA y(1);
  VectorXA x(1);
  x << a;
  lgammaexp_func(x,y);
  return(y[0]);
}

AScalar log1p(const AScalar& a) {
  
  static mb_atomic<log1p_cl> log1p_func("atomic_log1p");
  VectorXA y(1);
  VectorXA x(1);
  x << a;
  log1p_func(x,y);
  return(y[0]);
}

AScalar log1pexp(const AScalar& a) {
  
  static mb_atomic<log1pexp_cl> log1pexp_func("atomic_log1pexp");
  VectorXA y(1);
  VectorXA x(1);
  x << a;
  log1pexp_func(x,y);
  return(y[0]);
}

AScalar log1pmx(const AScalar& a) {
  
  static mb_atomic<log1pmx_cl> log1pmx_func("atomic_log1pmx");
  VectorXA y(1);
  VectorXA x(1);
  x << a;
  log1pmx_func(x,y);
  return(y[0]);
}

AScalar loginvlogit(const AScalar& a) {
  
  static mb_atomic<loginvlogit_cl> loginvlogit_func("atomic_loginvlogit");
  VectorXA y(1);
  VectorXA x(1);
  x << a;
  loginvlogit_func(x,y);
  return(y[0]);
}
    


#endif
