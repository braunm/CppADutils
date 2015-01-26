#ifndef __adtest
#define __adtest

#include <mb_base.h>
#include <except.h>
#include <utilfuncs.h>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Cholesky>
#include <cppad_atomic/atomics.h>
#include <MVN_AD.h>
#include <wish_AD.h>

using Eigen::MatrixBase;

class adtest  {

  typedef CppAD::AD<double> AScalar;
  typedef Eigen::Matrix<AScalar, Eigen::Dynamic, Eigen::Dynamic> MatrixXA;
  typedef Eigen::Matrix<AScalar, Eigen::Dynamic, 1> VectorXA;

 public:

  // constructor and eval_f are required by MB_Base
  
  adtest(const List&);

  template<typename Tpars>
    AScalar eval_f(const MatrixBase<Tpars>&); // required


  /* Called by optional plugins. The wrappers need to be members of
     MB_Base for the RCPP_MODULE mechanism to see them (I think).  
   */

  
  template<typename Tpars>
    AScalar eval_LL(const MatrixBase<Tpars>&);
  
  template<typename Tpars, typename TL>
    void eval_LLi(const MatrixBase<Tpars>&, const MatrixBase<TL>&);

  template<typename Tpars, typename TL>
    void eval_prior_i(const MatrixBase<Tpars>&, const MatrixBase<TL>&);

  template<typename Tpars>
    AScalar eval_prior(const MatrixBase<Tpars>&);  
 
  template<typename Tpars>
    AScalar eval_hyperprior(const MatrixBase<Tpars>&);

  int N, k, T; //must be public to use in plugins
 
 private:

 template<typename Tpars>
   void unwrap_params(const MatrixBase<Tpars>&);
  
  AScalar eval_LL();
  AScalar eval_prior();
  AScalar eval_hyperprior();
 
  template<typename TL>
    void eval_LLi(const MatrixBase<TL>&);

  template<typename TL>
    void eval_prior_i(const MatrixBase<TL>&);

 protected:
  
  // data and priors.

  std::vector<MatrixXA> X;
  std::vector<VectorXA> Y;

  MatrixXA Sb; // cov of beta
  VectorXA Mbeta;  // mean of b  
  MatrixXA Sbeta; // cov of b
  AScalar sdY;
  
  MatrixXA CholSbeta;
  AScalar logDetCholSbeta;
  MatrixXA CholSb;
  AScalar logDetCholSb;

  // Parameters

  MatrixXA B;
  VectorXA beta;

  
};

adtest::adtest(const List& params)
{

  using Rcpp::as;
  using Rcpp::List;
  using Eigen::MatrixXd;
  using Eigen::VectorXd;
  
  List & pars = static_cast<List &>(const_cast<List &>(params));

  Rcpp::ListOf<List> data = as<Rcpp::ListOf<List> >(pars["data"]);

  N = data.size();
  X.resize(N);
  Y.resize(N);

  for (size_t i=0; i<N; i++) {
    X[i] = as<MatrixXd>(data[i]["X"]).cast<AScalar>();
    Y[i] = as<MatrixXd>(data[i]["Y"]).cast<AScalar>();
  }
  k = X[1].cols();
  T = Y[1].size();
  
  List priors = as<List>(pars["priors"]);
  Sb = as<MatrixXd>(priors["sigma.b"]).cast<AScalar>();
  Mbeta = as<VectorXd>(priors["mean.beta"]).cast<AScalar>();
  Sbeta = as<MatrixXd>(priors["sigma.beta"]).cast<AScalar>();
  sdY = as<double>(priors["sdY"]);
  Eigen::LLT<MatrixXA> tmp(Sb);
  CholSb = tmp.matrixL();
  tmp.compute(Sbeta);
  CholSbeta = tmp.matrixL();
  
  logDetCholSb = CholSb.diagonal().array().log().sum();
  logDetCholSbeta = CholSbeta.diagonal().array().log().sum();

  // reserve space for parameters
  B.resize(k, N);
  beta.resize(k);
}

template<typename Tpars>
void adtest::unwrap_params(const MatrixBase<Tpars>& P)
{
  B = MatrixXA::Map(P.derived().data(), k, N);
  beta = P.derived().segment(k*N, k);
}

template<typename Tpars>
AScalar adtest::eval_f(const MatrixBase<Tpars>& P) {
  
  unwrap_params(P);
   
  AScalar LL = eval_LL();
  AScalar prior = eval_prior();
  AScalar hyperprior = eval_hyperprior();  
  AScalar f = LL + prior + hyperprior;
    
  return(f);
}

template<typename TL>
void adtest::eval_LLi(const MatrixBase<TL>& LLi_)
{ 

  MatrixBase<TL>& LLi = const_cast<MatrixBase<TL>& >(LLi_);
  LLi.setZero();
  
  for (size_t i=0; i<N; i++) {
    VectorXA mu = X[i] * B.col(i);
    for (size_t j=0; j<T; j++) {
      LLi(i) += dnorm_log(Y[i](j), mu(j), sdY);
    }
  }
}



AScalar adtest::eval_LL()
{ 
  VectorXA LLi(N); //vector to store individual-level likelihoods
  eval_LLi(LLi);
  AScalar LL = LLi.sum();
  return(LL);

}


template<typename TL>
void adtest::eval_prior_i(const MatrixBase<TL>& prior_i_) {
  
  MatrixBase<TL>& prior_i = const_cast<MatrixBase<TL>& >(prior_i_);

  // MVN prior on B

  prior_i.setConstant(-logDetCholSb - k * M_LN_SQRT_2PI);
  MatrixXA Z = CholSb.triangularView<Eigen::Lower>().solve(B.colwise() - beta);
  prior_i -= 0.5*((Z.array()*Z.array()).colwise().sum()).matrix().transpose();
}

AScalar adtest::eval_prior() {
  
  VectorXA log_prior_i(N);
  eval_prior_i(log_prior_i);

  AScalar log_prior = log_prior_i.sum();
  assert(my_finite(log_prior));
  return(log_prior);
}




AScalar adtest::eval_hyperprior() {

  VectorXA z = CholSbeta.triangularView<Eigen::Lower>().solve(beta - Mbeta);
  AScalar res = -k * M_LN_SQRT_2PI - logDetCholSbeta -0.5*z.squaredNorm();
  return(res);
}


template<typename Tpars>
AScalar adtest::eval_LL(const MatrixBase<Tpars>& P) {
  
  unwrap_params(P);
  typename Tpars::Scalar LL = eval_LL();
  
  return(LL);
 
}



template<typename Tpars>
AScalar adtest::eval_prior(const MatrixBase<Tpars>& P) {
  
  unwrap_params(P);
  typename Tpars::Scalar prior = eval_prior();
  
  return(prior);
 
}

template<typename Tpars>
AScalar adtest::eval_hyperprior(const MatrixBase<Tpars>& P) {
  
  unwrap_params(P);
  typename Tpars::Scalar hyperprior = eval_hyperprior();
  
  return(hyperprior);
 
}

template<typename Tpars, typename TL>
void adtest::eval_LLi(const MatrixBase<Tpars>& P,
		      const MatrixBase<TL>& LLi_) {
  
  MatrixBase<TL>& LLi = const_cast<MatrixBase<TL>& >(LLi_); 
  unwrap_params(P);
  eval_LLi(LLi);
}

template<typename Tpars, typename TL>
void adtest::eval_prior_i(const MatrixBase<Tpars>& P,
						const MatrixBase<TL>& log_prior_i_) {
 
  MatrixBase<TL>& log_prior_i = const_cast<MatrixBase<TL>& >(log_prior_i_); 
  unwrap_params(P);
  eval_prior_i(log_prior_i);

}

#endif
