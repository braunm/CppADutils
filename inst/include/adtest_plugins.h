double get_f_direct(const NumericVector& P_) {

  if (!tape_ready)
    throw MyException("tape not ready",__FILE__,__LINE__);

  
  VectorXA P = VectorXd::Map(P_.begin(), nvars()).template cast<AScalar>();
  AScalar fr = model->eval_f(P); // inherited from base class
  double res = Value(fr);
  return(res);
}

double get_LL(const NumericVector& P_) {
  if (!tape_ready)
    throw MyException("tape not ready",__FILE__,__LINE__);

  VectorXA P = VectorXd::Map(P_.begin(), nvars()).template cast<AScalar>();  
  AScalar fr = model->eval_LL(P); // inherited from base class
  double res = Value(fr);
  return(res);
}

NumericVector get_LLi(const NumericVector& P_) {
   if (!tape_ready)
    throw MyException("tape not ready",__FILE__,__LINE__);

   const int& N = model->N;
   VectorXA P = VectorXd::Map(P_.begin(), P_.size()).template cast<AScalar>();
   VectorXA LLi(N);
   model->eval_LLi(P, LLi);
   NumericVector res(N);
   for (size_t i=0; i<N; i++) {
     res(i) = CppAD::Value(LLi(i));
   }
   return(res);
}

NumericVector get_prior_i(const NumericVector& P_) {
   if (!tape_ready)
    throw MyException("tape not ready",__FILE__,__LINE__);

   const int& N = model->N;
   VectorXA P = VectorXd::Map(P_.begin(), P_.size()).template cast<AScalar>();
   VectorXA prior_i(N);
   model->eval_prior_i(P, prior_i);
   NumericVector res(N);
   for (size_t i=0; i<N; i++) {
     res(i) = Value(prior_i(i));
   }
   return(res);
}

