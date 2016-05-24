#include <Rcpp.h>
using Rcpp::NumericVector;
using Rcpp::NumericMatrix;






//' @title Inverse logit function
//' @description exp(x)/(1+exp(x))
//' @param X input numeric vector
//' @return exp(x)/(1+exp(x))
//' @export
//[[Rcpp::export]]
Rcpp::NumericVector invlogit(const NumericVector& X) {

  int k = X.size();
  NumericVector Y(k);
  for (int i=0; i<k; i++) {
    Y(i) = exp(X(i) - R::log1pexp(X(i)));
  }
  return(Y);
}

//' @title log(1+x)-x
//' @description log(1+x)-x, accurate even for small |x|
//' @param X input numeric vector
//' @return log(1+x)-x
//' @export
//[[Rcpp::export]]
Rcpp::NumericVector log1pmx(const NumericVector& X) {

  int k = X.size();
  NumericVector Y(k);
  for (int i=0; i<k; i++) {
    Y(i) = R::log1pmx(X(i));
  }
  return(Y);
}

//' @title Log inverse logit function
//' @description Log inverse logit, accurate even for very negative x
//' @param X input numeric vector
//' @return p = log(exp(x)/(1+exp(x)) = x-log(1+exp(x))
//' @export
//[[Rcpp::export]]
Rcpp::NumericVector loginvlogit(const NumericVector& X) {

  int k = X.size();
  NumericVector Y(k);
  for (int i=0; i<k; i++) {
    Y(i) = X(i) - R::log1pexp(X(i));
  }
  return(Y);
}

//' @title log(1+exp(x))
//' @description log(1+exp(x)), accurate even for very large X.
//' @param X input numeric vector
//' @return log(1+exp(x))
//' @export
//[[Rcpp::export]]
Rcpp::NumericVector log1pexp(const NumericVector& X) {

  int k = X.size();
  NumericVector Y(k);
  for (int i=0; i<k; i++) {
    Y(i) = R::log1pexp(X(i));
  }
  return(Y);
}

//' @title log(gamma(1+x))
//' @description log(gamma(1+x)), accurate even for very small |x|
//' @param X input numeric vector
//' @return log(gamma(1+x))
//' @export
//[[Rcpp::export]]
Rcpp::NumericVector lgamma1p(const NumericVector& X) {
  int k = X.size();
  NumericVector Y(k);
  for (int i=0; i<k; i++) {
    Y(i) = R::lgamma1p(X(i));
  }
  return(Y);
}


//' @title log(gamma(exp(x)))
//' @description log(gamma(exp(x))), accurate even for very small |x|
//' @param X input numeric vector
//' @return log(gamma(exp(x)))
//' @export
//[[Rcpp::export]]
Rcpp::NumericVector lgammaexp(const NumericVector& X) {
  int k = X.size();
  NumericVector Y(k);
  for (int i=0; i<k; i++) {
    Y(i) = R::lgamma1p(exp(X(i))) - X(i);
  }
  return(Y);
}


Rcpp::NumericVector lgamma(const NumericVector& X) {
  int k = X.size();
  NumericVector Y(k);
  for (int i=0; i<k; i++) {
    Y(i) = R::lgammafn(X(i));
  }
  return(Y);
}

Rcpp::NumericVector gamma(const NumericVector& X) {

  int k = X.size();
  NumericVector Y(k);
  for (int i=0; i<k; i++) {
    Y(i) = R::gammafn(X(i));
  }
  return(Y);
}


