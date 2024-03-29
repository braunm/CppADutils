% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{MatNorm_test}
\alias{MatNorm_test}
\alias{MatNorm_sparse_test}
\title{MatNorm_test}
\usage{
MatNorm_test(X_, M_, chol_U_, chol_V_, isPrec)

MatNorm_sparse_test(X_, M_, U_, V_, isPrec)
}
\arguments{
\item{X_}{numeric matrix (k x N)}

\item{M_}{numeric matrix (k x N)}

\item{chol_U_}{k x k lower triangular matrix of lower triangle of chol(U)}

\item{chol_V_}{N x N lower triangular matrix of lower triangle of chol(V)}

\item{isPrec}{TRUE if U and V are precision matrices.  FALSE is U and V are covariance matrices}

\item{U_}{k x k dgCMatrix (full, not Cholesky)}

\item{V_}{N x N dgCMatrix (full, not Cholesky)}
}
\value{
log pdf of matrix normal distribution

log pdf of matrix normal distribution
}
\description{
MatNorm_test
}
