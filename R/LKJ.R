## Test functions for LKJ prior


#' @title LKJ
#' @param Y vector of unconstrained parameters
#' @param eta LKJ parameter
#' @param K dimension (integer)
#' @return log LKJ pdf of Cholesky decomp, given unconstrained vector
lkj_R <- function(Y, eta, K) {
    W <- lkj_unwrap_R(Y,K)
    L <- W$L
    logjac <- W$logjac
    lkj_chol <- lkj_chol_logpdf_R(L, eta)
    res <- lkj_chol + logjac
    return(res)
}

#' @rdname LKJ
lkj_const_R <- function(eta, K) {
    res <- 0
    s <- 0
    for (i in 1:(K-1)) {
        diff <- K-i
        barg <- eta + 0.5*(diff-1)
        res <- res + diff*(2*lgamma(barg)-lgamma(2*barg))
        s <- s + diff*(2*eta-2+diff)
    }
    res <- res + s*log(2)
    return(res)
}

#' @rdname LKJ
lkj_unwrap_R <- function(Y, K) {

    Z <- matrix(0,K,K)
    W <- matrix(0,K,K)
    idx <- 1
    logjac <- 0
    W[1,1] <- 1
    for (i in 2:K) {
        sum_sqs <- 0
        for (j in 1:(i-1)) {
            Z[i,j] <- tanh(Y[idx])
            idx <- idx+1
            logjac <- logjac + log(1-Z[i,j]^2) + 0.5*log(1-sum_sqs)
            W[i,j] <- Z[i,j]*sqrt(1-sum_sqs)
            sum_sqs <- sum_sqs+W[i,j]^2
        }
        W[i,i] <- sqrt(1-sum_sqs)
    }
    res <- list(L=W, logjac=logjac)
    return(res)
}

#' @rdname LKJ
#' @param L lower triangular Cholesky decomposition
#' @return log LKJ pdf of Cholesky decomp, given unconstrained vector
lkj_chol_logpdf_R <- function(L, eta) {
    K <- NROW(L)
    q <- (K-4+2*eta)-seq(0,K-2)
    c <- lkj_const_R(eta, K)
    log_diags <- log(diag(L)[2:K])
    v <- q*log_diags
    res <- c + sum(v)
    return(res)
}


