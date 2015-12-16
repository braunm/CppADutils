rm(list=ls())

library(rstan)
library(plyr)
library(dplyr)
set.seed(123)

## R functions for LKJ

lkj_const <- function(eta, K) {
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

lkj_unwrap <- function(Y, K) {

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

lkj_chol_logpdf <- function(L, eta) {
    k <- NROW(L)
    q <- (K-4+2*eta)-seq(0,K-2)
    c <- lkj_const(eta, K)
    log_diags <- log(diag(L)[2:K])
    v <- q*log_diags
    res <- c + sum(v)
    return(res)
}

LKJ_R <- function(Y, eta, K) {
    W <- lkj_unwrap(Y,K)
    L <- W$L
    logjac <- W$logjac
    lkj_chol <- lkj_chol_logpdf(L, eta)
    res <- lkj_chol + logjac
    return(res)
}




K <- 5
eta <- 1

Y <- rnorm(choose(K,2))

stanmath <- chol_corr_matrix(Y,K)
LS <- stanmath$S
jS <- stanmath$lp

mine <- LKJ_test(Y, eta, K)
LM <- mine$L
jM <- mine$logjac


cn <- LKJ_const(eta, K)
cn_stan <- LKJ_stan_const(eta, K)
cn_R <- lkj_const(eta, K)

pdf1 <- mine$logpdf
pdf2 <- LKJ_chol_stan(LM, eta) + jS
pdf3 <- LKJ_R(Y, eta, K)

print(pdf1)
print(pdf2)
print(pdf3)


