library(rstan)
library(plyr)
library(dplyr)


K <- 8
eta <- 3

iter <- 3
warmup <- 1
draws <- iter-warmup

n <- K*(K-1)/2
Y <- rnorm(n)

stanmath <- chol_corr_matrix(Y,K)
LS <- stanmath$S
jS <- stanmath$lp

mine <- LKJ_unwrap(Y,K)
LM <- mine$chol_corr
jM <- mine$log_jac

cn <- LKJ_const(eta, K)
cn_stan <- LKJ_stan_const(eta, K)

pdf1 <- LKJ(Y, eta, K)
pdf2 <- LKJ_chol(LM, eta) + jM
pdf3 <- LKJ_chol_stan(LM, eta) + jS




