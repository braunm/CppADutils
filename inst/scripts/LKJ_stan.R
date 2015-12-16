library(rstan)
library(plyr)
library(dplyr)


K <- 5
eta <- 3

iter <- 3
warmup <- 1
draws <- iter-warmup

n <- K*(K-1)/2
Y <- rnorm(n)
Z <- tanh(Y)

idx <- 1
jz1 <- 0
jz2 <- 0
ZM <- matrix(0,K,K)
for (j in 1:(K-1)) {
    for (i in (j+1):K) {
        idx <- (j-1)*K + i - j*(j+1)/2
        ZM[i,j] <- Z[idx]
        jz1 <- jz1 + (K-j-1)*log(1-Z[idx]^2)/2
        jz2 <- jz2 + log(1-Z[idx]^2)
    }
}


A <- corr_matrix(Y,K)
B <- chol_corr_matrix(Y,K)

LL <- LKJ_unwrap(Y,K)
L <- LL$chol_corr
log_jac <- LL$log_jac

cn <- LKJ_const(eta, K)
cn_stan <- LKJ_stan_const(eta, K)

lkj_chol<- LKJ_chol(L, eta)
lkj_chol_stan <- LKJ_chol_stan(L, eta)



S <- tcrossprod(L)
lkj <- cn + (eta-1)*log(det(S))
lkj_stan <- LKJ_stan(S, eta)


data <- list(L=LL$chol_corr,K=K,eta=eta)


stop()

cat("Compiling\n")
mod <- stan_model(file="inst/scripts/LKJ_stan.stan",
                  model_name="LKJ",
                  verbose=FALSE)

cat("Sampling\n")
fit <- sampling(mod,
                data=data,
                chains=1,
                iter=iter,
                warmup=warmup,
                thin=1,
                algorithm="NUTS",
                verbose=TRUE,
                refresh=1)



La <- extract(fit, permuted=TRUE)
lp <- La$lp__


