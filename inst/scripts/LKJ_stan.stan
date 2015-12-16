data {
    int<lower=1> K;
    real<lower=0> eta;
    matrix[K,K] L;
}

parameters {
    real<lower=0,upper=1> u;
}

model {

    matrix[K,K] M;
    M <- tcrossprod(L);
    increment_log_prob(lkj_corr_log(M, eta));
}


