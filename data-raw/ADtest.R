if (!require(mvtnorm)) {
    stop("package mvtnorm required")
}
if (!require(devtools)) {
    stop("package devtools required")
}


set.seed(123)

k <- 2
T <- 50
N <- 5
mean.X <- 0
sd.X <- 1
sd.Y <- .2
beta <- seq(1,k, length=k)
sigma_b <- 2 * diag(k)
sigma_b[k,1] <- sigma_b[1,k] <- sigma_b[1,1]/5
B <- t(rmvnorm(N, beta, sigma_b))

m_beta<- rep(0,k)
sigma_beta <- diag(k)
sigma_beta[k,1] <- sigma_beta[1,k] <- sigma_beta[1,1]/10

data <- vector("list", length=N)

for (i in 1:N) {
    X <- matrix(rnorm(T*k, mean.X, sd.X), T, k)
    Y <- X %*% B[,i] + rnorm(T, 0, sd.Y)
    data[[i]] <- list(Y=Y, X=X)
}

ADtest <- list(data=data,
          priors=list(sigma.b=sigma_b,
              mean.beta=m_beta,
              sigma.beta=sigma_beta,
              sdY=sd.Y)
          )

devtools::use_data(ADtest)







