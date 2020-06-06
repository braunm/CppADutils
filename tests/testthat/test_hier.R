context("AD for hierarchical model")

test_that("hier1",{

    require(RcppEigen)
    require(numDeriv)
    require(mvtnorm)

    set.seed(123)
    data(ADtest)
    k <- NCOL(ADtest$data[[1]]$X)
    beta <- seq(1,k, length=k)
    sigma_b <- ADtest$priors$sigma.b
    B <- t(rmvnorm(length(ADtest$data), beta, sigma_b))
    V <- c(as.vector(B), beta)

    get_LLi <- function(D, b, s) {
        m <- D$X %*% b
        res <- sum(dnorm(D$Y, m, s, log=TRUE))
        return(res)
    }

    make_f <- function(D) {
        sigma_b <- D$priors$sigma.b
        m_beta <- D$priors$mean.beta
        sigma_beta <- D$priors$sigma.beta
        s <- D$priors$sdY
        N <- length(D$data)
        function(V) {
            LL <- 0
            pr <- 0
            b <- matrix(V[1:(N*k)], k, N)
            beta <- V[(N*k+1):length(V)]
            for (i in 1:N) {
                LL <- LL + get_LLi(D$data[[i]], b[,i], s)
                pr <- pr + dmvnorm(b[,i], beta, sigma_b, log=TRUE)
            }
            hyp <- dmvnorm(beta, m_beta, sigma_beta, log=TRUE)
            res <- LL + pr + hyp
            return(res)
        }
    }

    get_f <- make_f(ADtest)

    f1 <- get_f(V)
    g1 <- grad(get_f, V)
    h1 <- hessian(get_f, V)
    s1 <- drop0(h1, 1e-8)

    ## From C++ function using AD
    cl <- new("adtest", ADtest)
    cl$record.tape(V)
    cl$init.hessian(V)
    f2 <- cl$get.f(V)
    g2 <- cl$get.df(V)
    s2 <- cl$get.hessian.sparse(V)
    h2 <- drop0(cl$get.hessian(V))

    expect_equal(f1, f2)
    expect_equal(g1, g2)
    expect_equal(s1, h2)
    expect_equal(s1, s2)
})





