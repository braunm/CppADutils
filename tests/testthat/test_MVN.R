context("MVN")

test_that("MVN", {

    require(mvtnorm)
    require(Matrix)
    set.seed(123)
    
    N <- 8
    k1 <- 2
    k2 <- 4
    k <- k1*k2
    S <- Matrix(0, k, k)
    
    for (i in 1:k2) {
        r <- ((i-1)*k1+1):(i*k1)
        S[r, r] <- rWishart(1, k1+4, diag(k1))[,,1]
    }
    
    S <- as(S,"dgCMatrix")
    D <- as(S, "matrix")
    
    Dinv <- solve(D)
    mu1 <- rnorm(k)
    X <- rmvnorm(N, mu1, D)
    mu2 <- rmvnorm(N, mu1, 2*diag(k))
    
    d1 <- dmvnorm(X, mu1, D, log=TRUE)
    d2 <- dmvnorm(X - mu2, rep(0,k), D, log=TRUE)
    
    m1 <- MVN_test(t(X), as.matrix(mu1), D, FALSE)
    m2 <- MVN_test(t(X), t(mu2), D, FALSE)
    
    e1 <- dmvnorm(X, mu1, Dinv, log=TRUE)
    e2 <- dmvnorm(X - mu2, rep(0,k), Dinv, log=TRUE)
    
    f1 <- MVN_test(t(X), as.matrix(mu1), D, TRUE)
    f2 <- MVN_test(t(X), t(mu2), D, TRUE)
    
    p1 <- Sparse_MVN_test(t(X), as.matrix(mu1), S, FALSE)
    p2 <- Sparse_MVN_test(t(X), t(mu2), S, FALSE)
    
    z1 <- Sparse_MVN_test(t(X), as.matrix(mu1), S, TRUE)
    z2 <- Sparse_MVN_test(t(X), t(mu2), S, TRUE)
    
    expect_equal(d1, m1)
    expect_equal(d2, m2)
    expect_equal(e1, f1)
    expect_equal(e2, f2)
    expect_equal(d1, p1)
    expect_equal(d2, p2)
    expect_equal(e1, z1)
    expect_equal(e2, z2)  
})
