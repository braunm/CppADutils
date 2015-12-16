## test for matrix normal

context("MatNorm")

test_that("MatNorm", {

    library(Matrix)
    library(mvtnorm)

    MatNorm_R <- function(Y, M, chol_U, chol_V, isPrec) {

        S <- tcrossprod(kronecker(chol_V, chol_U))
        if (isPrec) {
            C <- solve(S)
        } else {
            C <- S
        }
        res <- dmvnorm(as.vector(Y), as.vector(M), C, log=TRUE)
        return(res)
    }


    k <- 4
    N <- 7

    X <- matrix(rnorm(k*N),k,N)
    M <- matrix(rnorm(k*N),k,N)

    U <- rWishart(1, k+4, 2*diag(k))[,,1]
    V <- rWishart(1, N+4, 3*diag(N))[,,1]
    U[k,1] <- 0.3
    U[1,k] <- U[k,1]
    V[N,1] <- -0.2
    V[1,N] <- V[N,1]

    chol_U <- t(chol(U))
    chol_V <- t(chol(V))

    US <- as(U, "dgCMatrix")
    VS <- as(V, "dgCMatrix")

    PR1 <- MatNorm_R(X, M, chol_U, chol_V, TRUE)
    PC1 <- MatNorm_test(X, M, chol_U, chol_V, TRUE)
    SC1 <- MatNorm_sparse_test(X, M, US, VS, TRUE)

    expect_equal(PR1, PC1)
    expect_equal(PR1, SC1)
    expect_equal(PC1, SC1)

    PR2 <- MatNorm_R(X, M, chol_U, chol_V, FALSE)
    PC2 <- MatNorm_test(X, M, chol_U, chol_V, FALSE)
    SC2 <- MatNorm_sparse_test(X, M, US, VS, FALSE)

    expect_equal(PR2, PC2)
    expect_equal(PR2, SC2)
    expect_equal(PC2, SC2)



})


