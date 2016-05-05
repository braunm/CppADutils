context("Wishart")

test_that("Wishart", {

    log_mvg <- function(v,p) {
        res <- p*(p-1)*log(pi)/4
        for (j in 1:p) {
            res <- res + lgamma(v+(1-j)/2)
        }
        return(res)
    }

    dwish <- function(X, v, S) {
        k <- NROW(X)
        res <- (v-k-1)*log(det(X))/2
        res <- res - sum(diag(solve(S,X)))/2
        res <- res - v*k*log(2)/2 - v*log(det(S))/2
        res <- res - log_mvg(v/2, k)
        return(res)
    }

    diwish <- function(X, v, S) {
        k <- NROW(X)
        res <- v*log(det(S))/2
        res <- res - sum(diag(S %*% solve(X)))/2
        res <- res - v*k*log(2)/2 - (v+k+1)*log(det(X))/2
        res <- res - log_mvg(v/2, k)
        return(res)
    }

    set.seed(123)

    k <- 4
    v <- k+4
    S <- rWishart(1, v, diag(k))[,,1]
    X <- rWishart(1, v, S)[,,1]
    Sinv <- solve(S)

    d1 <- dwish(X, v, S)
    d2 <- diwish(X, v, S)

    m1 <- Wish_test(X, v, S)
    m2 <- Inv_Wish_test(X, v, S)

    expect_equal(d1, m1)
    expect_equal(d2, m2)

})
