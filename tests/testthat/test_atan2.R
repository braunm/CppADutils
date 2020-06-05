context("atan2a")

test_that("atan2a",{
    require(numDeriv)
    fn <- function(x,y) {
        res <- atan2(x,y)
        return(res)
    }

    R_func <- function(x) {
        n <- length(x)/2
        res <- 0
        for (i in 1:n) {
            res <- res + fn(x[2*(i-1)+1], x[2*i])
        }
        return(res^2)
    }

    x <- c(.2,.3,.4,.5)
    R_val <- R_func(x)
    R_grad <- grad(R_func, x, method.args=list(r=8))
    R_hess <- hessian(R_func, x, method.args=list(r=8))
    R_hess_spLT <- tril(drop0(R_hess, 1e-8))

    c_list <- cppad_atan2a(x)
    c_val <- c_list$val
    c_grad <- c_list$grad
    c_hess_dense <- c_list$hess.dense
    c_hess_sp <- c_list$hess.sp
    c_hess_spLT <- tril(drop0(c_list$hess.spLT, 1e-8))

    expect_equal(c_val, R_val)
    expect_equal(c_grad, R_grad, tolerance=1e-7)
    expect_equal(drop0(c_hess_dense), drop0(R_hess, 1e-8), tolerance=1e-8)
    expect_equal(c_hess_dense, c_hess_sp)
    expect_equal(c_hess_spLT, tril(drop0(c_hess_sp, 1e-8)))
    expect_equal(c_hess_spLT, R_hess_spLT, tolerance=1e-8)
})



context("atan2b")

test_that("atan2b",{
    require(numDeriv)
##    require(CppADutils)
    fn1 <- function(x,y) {
        res1 <- atan2(x,y)
        return(res1)
    }


    fn2 <- function(a, b) {
        res2 <- a^b
        return(res2)
    }


    R_func <- function(x) {
        n <- length(x)/2
        res1 <- 0
        res2 <- 0
        for (i in 1:n) {
            res1 <- res1 + fn1(x[2*(i-1)+1], x[2*i])
            res2 <- res2 + fn2(x[2*(i-1)+1], x[2*i])
        }
        return(res1*res2)
    }

    x <- c(.2, .3, .4, .5)
    R_val <- R_func(x)
    R_grad <- grad(R_func, x, method.args=list(r=8))
    R_hess <- hessian(R_func, x, method.args=list(r=8))
    R_hess_spLT <- tril(drop0(R_hess, 1e-8))

    c_list <- cppad_atan2b(x)
    c_val <- c_list$val
    c_grad <- c_list$grad
    c_hess_dense <- c_list$hess.dense
    c_hess_sp <- c_list$hess.sp
    c_hess_spLT <- tril(drop0(c_list$hess.spLT, 1e-8))

    expect_equal(c_val, R_val)
    expect_equal(c_grad, R_grad, tolerance=1e-7)
    expect_equal(drop0(c_hess_dense), drop0(R_hess, 1e-8), tolerance=1e-8)
    expect_equal(c_hess_dense, c_hess_sp)
    expect_equal(c_hess_spLT, tril(drop0(c_hess_sp, 1e-8)))
    expect_equal(c_hess_spLT, R_hess_spLT, tolerance=1e-8)
})

