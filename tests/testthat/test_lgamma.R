context("lgamma1")

test_that("lgamma1",{
    require(numDeriv)
    R_func <- function(x) {
        n <- length(x)/2
        res <- 0
        for (i in 1:n) {
            res <- res + lgamma(x[2*(i-1)+1]) * lgamma(x[2*i])
        }
        return(res)
    }

    x <- c(.4, 2.2, 4.4, 3.3, 5, 1.5)

    R_val <- R_func(x)
    R_grad <- grad(R_func, x)
    R_hess <- hessian(R_func, x, method.args=list(r=6))
    R_hess_spLT <- tril(drop0(R_hess, 1e-9))
    
    c_list <- cppad_lgamma1(x)
    c_val <- c_list$val
    c_grad <- c_list$grad
    c_hess_dense <- c_list$hess.dense
    c_hess_sp <- c_list$hess.sp
    c_hess_spLT <- tril(c_list$hess.spLT)
 
    expect_equal(c_val, R_val)
    expect_equal(c_grad, R_grad)
    expect_equal(c_hess_dense, R_hess)
    expect_equal(c_hess_dense, c_hess_sp)
    expect_equal(c_hess_spLT, tril(drop0(c_hess_sp, 1e-9)))
    expect_equal(c_hess_spLT, R_hess_spLT)   
})

context("lgammaLogExp")

test_that("lgammaLogExp",{

    R_func <- function(x) {
        n <- length(x)/2
        res <- 0
        for (i in 1:n) {
            a <- x[2*(i-1)+1]
            b <- x[2*i]
            res <- res + (lgamma(log(a)) * lgamma(exp(b)))
        }
        return(res)
    }

    x <- c(1.1, 1.2, 1.3, 1.4, 1.5, 1.6)

    R_val <- R_func(x)
    R_grad <- grad(R_func, x)
    R_hess <- hessian(R_func, x, method.args=list(r=6))
    R_hess_spLT <- tril(drop0(R_hess, 1e-9))
    
    c_list <- cppad_lgammaLogExp(x)
    c_val <- c_list$val
    c_grad <- c_list$grad
    c_hess_dense <- c_list$hess.dense
    c_hess_sp <- c_list$hess.sp
    c_hess_spLT <- tril(c_list$hess.spLT,1e-9)

    expect_equal(c_val, R_val)
    expect_equal(c_grad, R_grad)
    expect_equal(c_hess_dense, R_hess)
    expect_equal(c_hess_dense, c_hess_sp)
    expect_equal(c_hess_spLT, tril(drop0(c_hess_sp, 1e-9)))
    expect_equal(c_hess_spLT, R_hess_spLT)   
})

    
