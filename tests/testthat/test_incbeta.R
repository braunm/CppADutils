context("incbeta")

test_that("incbeta",{
    require(numDeriv)
    R_func <- function(x) {
        n <- length(x)/3
        res <- 0
        for (i in 1:n) {
            res <- res + pbeta(x[3*(i-1)+1], x[3*(i-1)+2], x[3*i])
        }
        return(res*res)
    }

    x <- c(0.8, 2.3, 4.4, 0.3, 2.0, 2.1)

    R_val <- R_func(x)
    R_grad <- grad(R_func, x, method.args=list(r=8))
    R_hess <- hessian(R_func, x, method.args=list(r=8))
    R_hess_spLT <- tril(drop0(R_hess, 1e-8))
    
    c_list <- cppad_incbeta(x)
    c_val <- c_list$val
    c_grad <- c_list$grad
    c_hess_dense <- c_list$hess.dense
    c_hess_sp <- c_list$hess.sp
    c_hess_spLT <- tril(c_list$hess.spLT)
 
    expect_equal(c_val, R_val)
    expect_equal(c_grad, R_grad, tolerance=1e-5)
    expect_equal(c_hess_dense, R_hess, tolerance=1e-6)
    expect_equal(c_hess_dense, c_hess_sp)
    expect_equal(c_hess_spLT, tril(drop0(c_hess_sp, 1e-8)))
    expect_equal(c_hess_spLT, R_hess_spLT, tolerance=1e-6)   
})




context("incbeta2")
test_that("incbeta2",{
    require(numDeriv)

    halft <- function(z, v, s) {
        res <- lgamma((v+1)/2) - lgamma(v/2)
        res <- res - 0.5*((v+1)*log(1+z^2/(v*s^2)) + log(pi * v)) - log(s)
        res <- res + log(2)
        return(res)
    }
   
    R_func <- function(x) {
        n <- length(x)/3
        res1 <- 0
##        res2 <- 0
        for (i in 1:n) {
            res1 <- res1 + pbeta(x[3*(i-1)+1], x[3*(i-1)+2], x[3*i])
##            res2 <- res2 + halft(x[3*(i-1)+1], x[3*(i-1)+2], x[3*i])
        }
        return(res1*res1)
    }

    x <- c(0.8, 2.3, 4.4, 0.3, 2.0, 2.1)

    R_val <- R_func(x)
    R_grad <- grad(R_func, x, method.args=list(r=8))
    R_hess <- hessian(R_func, x, method.args=list(r=8))
    R_hess_spLT <- tril(drop0(R_hess, 1e-8))
    
    c_list <- cppad_incbeta2(x)
    c_val <- c_list$val
    c_grad <- c_list$grad
    c_hess_dense <- c_list$hess.dense
    c_hess_sp <- c_list$hess.sp
    c_hess_spLT <- tril(c_list$hess.spLT)
 
    expect_equal(c_val, R_val)
    expect_equal(c_grad, R_grad, tolerance=1e-5)
    expect_equal(c_hess_dense, R_hess, tolerance=1e-6)
    expect_equal(c_hess_dense, c_hess_sp)
    expect_equal(c_hess_spLT, tril(drop0(c_hess_sp, 1e-8)))
    expect_equal(c_hess_spLT, R_hess_spLT, tolerance=1e-6)   
})
    
