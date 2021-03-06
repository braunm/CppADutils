context("dt_log")

test_that("dt_log",{
    require(numDeriv)
    fn <- function(z, v, s) {
        res <- lgamma((v+1)/2) - lgamma(v/2)
        res <- res - 0.5*((v+1)*log(1+z^2/(v*s^2)) + log(pi * v)) - log(s)
        return(res)
    }
    
    R_func <- function(x) {
        n <- length(x)/3
        res <- 0
        for (i in 1:n) {
            res <- res + fn(x[3*(i-1)+1], x[3*(i-1)+2], x[3*i])
        }
        return(res)
    }

    x <- c(0.8, 4, 1, 0.3, 5, 1, -1, 2, 1.5)

    R_val <- R_func(x)
    R_grad <- grad(R_func, x, method.args=list(r=8))
    R_hess <- hessian(R_func, x, method.args=list(r=8))
    R_hess_spLT <- tril(drop0(R_hess, 1e-8))
    
    c_list <- cppad_dt_log(x)
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
    
