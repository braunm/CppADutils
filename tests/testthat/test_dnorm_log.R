context("dnorm_log")

test_that("dnorm_log",{
    require(numDeriv)
    R_func <- function(x) {
        n <- length(x)/3
        res <- 0
        for (i in 1:n) {
            res <- res + dnorm(x[3*(i-1)+1], x[3*(i-1)+2], x[3*i], log=TRUE)
        }
        return(res*res)
    }

    x <- c(0.8, 2.3, 4.4, 0.3, 2.0, 2.1, -4, -2, 3)

    R_val <- R_func(x)
    R_grad <- grad(R_func, x, method.args=list(r=8))
    R_hess <- hessian(R_func, x, method.args=list(r=8))
    R_hess_spLT <- tril(drop0(R_hess, 1e-8))
    
    c_list <- cppad_dnorm_log(x)
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
    
