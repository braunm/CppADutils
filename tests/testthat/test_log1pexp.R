context("log1pexp")

test_that("log1pexp",{
    require(numDeriv)
    R_func <- function(x) {
        n <- length(x)/2
        res <- 0
        for (i in 1:n) {
            res <- res + log1pexp(x[2*(i-1)+1]) * lgamma(log1pexp(x[2*i]))
        }
        return(res)
    }

    x <- c(.4, 2.2, 4.4, 3.3, 5, 1.5)

    R_val <- R_func(x)
    R_grad <- grad(R_func, x)
    R_hess <- hessian(R_func, x, method.args=list(r=8))
    R_hess_spLT <- tril(drop0(R_hess, 1e-8))
    
    c_list <- cppad_log1pexp(x)
    c_val <- c_list$val
    c_grad <- c_list$grad
    c_hess_dense <- c_list$hess.dense
    c_hess_sp <- c_list$hess.sp
    c_hess_spLT <- tril(c_list$hess.spLT)
 
    expect_equal(c_val, R_val)
    expect_equal(c_grad, R_grad)
    expect_equal(c_hess_dense, R_hess)
    expect_equal(c_hess_dense, c_hess_sp)
    expect_equal(c_hess_spLT, tril(drop0(c_hess_sp, 1e-8)))
    expect_equal(c_hess_spLT, R_hess_spLT)   
})
    
