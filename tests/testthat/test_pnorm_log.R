context("pnorm_log")

test_that("pnorm_log",{

    fn <- function(z, m, s) {
        res <- pnorm(z, m, s, log=TRUE, lower.tail=TRUE)
        return(res)
    }
    
    R_func <- function(x) {
        n <- length(x)/3
        res <- 0
        for (i in 1:n) {
            res <- res + fn(x[3*(i-1)+1], x[3*(i-1)+2], x[3*i])
        }
        return(res*res)
    }

    x <- c(1.5, 1.2, 2, 0, .5, 1, -1, -1.2, 1.5)

    R_val <- R_func(x)
    R_grad <- grad(R_func, x, method.args=list(r=8))
    R_hess <- hessian(R_func, x, method.args=list(r=6))
    R_hess_spLT <- tril(drop0(R_hess, 1e-7))
    
    c_list <- cppad_pnorm_log(x)
    c_val <- c_list$val
    c_grad <- c_list$grad
    c_hess_dense <- c_list$hess.dense
    c_hess_sp <- c_list$hess.sp
    c_hess_spLT <- tril(c_list$hess.spLT)
 
    expect_equal(c_val, R_val,tolerance=5e-7)
    expect_equal(c_grad, R_grad, tolerance=1e-5)
    expect_equivalent(drop0(c_hess_dense), drop0(R_hess,1e-8))
    expect_equal(c_hess_dense, c_hess_sp)
    expect_equal(c_hess_spLT, tril(drop0(c_hess_sp, 1e-8)),tolerance=1e-7)
    expect_equal(c_hess_spLT, R_hess_spLT, tolerance=1e-5)   
})

context("pnorm")

test_that("pnorm",{

    fn <- function(z, m, s) {
        res <- pnorm(z, m, s, log=FALSE, lower.tail=TRUE)
        return(res)
    }
    
    R_func <- function(x) {
        n <- length(x)/3
        res <- 0
        for (i in 1:n) {
            res <- res + fn(x[3*(i-1)+1], x[3*(i-1)+2], x[3*i])
        }
        return(res*res)
    }

    x <- c(1.5, 1.2, 2, 0, .5, 1, -1, -1.2, 1.5)

    R_val <- R_func(x)
    R_grad <- grad(R_func, x, method.args=list(r=8))
    R_hess <- hessian(R_func, x, method.args=list(r=6))
    R_hess_spLT <- tril(drop0(R_hess, 1e-7))
    
    c_list <- cppad_pnorm(x)
    c_val <- c_list$val
    c_grad <- c_list$grad
    c_hess_dense <- c_list$hess.dense
    c_hess_sp <- c_list$hess.sp
    c_hess_spLT <- tril(c_list$hess.spLT)
 
    expect_equal(c_val, R_val,tolerance=5e-7)
    expect_equal(c_grad, R_grad, tolerance=1e-5)
    expect_equivalent(drop0(c_hess_dense), drop0(R_hess,1e-8))
    expect_equal(c_hess_dense, c_hess_sp)
    expect_equal(c_hess_spLT, tril(drop0(c_hess_sp, 1e-8)),tolerance=1e-7)
    expect_equal(c_hess_spLT, R_hess_spLT, tolerance=1e-5)   
})
    
