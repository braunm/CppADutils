context("dnorm_log")


test_that("dnorm_log1",{

    R_dnorm_log1 <- function(x) {dnorm(x[1], x[2], x[3], log=TRUE)}
    x <- c(2.9, 2.5, 1.1)

    R_val <- R_dnorm_log1(x)
    R_grad <- grad(R_dnorm_log1, x)
    R_hess <- hessian(R_dnorm_log1, x, method.args=list(r=8))
    R_hess_spLT <- drop0(tril(R_hess), 1e-10)
    
    c_list <- cppad_dnorm_log1(x)
    c_val <- c_list$val
    c_grad <- c_list$grad
    c_hess_dense <- c_list$hess.dense
    c_hess_sp <- c_list$hess.sp
    c_hess_spLT <- tril(c_list$hess.spLT)
    
    expect_equal(c_val, R_val)
    expect_equal(c_grad, R_grad)
    expect_equal(c_hess_dense, R_hess)
    expect_equal(c_hess_dense, c_hess_sp)
    expect_equal(c_hess_spLT, tril(drop0(c_hess_sp, 1e-10)))
    expect_equal(c_hess_spLT, R_hess_spLT)  
})


test_that("dnorm_log_sum",{

    R_dnorm_log_sum <- function(x) {
        n <- floor(length(x)/3)
        res <- 0.0
        for (i in 1:n) {
            res <- res + dnorm(x[3*(i-1)+1], x[3*(i-1)+2], x[3*i], log=TRUE)
        }
        return(res)
    }
    x <- c(2.9, 2.5, 1.1, 0, 0, 1, -.5, .1, 2)


    R_val <- R_dnorm_log_sum(x)
    R_grad <- grad(R_dnorm_log_sum, x)
    R_hess <- hessian(R_dnorm_log_sum, x, method.args=list(r=8))
    R_hess_spLT <- tril(drop0(R_hess, 1e-7))
    
    c_list <- cppad_dnorm_log_sum(x)
    c_val <- c_list$val
    c_grad <- c_list$grad
    c_hess_dense <- c_list$hess.dense
    c_hess_sp <- c_list$hess.sp
    c_hess_spLT <- tril(drop0(c_list$hess.spLT,1e-10))

    expect_equal(c_val, R_val)
    expect_equal(c_grad, R_grad)
    expect_equal(c_hess_dense, R_hess)
    expect_equal(c_hess_dense, c_hess_sp)
    expect_equal(c_hess_spLT, tril(drop0(c_hess_sp, 1e-9)))
    expect_equal(c_hess_spLT, R_hess_spLT)   
})

    
