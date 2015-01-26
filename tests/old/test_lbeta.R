context("lbeta")


test_that("lbeta1",{

    R_lbeta1 <- function(x) {lbeta(x[1], x[2])}
    x <- c(.4, 2.2)

    R_lbeta1_val <- R_lbeta1(x)
    R_lbeta1_grad <- grad(R_lbeta1, x)
    R_lbeta1_hess <- hessian(R_lbeta1, x, method.args=list(r=6))
    R_lbeta1_hess_spLT <- drop0(tril(R_lbeta1_hess), 1e-10)
    
    c_lbeta1_list <- cppad_lbeta1(x)
    c_lbeta1_val <- c_lbeta1_list$val
    c_lbeta1_grad <- c_lbeta1_list$grad
    c_lbeta1_hess_dense <- c_lbeta1_list$hess.dense
    c_lbeta1_hess_sp <- c_lbeta1_list$hess.sp
    c_lbeta1_hess_spLT <- tril(c_lbeta1_list$hess.spLT)
 
    expect_equal(c_lbeta1_val, R_lbeta1_val)
    expect_equal(c_lbeta1_grad, R_lbeta1_grad)
    expect_equal(c_lbeta1_hess_dense, R_lbeta1_hess)
    expect_equal(c_lbeta1_hess_dense, c_lbeta1_hess_sp)
    expect_equal(c_lbeta1_hess_spLT, tril(drop0(c_lbeta1_hess_sp)))
    expect_equal(c_lbeta1_hess_spLT, R_lbeta1_hess_spLT)   
})


test_that("lbeta_sum",{

    R_lbeta_sum <- function(x) {
        n <- floor(length(x)/2)
        res <- 0.0
        for (i in 1:n) {
            res <- res + lbeta(x[2*(i-1)+1], x[2*i])
        }
        return(res)
    }
    x <- seq(2:7) + 0.1

    R_lbeta_sum_val <- R_lbeta_sum(x)
    R_lbeta_sum_grad <- grad(R_lbeta_sum, x)
    R_lbeta_sum_hess <- hessian(R_lbeta_sum, x, method.args=list(r=6))
    R_lbeta_sum_hess_spLT <- drop0(tril(R_lbeta_sum_hess), 1e-10)
    
    c_lbeta_sum_list <- cppad_lbeta_sum(x)
    c_lbeta_sum_val <- c_lbeta_sum_list$val
    c_lbeta_sum_grad <- c_lbeta_sum_list$grad
    c_lbeta_sum_hess_dense <- c_lbeta_sum_list$hess.dense
    c_lbeta_sum_hess_sp <- c_lbeta_sum_list$hess.sp
    c_lbeta_sum_hess_spLT <- tril(c_lbeta_sum_list$hess.spLT)

    expect_equal(c_lbeta_sum_val, R_lbeta_sum_val)
    expect_equal(c_lbeta_sum_grad, R_lbeta_sum_grad)
    expect_equal(c_lbeta_sum_hess_dense, R_lbeta_sum_hess)
    expect_equal(c_lbeta_sum_hess_dense, c_lbeta_sum_hess_sp)
    expect_equal(c_lbeta_sum_hess_spLT, tril(drop0(c_lbeta_sum_hess_sp, 1e-10)))
    expect_equal(c_lbeta_sum_hess_spLT, R_lbeta_sum_hess_spLT)   
})

    
