context("inc_beta_reg")


test_that("inc_beta_reg1",{

    R_inc_beta_reg1 <- function(x) {inc_beta_reg(x[1], x[2], x[3])}
    x <- c(.4, 2.2, 5.5)

    R_inc_beta_reg1_val <- R_inc_beta_reg1(x)
    R_inc_beta_reg1_grad <- grad(R_inc_beta_reg1, x)
    R_inc_beta_reg1_hess <- hessian(R_inc_beta_reg1, x, method.args=list(r=6))
    R_inc_beta_reg1_hess_spLT <- drop0(tril(R_inc_beta_reg1_hess), 1e-10)
    
    c_inc_beta_reg1_list <- cppad_inc_beta_reg1(x)
    c_inc_beta_reg1_val <- c_inc_beta_reg1_list$val
    c_inc_beta_reg1_grad <- c_inc_beta_reg1_list$grad
    c_inc_beta_reg1_hess_dense <- c_inc_beta_reg1_list$hess.dense
    c_inc_beta_reg1_hess_sp <- c_inc_beta_reg1_list$hess.sp
    c_inc_beta_reg1_hess_spLT <- tril(c_inc_beta_reg1_list$hess.spLT)
 
    expect_equal(c_inc_beta_reg1_val, R_inc_beta_reg1_val)
    expect_equal(c_inc_beta_reg1_grad, R_inc_beta_reg1_grad)
    expect_equal(c_inc_beta_reg1_hess_dense, R_inc_beta_reg1_hess)
    expect_equal(c_inc_beta_reg1_hess_dense, c_inc_beta_reg1_hess_sp)
    expect_equal(c_inc_beta_reg1_hess_spLT, tril(drop0(c_inc_beta_reg1_hess_sp)))
    expect_equal(c_inc_beta_reg1_hess_spLT, R_inc_beta_reg1_hess_spLT)   
})


test_that("inc_beta_reg_sum",{

    R_inc_beta_reg_sum <- function(x) {
        n <- floor(length(x)/3)
        res <- 0.0
        for (i in 1:n) {
            res <- res + inc_beta_reg(x[3*(i-1)+1], x[3*(i-1)+2], x[3*i])
        }
        return(res)
    }
    x <- x <- c(.4, 2.2, 5.5, .02, .01, 12, .9, 23, 23)

    R_inc_beta_reg_sum_val <- R_inc_beta_reg_sum(x)
    R_inc_beta_reg_sum_grad <- grad(R_inc_beta_reg_sum, x)
    R_inc_beta_reg_sum_hess <- hessian(R_inc_beta_reg_sum, x, method.args=list(r=6))
    R_inc_beta_reg_sum_hess_spLT <- drop0(tril(R_inc_beta_reg_sum_hess), 1e-10)
    
    c_inc_beta_reg_sum_list <- cppad_inc_beta_reg_sum(x)
    c_inc_beta_reg_sum_val <- c_inc_beta_reg_sum_list$val
    c_inc_beta_reg_sum_grad <- c_inc_beta_reg_sum_list$grad
    c_inc_beta_reg_sum_hess_dense <- c_inc_beta_reg_sum_list$hess.dense
    c_inc_beta_reg_sum_hess_sp <- c_inc_beta_reg_sum_list$hess.sp
    c_inc_beta_reg_sum_hess_spLT <- tril(c_inc_beta_reg_sum_list$hess.spLT)

    expect_equal(c_inc_beta_reg_sum_val, R_inc_beta_reg_sum_val)
    expect_equal(c_inc_beta_reg_sum_grad, R_inc_beta_reg_sum_grad)
    expect_equal(c_inc_beta_reg_sum_hess_dense, R_inc_beta_reg_sum_hess)
    expect_equal(c_inc_beta_reg_sum_hess_dense, c_inc_beta_reg_sum_hess_sp)
    expect_equal(c_inc_beta_reg_sum_hess_spLT, tril(drop0(c_inc_beta_reg_sum_hess_sp, 1e-10)))
    expect_equal(c_inc_beta_reg_sum_hess_spLT, R_inc_beta_reg_sum_hess_spLT)   
})

    
