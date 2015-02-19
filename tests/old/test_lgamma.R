context("lgamma")


test_that("lgamma_sum",{

    R_lgamma_sum <- function(x) {sum(lgamma(x))}
    x <- c(1.1, 2.2, 3.3)

    R_lgamma_sum_val <- R_lgamma_sum(x)
    R_lgamma_sum_grad <- grad(R_lgamma_sum, x)
    R_lgamma_sum_hess <- hessian(R_lgamma_sum, x)
    R_lgamma_sum_hess_spLT <- drop0(tril(R_lgamma_sum_hess), 1e-10)
    
    c_lgamma_sum_list <- cppad_lgamma_sum(x)
    c_lgamma_sum_val <- c_lgamma_sum_list$val
    c_lgamma_sum_grad <- c_lgamma_sum_list$grad
    c_lgamma_sum_hess_dense <- c_lgamma_sum_list$hess.dense
    c_lgamma_sum_hess_sp <- c_lgamma_sum_list$hess.sp
    c_lgamma_sum_hess_spLT <- tril(c_lgamma_sum_list$hess.spLT)
  
    expect_equal(c_lgamma_sum_val, R_lgamma_sum_val)
    expect_equal(c_lgamma_sum_grad, R_lgamma_sum_grad)
    expect_equal(c_lgamma_sum_hess_dense, R_lgamma_sum_hess)
    expect_equal(c_lgamma_sum_hess_dense, c_lgamma_sum_hess_sp)
    expect_equal(c_lgamma_sum_hess_spLT, tril(drop0(c_lgamma_sum_hess_sp)))
    expect_equal(c_lgamma_sum_hess_spLT, R_lgamma_sum_hess_spLT)   
})


test_that("lgamma_prod",{

    R_lgamma_prod <- function(x) {prod(lgamma(x))}
    x <- c(1.1, 2.2, 3.3)

    R_lgamma_prod_val <- R_lgamma_prod(x)
    R_lgamma_prod_grad <- grad(R_lgamma_prod, x)
    R_lgamma_prod_hess <- hessian(R_lgamma_prod, x)
    R_lgamma_prod_hess_spLT <- drop0(tril(R_lgamma_prod_hess), 1e-10)
    
    c_lgamma_prod_list <- cppad_lgamma_prod(x)
    c_lgamma_prod_val <- c_lgamma_prod_list$val
    c_lgamma_prod_grad <- c_lgamma_prod_list$grad
    c_lgamma_prod_hess_dense <- c_lgamma_prod_list$hess.dense
    c_lgamma_prod_hess_sp <- c_lgamma_prod_list$hess.sp
    c_lgamma_prod_hess_spLT <- tril(c_lgamma_prod_list$hess.spLT)

    expect_equal(c_lgamma_prod_val, R_lgamma_prod_val)
    expect_equal(c_lgamma_prod_grad, R_lgamma_prod_grad)
    expect_equal(c_lgamma_prod_hess_dense, R_lgamma_prod_hess)
    expect_equal(c_lgamma_prod_hess_dense, c_lgamma_prod_hess_sp)
    expect_equal(c_lgamma_prod_hess_spLT, tril(drop0(c_lgamma_prod_hess_sp)))
    expect_equal(c_lgamma_prod_hess_spLT, R_lgamma_prod_hess_spLT)
})


    
