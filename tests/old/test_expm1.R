context("expm1")


test_that("expm1_sum",{

    R_expm1_sum <- function(x) {sum(expm1(x))}
    x <- c(1.1, 2.2, 3.3)

    R_expm1_sum_val <- R_expm1_sum(x)
    R_expm1_sum_grad <- grad(R_expm1_sum, x)
    R_expm1_sum_hess <- hessian(R_expm1_sum, x, method.args=list(r=6))
    R_expm1_sum_hess_spLT <- drop0(tril(R_expm1_sum_hess), 1e-6)
    
    c_expm1_sum_list <- cppad_expm1_sum(x)
    c_expm1_sum_val <- c_expm1_sum_list$val
    c_expm1_sum_grad <- c_expm1_sum_list$grad
    c_expm1_sum_hess_dense <- c_expm1_sum_list$hess.dense
    c_expm1_sum_hess_sp <- c_expm1_sum_list$hess.sp
    c_expm1_sum_hess_spLT <- tril(c_expm1_sum_list$hess.spLT)

    expect_equal(c_expm1_sum_val, R_expm1_sum_val)
    expect_equal(c_expm1_sum_grad, R_expm1_sum_grad)
    expect_equal(c_expm1_sum_hess_dense, R_expm1_sum_hess)
    expect_equal(c_expm1_sum_hess_dense, c_expm1_sum_hess_sp)
    expect_equal(c_expm1_sum_hess_spLT, tril(drop0(c_expm1_sum_hess_sp)))
    expect_equal(c_expm1_sum_hess_spLT, R_expm1_sum_hess_spLT)   
})


test_that("expm1_prod",{

    R_expm1_prod <- function(x) {prod(expm1(x))}
    x <- c(1.1, 2.2, 3.3)

    R_expm1_prod_val <- R_expm1_prod(x)
    R_expm1_prod_grad <- grad(R_expm1_prod, x)
    R_expm1_prod_hess <- hessian(R_expm1_prod, x, method.args=list(r=6))
    R_expm1_prod_hess_spLT <- drop0(tril(R_expm1_prod_hess), 1e-6)
    
    c_expm1_prod_list <- cppad_expm1_prod(x)
    c_expm1_prod_val <- c_expm1_prod_list$val
    c_expm1_prod_grad <- c_expm1_prod_list$grad
    c_expm1_prod_hess_dense <- c_expm1_prod_list$hess.dense
    c_expm1_prod_hess_sp <- c_expm1_prod_list$hess.sp
    c_expm1_prod_hess_spLT <- tril(c_expm1_prod_list$hess.spLT)

    expect_equal(c_expm1_prod_val, R_expm1_prod_val)
    expect_equal(c_expm1_prod_grad, R_expm1_prod_grad)
    expect_equal(c_expm1_prod_hess_dense, R_expm1_prod_hess)
    expect_equal(c_expm1_prod_hess_dense, c_expm1_prod_hess_sp)
    expect_equal(c_expm1_prod_hess_spLT, tril(drop0(c_expm1_prod_hess_sp)))
    expect_equal(c_expm1_prod_hess_spLT, R_expm1_prod_hess_spLT)
})


    
