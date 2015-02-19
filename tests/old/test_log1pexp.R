context("log1pexp")


test_that("log1pexp_sum",{

    R_log1pexp_sum <- function(x) {sum(log1p(exp(x)))}
    x <- c(1.1, 2.2, 3.3)

    R_log1pexp_sum_val <- R_log1pexp_sum(x)
    R_log1pexp_sum_grad <- grad(R_log1pexp_sum, x)
    R_log1pexp_sum_hess <- hessian(R_log1pexp_sum, x)
    R_log1pexp_sum_hess_spLT <- drop0(tril(R_log1pexp_sum_hess), 1e-10)
    
    c_log1pexp_sum_list <- cppad_log1pexp_sum(x)
    c_log1pexp_sum_val <- c_log1pexp_sum_list$val
    c_log1pexp_sum_grad <- c_log1pexp_sum_list$grad
    c_log1pexp_sum_hess_dense <- c_log1pexp_sum_list$hess.dense
    c_log1pexp_sum_hess_sp <- c_log1pexp_sum_list$hess.sp
    c_log1pexp_sum_hess_spLT <- tril(c_log1pexp_sum_list$hess.spLT)

    expect_equal(c_log1pexp_sum_val, R_log1pexp_sum_val)
    expect_equal(c_log1pexp_sum_grad, R_log1pexp_sum_grad)
    expect_equal(c_log1pexp_sum_hess_dense, R_log1pexp_sum_hess)
    expect_equal(c_log1pexp_sum_hess_dense, c_log1pexp_sum_hess_sp)
    expect_equal(c_log1pexp_sum_hess_spLT, tril(drop0(c_log1pexp_sum_hess_sp)))
    expect_equal(c_log1pexp_sum_hess_spLT, R_log1pexp_sum_hess_spLT)   
})


test_that("log1pexp_prod",{

    R_log1pexp_prod <- function(x) {prod(log1p(exp(x)))}
    x <- c(1.1, 2.2, 3.3)

    R_log1pexp_prod_val <- R_log1pexp_prod(x)
    R_log1pexp_prod_grad <- grad(R_log1pexp_prod, x)
    R_log1pexp_prod_hess <- hessian(R_log1pexp_prod, x)
    R_log1pexp_prod_hess_spLT <- drop0(tril(R_log1pexp_prod_hess), 1e-10)
    
    c_log1pexp_prod_list <- cppad_log1pexp_prod(x)
    c_log1pexp_prod_val <- c_log1pexp_prod_list$val
    c_log1pexp_prod_grad <- c_log1pexp_prod_list$grad
    c_log1pexp_prod_hess_dense <- c_log1pexp_prod_list$hess.dense
    c_log1pexp_prod_hess_sp <- c_log1pexp_prod_list$hess.sp
    c_log1pexp_prod_hess_spLT <- tril(c_log1pexp_prod_list$hess.spLT)

    expect_equal(c_log1pexp_prod_val, R_log1pexp_prod_val)
    expect_equal(c_log1pexp_prod_grad, R_log1pexp_prod_grad)
    expect_equal(c_log1pexp_prod_hess_dense, R_log1pexp_prod_hess)
    expect_equal(c_log1pexp_prod_hess_dense, c_log1pexp_prod_hess_sp)
    expect_equal(c_log1pexp_prod_hess_spLT, tril(drop0(c_log1pexp_prod_hess_sp)))
    expect_equal(c_log1pexp_prod_hess_spLT, R_log1pexp_prod_hess_spLT)
})


    
