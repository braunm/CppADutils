context("loginvlogit")


test_that("loginvlogit_sum",{

    R_loginvlogit_sum <- function(x) {sum(x-log1pexp(x))}
    x <- c(1.1, 2.2, 3.3)

    R_loginvlogit_sum_val <- R_loginvlogit_sum(x)
    R_loginvlogit_sum_grad <- grad(R_loginvlogit_sum, x)
    R_loginvlogit_sum_hess <- hessian(R_loginvlogit_sum, x, method.args=list(r=6))
    R_loginvlogit_sum_hess_spLT <- drop0(tril(R_loginvlogit_sum_hess), 1e-10)
    
    c_loginvlogit_sum_list <- cppad_loginvlogit_sum(x)
    c_loginvlogit_sum_val <- c_loginvlogit_sum_list$val
    c_loginvlogit_sum_grad <- c_loginvlogit_sum_list$grad
    c_loginvlogit_sum_hess_dense <- c_loginvlogit_sum_list$hess.dense
    c_loginvlogit_sum_hess_sp <- c_loginvlogit_sum_list$hess.sp
    c_loginvlogit_sum_hess_spLT <- tril(c_loginvlogit_sum_list$hess.spLT)

    expect_equal(c_loginvlogit_sum_val, R_loginvlogit_sum_val)
    expect_equal(c_loginvlogit_sum_grad, R_loginvlogit_sum_grad)
    expect_equal(c_loginvlogit_sum_hess_dense, R_loginvlogit_sum_hess)
    expect_equal(c_loginvlogit_sum_hess_dense, c_loginvlogit_sum_hess_sp)
    expect_equal(c_loginvlogit_sum_hess_spLT, tril(drop0(c_loginvlogit_sum_hess_sp)))
    expect_equal(c_loginvlogit_sum_hess_spLT, R_loginvlogit_sum_hess_spLT)   
})


test_that("loginvlogit_prod",{

    R_loginvlogit_prod <- function(x) {prod(x-log1pexp(x))}
    x <- c(1.1, 2.2, 3.3)

    R_loginvlogit_prod_val <- R_loginvlogit_prod(x)
    R_loginvlogit_prod_grad <- grad(R_loginvlogit_prod, x)
    R_loginvlogit_prod_hess <- hessian(R_loginvlogit_prod, x, method.args=list(r=6))
    R_loginvlogit_prod_hess_spLT <- drop0(tril(R_loginvlogit_prod_hess), 1e-10)
    
    c_loginvlogit_prod_list <- cppad_loginvlogit_prod(x)
    c_loginvlogit_prod_val <- c_loginvlogit_prod_list$val
    c_loginvlogit_prod_grad <- c_loginvlogit_prod_list$grad
    c_loginvlogit_prod_hess_dense <- c_loginvlogit_prod_list$hess.dense
    c_loginvlogit_prod_hess_sp <- c_loginvlogit_prod_list$hess.sp
    c_loginvlogit_prod_hess_spLT <- tril(c_loginvlogit_prod_list$hess.spLT)

    expect_equal(c_loginvlogit_prod_val, R_loginvlogit_prod_val)
    expect_equal(c_loginvlogit_prod_grad, R_loginvlogit_prod_grad)
    expect_equal(c_loginvlogit_prod_hess_dense, R_loginvlogit_prod_hess)
    expect_equal(c_loginvlogit_prod_hess_dense, c_loginvlogit_prod_hess_sp)
    expect_equal(c_loginvlogit_prod_hess_spLT, tril(drop0(c_loginvlogit_prod_hess_sp)))
    expect_equal(c_loginvlogit_prod_hess_spLT, R_loginvlogit_prod_hess_spLT)
})


    
