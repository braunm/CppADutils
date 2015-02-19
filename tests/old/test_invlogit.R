context("invlogit")


test_that("invlogit_sum",{

    R_invlogit_sum <- function(x) {sum(exp(x-log1pexp(x)))}
    x <- c(1.1, 2.2, 3.3)

    R_invlogit_sum_val <- R_invlogit_sum(x)
    R_invlogit_sum_grad <- grad(R_invlogit_sum, x)
    R_invlogit_sum_hess <- hessian(R_invlogit_sum, x, method.args=list(r=6))
    R_invlogit_sum_hess_spLT <- drop0(tril(R_invlogit_sum_hess), 1e-10)
    
    c_invlogit_sum_list <- cppad_invlogit_sum(x)
    c_invlogit_sum_val <- c_invlogit_sum_list$val
    c_invlogit_sum_grad <- c_invlogit_sum_list$grad
    c_invlogit_sum_hess_dense <- c_invlogit_sum_list$hess.dense
    c_invlogit_sum_hess_sp <- c_invlogit_sum_list$hess.sp
    c_invlogit_sum_hess_spLT <- tril(c_invlogit_sum_list$hess.spLT)
 
    expect_equal(c_invlogit_sum_val, R_invlogit_sum_val)
    expect_equal(c_invlogit_sum_grad, R_invlogit_sum_grad)
    expect_equal(c_invlogit_sum_hess_dense, R_invlogit_sum_hess)
    expect_equal(c_invlogit_sum_hess_dense, c_invlogit_sum_hess_sp)
    expect_equal(c_invlogit_sum_hess_spLT, tril(drop0(c_invlogit_sum_hess_sp)))
    expect_equal(c_invlogit_sum_hess_spLT, R_invlogit_sum_hess_spLT)   
})


test_that("invlogit_prod",{

    R_invlogit_prod <- function(x) {prod(exp(x-log1pexp(x)))}
    x <- c(1.1, 2.2, 3.3)

    R_invlogit_prod_val <- R_invlogit_prod(x)
    R_invlogit_prod_grad <- grad(R_invlogit_prod, x)
    R_invlogit_prod_hess <- hessian(R_invlogit_prod, x, method.args=list(r=6))
    R_invlogit_prod_hess_spLT <- drop0(tril(R_invlogit_prod_hess), 1e-10)
    
    c_invlogit_prod_list <- cppad_invlogit_prod(x)
    c_invlogit_prod_val <- c_invlogit_prod_list$val
    c_invlogit_prod_grad <- c_invlogit_prod_list$grad
    c_invlogit_prod_hess_dense <- c_invlogit_prod_list$hess.dense
    c_invlogit_prod_hess_sp <- c_invlogit_prod_list$hess.sp
    c_invlogit_prod_hess_spLT <- tril(c_invlogit_prod_list$hess.spLT)

    expect_equal(c_invlogit_prod_val, R_invlogit_prod_val)
    expect_equal(c_invlogit_prod_grad, R_invlogit_prod_grad)
    expect_equal(c_invlogit_prod_hess_dense, R_invlogit_prod_hess)
    expect_equal(c_invlogit_prod_hess_dense, c_invlogit_prod_hess_sp)
    expect_equal(c_invlogit_prod_hess_spLT, tril(drop0(c_invlogit_prod_hess_sp)))
    expect_equal(c_invlogit_prod_hess_spLT, R_invlogit_prod_hess_spLT)
})


    
