context("lgammaexp")


test_that("lgammaexp_sum",{

    R_lgammaexp_sum <- function(x) {sum(lgammaexp(x))}
    x <- c(-2, -1, 3.2, 10)

    R_lgammaexp_sum_val <- R_lgammaexp_sum(x)
    R_lgammaexp_sum_grad <- grad(R_lgammaexp_sum, x)
    R_lgammaexp_sum_hess <- hessian(R_lgammaexp_sum, x, method.args=list(r=6))
    R_lgammaexp_sum_hess_spLT <- drop0(tril(R_lgammaexp_sum_hess), 1e-6)
    
    c_lgammaexp_sum_list <- cppad_lgammaexp_sum(x)
    c_lgammaexp_sum_val <- c_lgammaexp_sum_list$val
    c_lgammaexp_sum_grad <- c_lgammaexp_sum_list$grad
    c_lgammaexp_sum_hess_dense <- c_lgammaexp_sum_list$hess.dense
    c_lgammaexp_sum_hess_sp <- c_lgammaexp_sum_list$hess.sp
    c_lgammaexp_sum_hess_spLT <- tril(c_lgammaexp_sum_list$hess.spLT)

    expect_equal(c_lgammaexp_sum_val, R_lgammaexp_sum_val)
    expect_equal(c_lgammaexp_sum_grad, R_lgammaexp_sum_grad)
    expect_equal(c_lgammaexp_sum_hess_dense, R_lgammaexp_sum_hess)
    expect_equal(c_lgammaexp_sum_hess_dense, c_lgammaexp_sum_hess_sp)
    expect_equal(c_lgammaexp_sum_hess_spLT, tril(drop0(c_lgammaexp_sum_hess_sp)))
    expect_equal(c_lgammaexp_sum_hess_spLT, R_lgammaexp_sum_hess_spLT)   
})


test_that("lgammaexp_prod",{

    R_lgammaexp_prod <- function(x) {prod(lgammaexp(x))}
    x <- c(-2, -1, 3.2, 10)

    R_lgammaexp_prod_val <- R_lgammaexp_prod(x)
    R_lgammaexp_prod_grad <- grad(R_lgammaexp_prod, x)
    R_lgammaexp_prod_hess <- hessian(R_lgammaexp_prod, x, method.args=list(r=6))
    R_lgammaexp_prod_hess_spLT <- drop0(tril(R_lgammaexp_prod_hess), 1e-6)
    
    c_lgammaexp_prod_list <- cppad_lgammaexp_prod(x)
    c_lgammaexp_prod_val <- c_lgammaexp_prod_list$val
    c_lgammaexp_prod_grad <- c_lgammaexp_prod_list$grad
    c_lgammaexp_prod_hess_dense <- c_lgammaexp_prod_list$hess.dense
    c_lgammaexp_prod_hess_sp <- c_lgammaexp_prod_list$hess.sp
    c_lgammaexp_prod_hess_spLT <- tril(c_lgammaexp_prod_list$hess.spLT)

    expect_equal(c_lgammaexp_prod_val, R_lgammaexp_prod_val)
    expect_equal(c_lgammaexp_prod_grad, R_lgammaexp_prod_grad)
    expect_equal(c_lgammaexp_prod_hess_dense, R_lgammaexp_prod_hess)
    expect_equal(c_lgammaexp_prod_hess_dense, c_lgammaexp_prod_hess_sp)
    expect_equal(c_lgammaexp_prod_hess_spLT, tril(drop0(c_lgammaexp_prod_hess_sp)))
    expect_equal(c_lgammaexp_prod_hess_spLT, R_lgammaexp_prod_hess_spLT)
})


    
