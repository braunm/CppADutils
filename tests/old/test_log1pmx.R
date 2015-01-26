context("log1pmx")


test_that("log1pmx_sum",{

    R_log1pmx_sum <- function(x) {sum(log1p(x)-x)}
    x <- c(1.1, 2.2, 3.3)

    R_log1pmx_sum_val <- R_log1pmx_sum(x)
    R_log1pmx_sum_grad <- 1/(1+x)-1
    R_log1pmx_sum_hess <- hessian(R_log1pmx_sum, x, method.args=list(r=6))
    R_log1pmx_sum_hess_spLT <- drop0(tril(R_log1pmx_sum_hess), 1e-10)
    
    c_log1pmx_sum_list <- cppad_log1pmx_sum(x)
    c_log1pmx_sum_val <- c_log1pmx_sum_list$val
    c_log1pmx_sum_grad <- c_log1pmx_sum_list$grad
    c_log1pmx_sum_hess_dense <- c_log1pmx_sum_list$hess.dense
    c_log1pmx_sum_hess_sp <- c_log1pmx_sum_list$hess.sp
    c_log1pmx_sum_hess_spLT <- tril(c_log1pmx_sum_list$hess.spLT)
 
    expect_equal(c_log1pmx_sum_val, R_log1pmx_sum_val)
    expect_equal(c_log1pmx_sum_grad, R_log1pmx_sum_grad)
    expect_equal(c_log1pmx_sum_hess_dense, R_log1pmx_sum_hess)
    expect_equal(c_log1pmx_sum_hess_dense, c_log1pmx_sum_hess_sp)
    expect_equal(c_log1pmx_sum_hess_spLT, tril(drop0(c_log1pmx_sum_hess_sp)))
    expect_equal(c_log1pmx_sum_hess_spLT, R_log1pmx_sum_hess_spLT)   
 })


test_that("log1pmx_prod",{

    R_log1pmx_prod <- function(x) {prod(log1p(x)-x)}
    x <- c(1.1, 2.2, 3.3)

    R_log1pmx_prod_val <- R_log1pmx_prod(x)
    R_log1pmx_prod_grad <- -x*R_log1pmx_prod_val/((log1p(x)-x)*(1+x))
    R_log1pmx_prod_hess <- hessian(R_log1pmx_prod, x, method.args=list(r=6))
    R_log1pmx_prod_hess_spLT <- drop0(tril(R_log1pmx_prod_hess), 1e-10)
    
    c_log1pmx_prod_list <- cppad_log1pmx_prod(x)
    c_log1pmx_prod_val <- c_log1pmx_prod_list$val
    c_log1pmx_prod_grad <- c_log1pmx_prod_list$grad
    c_log1pmx_prod_hess_dense <- c_log1pmx_prod_list$hess.dense
    c_log1pmx_prod_hess_sp <- c_log1pmx_prod_list$hess.sp
    c_log1pmx_prod_hess_spLT <- tril(c_log1pmx_prod_list$hess.spLT)
  
    expect_equal(c_log1pmx_prod_val, R_log1pmx_prod_val)
    expect_equal(c_log1pmx_prod_grad, R_log1pmx_prod_grad)
    expect_equal(c_log1pmx_prod_hess_dense, R_log1pmx_prod_hess)
    expect_equal(c_log1pmx_prod_hess_dense, c_log1pmx_prod_hess_sp)
    expect_equal(c_log1pmx_prod_hess_spLT, tril(drop0(c_log1pmx_prod_hess_sp)))
    expect_equal(c_log1pmx_prod_hess_spLT, R_log1pmx_prod_hess_spLT)
})


    
