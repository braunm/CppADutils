context("comp1")

test_that("comp1_sum",{

    x <- c(2.1, 3.2, 4.3, 1.9)   
    R_comp1_sum <- function(x) {log1p(x[1]*x[4]^2) + sum(expm1(lgamma1p(x[1:3])))}
 
    R_comp1_sum_val <- R_comp1_sum(x)
    R_comp1_sum_grad <- grad(R_comp1_sum, x)
    R_comp1_sum_hess <- hessian(R_comp1_sum, x, method.args=list(r=6))
    R_comp1_sum_hess_spLT <- drop0(tril(R_comp1_sum_hess), 1e-10)
    
    c_comp1_sum_list <- cppad_comp1_sum(x)
    c_comp1_sum_val <- c_comp1_sum_list$val
    c_comp1_sum_grad <- c_comp1_sum_list$grad
    c_comp1_sum_hess_dense <- c_comp1_sum_list$hess.dense
    c_comp1_sum_hess_sp <- c_comp1_sum_list$hess.sp
    c_comp1_sum_hess_spLT <- tril(c_comp1_sum_list$hess.spLT)
 
    expect_equal(c_comp1_sum_val, R_comp1_sum_val)
    expect_equal(c_comp1_sum_grad, R_comp1_sum_grad)
    expect_equal(c_comp1_sum_hess_dense, R_comp1_sum_hess)
    expect_equal(c_comp1_sum_hess_dense, c_comp1_sum_hess_sp)
     expect_equal(c_comp1_sum_hess_spLT, tril(drop0(c_comp1_sum_hess_sp)))
    expect_equal(c_comp1_sum_hess_spLT, R_comp1_sum_hess_spLT)   
})


test_that("comp1_prod",{

    x <- c(2.1, 3.2, 4.3, 1.9)
    R_comp1_prod <- function(x) {log1p(x[4]^2) + prod(expm1(lgamma1p(x[1:3])))}

    R_comp1_prod_val <- R_comp1_prod(x)
    R_comp1_prod_grad <- grad(R_comp1_prod, x)
    R_comp1_prod_hess <- hessian(R_comp1_prod, x, method.args=list(r=6))
    R_comp1_prod_hess_spLT <- drop0(tril(R_comp1_prod_hess), 1e-10)
    
    c_comp1_prod_list <- cppad_comp1_prod(x)
    c_comp1_prod_val <- c_comp1_prod_list$val
    c_comp1_prod_grad <- c_comp1_prod_list$grad
    c_comp1_prod_hess_dense <- c_comp1_prod_list$hess.dense
    c_comp1_prod_hess_sp <- c_comp1_prod_list$hess.sp
    c_comp1_prod_hess_spLT <- tril(c_comp1_prod_list$hess.spLT)
  
    expect_equal(c_comp1_prod_val, R_comp1_prod_val)
    expect_equal(c_comp1_prod_grad, R_comp1_prod_grad)
    expect_equal(c_comp1_prod_hess_dense, R_comp1_prod_hess)
    expect_equal(c_comp1_prod_hess_dense, c_comp1_prod_hess_sp)
    expect_equal(c_comp1_prod_hess_spLT, tril(drop0(c_comp1_prod_hess_sp)))
    expect_equal(c_comp1_prod_hess_spLT, R_comp1_prod_hess_spLT)

})

test_that("comp1b",{

    x <- c(2.1, 3.2, 4.3, 1.9)   
    R_comp1b <- function(x) {
        a <- log1p(x[1]*x[4]^2) + sum(expm1(lgamma1p(x[1:3])))
        return(invlogit(0.4)*a)
    }
 
    R_comp1b_val <- R_comp1b(x)
    R_comp1b_grad <- grad(R_comp1b, x)
    R_comp1b_hess <- hessian(R_comp1b, x, method.args=list(r=6))
    R_comp1b_hess_spLT <- drop0(tril(R_comp1b_hess), 1e-10)
    
    c_comp1b_list <- cppad_comp1b(x)
    c_comp1b_val <- c_comp1b_list$val
    c_comp1b_grad <- c_comp1b_list$grad
    c_comp1b_hess_dense <- c_comp1b_list$hess.dense
    c_comp1b_hess_sp <- c_comp1b_list$hess.sp
    c_comp1b_hess_spLT <- tril(c_comp1b_list$hess.spLT)
  
    expect_equal(c_comp1b_val, R_comp1b_val)
    expect_equal(c_comp1b_grad, R_comp1b_grad)
    expect_equal(c_comp1b_hess_dense, R_comp1b_hess)
    expect_equal(c_comp1b_hess_dense, c_comp1b_hess_sp)
    expect_equal(c_comp1b_hess_spLT, tril(drop0(c_comp1b_hess_sp)))
    expect_equal(c_comp1b_hess_spLT, R_comp1b_hess_spLT)   
})


    
