context("dnormTrunc0_log")

## To improve accuracy of the numerical estimate, try the complex
## step method, using pracma::erfz to get complex pnorm

test_that("dnormTrunc0_log",{
  require(numDeriv)

  ## M_LN_SQRT_2PI <- 0.9189385332046727805633
  ## M_LN2 <-   0.693147180559945309417232121458


  ## pnorm_complex <- function(x, m, s)  0.5 * (1 + pracma::erfz((x - m) / (s * sqrt(2))))
  ## log_dnorm_complex <- function(x, m, s) -log(s * sqrt(2 * pi)) - 0.5 * ((x - m) / s)^2

  fn <- function(z, m, s) {
    res <- dnorm(z, m, s, log=TRUE) - pnorm(0, m, s, log=TRUE, lower.tail=FALSE)
    ##    res <- log_dnorm_complex(z, m, s) - log(1.0 - pnorm_complex(0, m, s))
    return(res)
  }

  ## CP_err calls the CppAD erf function, which is not correct
  ## pracma::erf and CppADutils::erf calls the std erf, which is
  ## But to test AD, using the CppAD erf, for now

  ##   fn <- function(z, m, s) {

  ##   ##  R <- pracma::erf(-sqrt(0.5) * m / s)
  ## ##    R <- CP_erf(-sqrt(0.5) * m / s)
  ##     R <- MB_erf(-sqrt(0.5) * m / s)
  ##     D <- -M_LN_SQRT_2PI - log(s) - (z - m)^2 / (2 * s * s)
  ##     res <- D - log(1 - R) + M_LN2
  ##     return(res)

  ##   }

  R_func <- function(x) {
    n <- length(x)/3
    res <- 0
    for (i in 1:n) {
      res <- res + fn(x[3*(i-1)+1], x[3*(i-1)+2], x[3*i])
    }
    return(res^2)
  }


  e1 <- pracma::erf(.5)
  e2 <- CP_erf(.5)
  e3 <- MB_erf(.5)



  x <- c(0.8, 2.3, 4.4, 0.3, 2.0, 2.1, -1, .2, 3)

  R_val <- R_func(x)
  R_grad <- grad(R_func, x)
  R_hess <- hessian(R_func, x)
  R_hess_spLT <- tril(drop0(R_hess, 1e-14))

  c_list <- cppad_dnormTrunc0_log(x)
  c_val <- c_list$val
  c_grad <- c_list$grad
  c_hess_dense <- c_list$hess.dense
  c_hess_sp <- c_list$hess.sp
  c_hess_spLT <- tril(c_list$hess.spLT)

  expect_equal(c_val, R_val, tolerance=1e-8)
  expect_equal(c_grad, R_grad, tolerance=1e-8)
  expect_equal(c_hess_dense, R_hess, tolerance=1e-8)
  expect_equal(c_hess_dense, c_hess_sp)
  expect_equal(c_hess_spLT, tril(drop0(c_hess_sp, 1e-8)))
  expect_equal(c_hess_spLT, R_hess_spLT, tolerance=1e-8)
})


