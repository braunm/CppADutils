context("LKJ")

test_that("LKJ", {

    set.seed(123)

    K <- 6
    eta <- 2.1

    cn <- LKJ_const(eta, K)
    cn_R <- lkj_const_R(eta, K)
    expect_equal(cn, cn_R)

    Y <- rnorm(choose(K,2))
    R <- lkj_unwrap_R(Y, K)
    LR <- R$L
    pdfR <- lkj_R(Y, eta, K)

    C <- LKJ_test(Y, eta, K)
    LC <- C$L
    pdfC <- C$logpdf

    expect_equal(LR, LC)
    expect_equal(pdfR, pdfC)

})
