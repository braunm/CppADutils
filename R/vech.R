#' @title vech
#' @param X vector
#' @param M matrix
#' @return vector or matrix
#' @export
vech <- function(M) {
    v <- NULL
    k <- NROW(M)
    stopifnot(k==NCOL(M))
    for (i in 1:k) {
        v <- c(v,M[i:k,i])
    }
    return(v)
}

#' @rdname vech
#' @export
subvech <- function(M) {
    v <- NULL
    k <- NROW(M)
    stopifnot(k==NCOL(M))
    for (i in 1:(k-1)) {
        v <- c(v,M[(i+1):k,i])
    }
    return(v)
}

## #' @rdname vech
## #' @export
## invvech <- function(V) {
##     k <- length(V)
##     n <- (sqrt(8*k+1)-1)/2
##     stopifnot(is.integer(n))
##     M <- matrix(NA,n,n)
##     for (i in 1:n) {
##         st <- sum(1:i)
##         M[i:n,i] <- V[sum(1:i),

##     return(M)
## }

## #' @rdname logit
## #' @export
## inv.logit <- function(x) {
##     exp(log.inv.logit)
## }
