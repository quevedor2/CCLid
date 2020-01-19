#' checkMse
#' @description Finds the MSE for a NMF fit given a rank of 0 to 2.
#' If the MSE is lowest for rank 0, then this suggests that the 
#' design matrix is complete and perfectly explains the observed data.
#' If the MSE for a rank 1 or rank 2 is lowest, this suggests that
#' the design matrix partially explains the observed data.
#' 
#' Upon examining the H matrix, we can see the VERY rough proportion 
#' of the data in rank k that is unexplained from the design matrix.
#'
#' @param M Input matrix
#' @param M0 Design matrix
#'
#' @return
#' @export
checkMse <- function(M, M0){
  z.k <- sapply(setNames(0:2, 0:2), function(k){
    z <- nnmf(M, k = k, check.k = FALSE, init=list(W0 = M0));
    tail(z$mse,1)
  })
  return(z.k)
}
