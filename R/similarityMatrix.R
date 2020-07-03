#' similarityMatrix
#' @description Creates a similarity/concordance matrix given a
#' SNP x Sample matrix and a method of concordance
#' 
#' @param mat Matrix of BAFs
#' @param method Either 'cor', 'jaccard', or 'euclidean'
#' @importFrom stats dist
#' @importFrom stats cor
#' 
#' @export
similarityMatrix <- function(mat, method){
  if(any(is.na(mat))){
    D <- switch(method,
                "euclidean"=dist(t(mat)),
                stop("Only 'euclidean' is set up for matrices with NA values"))
    D <- as.matrix(D)
  } else {
    D <- switch(method,
                "cor"=cor(mat),
                "jaccard"=philentropy::distance(t(mat), method='jaccard'),
                "euclidean"=Rfast::Dist(t(mat), method='euclidean'))
  }
  
  colnames(D) <- rownames(D) <- colnames(mat)
  return(D)
}