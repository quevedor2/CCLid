#' trainLogit
#' @description Trains a logistic regression
#' 
#' @param balanced Balanced groups
#' @param ... Extra param for .createFormula()
#' @importFrom stats as.formula
#' @importFrom stats glm
#' @importFrom stats binomial
#' @export
trainLogit <- function(balanced, ...){
  fs <- .createFormula(...)
  models <- lapply(fs, function(f){
    glm(f, family=binomial(link='logit'),data=balanced)
  })
  return(models)
}