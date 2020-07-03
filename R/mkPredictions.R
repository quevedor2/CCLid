#' mkPredictions
#' @description Uses a trained model to predict on a pred data structure
#' 
#' @param pred Prediction data
#' @param models logistic reression model
#' @importFrom stats setNames
#' @importFrom stats predict
#' @importFrom stats p.adjust
#' @export
#'
mkPredictions <- function(pred, models){
  if(is.null(names(models))){
    names(models) <- sapply(models, function(m){
      paste(colnames(m$model)[-1], collapse="_")
    })
  }
  
  fits <- sapply(setNames(names(models),
                          names(models)), function(m){
                            model <- models[[m]]
                            
                            predict(model, newdata=pred[,m,drop=FALSE], type='response')
                          })
  
  p.fits <- apply(fits, 2, function(pf){
    pf <- cut(pf, breaks = c(0,0.5,1))
    levels(pf) <- c('M', 'NM')
    as.character(pf)
  })
  colnames(fits) <- paste0(colnames(fits), ".fit")
  colnames(p.fits) <- paste0(colnames(p.fits), ".p.fit")
  
  pred <- do.call(cbind, list(pred, fits, p.fits))
  
  pred$z <- .zval(pred$baf)
  pred$p <- .zval(pred$baf, p=TRUE)
  pred$q <- p.adjust(pred$p, method='fdr')
  pred$Var1 <- as.character(pred$Var1)
  pred$Var2 <- as.character(pred$Var2)
  
  return(pred)
}