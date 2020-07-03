#' assemblePredDat
#' @description Creates a melted data-structure designed for predictions,
#' much in the same vain as the training dataset: Sample X, sample 2, value
#' 
#' @param D.vals Group values (M/NM)
#' @param known.class Default = FALSE, hardcodes a M/NM for demo data
#'
#' @export
assemblePredDat <- function(D.vals, known.class=FALSE){
  pred <- lapply(names(D.vals), function(ana, known.class=FALSE){
    ana.df <- do.call(rbind, D.vals[[ana]])
    colnames(ana.df)[ncol(ana.df)] <- ana
    if(known.class){
      ana.df$class <- rep(c("M", "NM"), sapply(D.vals[[ana]], nrow))
    }
    return(ana.df)
  }, known.class=known.class)
  
  pred <- Reduce(function(x,y) merge(x,y,by=c('Var1', 'Var2', 'class')), pred)
  return(pred)
}