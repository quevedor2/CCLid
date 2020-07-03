#' balanceGrps
#' @description Part of the process needed to balance the groups for a logistic
#' regression. It undersamples the non-matching group to be comparable with the
#' matching group
#' 
#' @param D.vals Contains the matching/nonmatching groups
#'
#' @export
#'
balanceGrps <- function(D.vals){
  m.vals <- Reduce(function(x,y) merge(x,y, by=c("Var1", "Var2")), lapply(D.vals, function(i) i$M))
  nm.vals <- Reduce(function(x,y) merge(x,y, by=c("Var1", "Var2")), lapply(D.vals, function(i) i$NM))
  colnames(nm.vals) <- colnames(m.vals) <- c("Var1", "Var2", names(D.vals))
  
  set.seed(12)
  sample.nm <- nm.vals[sample(1:nrow(nm.vals), nrow(m.vals)),] ## Balance groups
  balanced <- data.frame("id.raw"=c(rep("Match", nrow(m.vals)),
                                    rep("Nonmatch", nrow(sample.nm))))
  for(d.type in names(D.vals)){
    balanced$tmp <- c(m.vals[,d.type], sample.nm[,d.type])
    colnames(balanced)[ncol(balanced)] <- d.type
  }
  balanced$id <- as.integer(factor(balanced$id.raw)) - 1
  return(balanced)
}
