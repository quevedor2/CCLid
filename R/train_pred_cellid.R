#' z-statistic
#'
#' @param dat 
#'
#' @return
.zval <- function(dat, p=FALSE){
  z <- (dat - mean(dat)) / sd(dat)
  if(p) z <- 2*pnorm(-abs(z))
  return(z)
}

#' assignGrpIDs
#' @description Creates a standardized set of group names for 
#' samples using a reference
#' 
#' @param mat 
#' @param meta.df 
#'
#' @return
#' @export
assignGrpIDs <- function(mat, meta.df){
  new.ids <- sapply(colnames(mat), function(i){
    i <- gsub(".[xy]$", "", i)
    cidx <- grep(i, meta.df, ignore.case = T)
    cidx <- cidx[length(cidx)]
    ridx <- grep(paste0("^", i, "(.cel)?$"), meta.df[,cidx], ignore.case = T)
    
    if(length(ridx) > 0){
      paste0(colnames(meta.df)[cidx], "_", meta.df[ridx,]$ID)
    } else {
      i
    }
  })
  return(new.ids)
}

#' splitConcordanceVals
#' @description Melts the data structure into Sample 1, Sample 2, concordance
#' and separates into Match/nonmatch if a reference metadata is given
#' 
#' @param dm 
#' @param meta.df 
#'
#' @return
#' @export
splitConcordanceVals <- function(dm, meta.df){
  dr.nm <- dm
  
  .mPosToVec <- function(x, nrw){
    as.integer(sapply(nrw * (x-1), function(i) i+x))
  }
  
  n.rows <- nrow(dm)
  
  all.idx <- sapply(meta.df$ID, function(i) grep(paste0("_", i, "$"), x=colnames(dr.nm)))
  all.idx.v <- unlist(sapply(all.idx, .mPosToVec, nrw=n.rows))
  
  m.vals <- lapply(all.idx, function(idx){
    .meltDf(dr.nm[idx,idx,drop=FALSE])
  })
  dr.nm[all.idx.v] <- NA
  
  ## Old, slightly slower method:
  # m.vals <- list()
  # for(i in meta.df$ID){
  #   idx <- grep(paste0("_", i, "$"), colnames(dr.nm))
  #   m.val <- .meltDf(dr.nm[idx,idx,drop=FALSE])
  #   dr.nm[idx,idx] <- NA
  #   m.vals[[i]] <- m.val
  # }
  m.vals <- do.call(rbind, m.vals)
  nm.vals <- .meltDf(dr.nm)
  
  return(list("M"=m.vals,
              "NM"=nm.vals))
}

#' plotHist
#' @description Plots the histogram for matching/non-matching dataset
#' @param D 
#'
#' @return
#' @export
plotHist <- function(D){
  nm.vals <- D$NM
  m.vals <- D$M
  
  all.vals <- c(nm.vals$value, m.vals$value)
  max.x <- ceiling(max(all.vals) + abs(quantile(all.vals, 0.05)))
  min.x <- floor(min(all.vals) - abs(quantile(all.vals, 0.05)))
  spc <- (max.x - min.x) / 40
  hist(m.vals$value, col=alpha("green",0.5), border='white', 
       breaks=seq(min.x, max.x, by=spc), prob = TRUE,  
       xlab = "euclidean-distance", main = "GDSC-CCLE")
  lines(density(m.vals$value), lwd = 2, col = "darkgreen")
  
  hist(nm.vals$value, col=alpha("grey", 0.50), border='white', 
       breaks=seq(min.x, max.x, by=spc), prob = TRUE, add=T)
  lines(density(nm.vals$value), lwd = 2, col = "black")
}

#' balanceGrps
#' @description Part of the process needed to balance the groups for a logistic
#' regression. It undersamples the non-matching group to be comparable with the
#' matching group
#' 
#' @param D.vals 
#'
#' @return
#' @export
#'
#' @examples
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

#' similarityMatrix
#' @description Creates a similarity/concordance matrix given a
#' SNP x Sample matrix and a method of concordance
#' 
#' @param mat 
#' @param method Either 'cor', 'jaccard', or 'euclidean'
#'
#' @return
#' @export
similarityMatrix <- function(mat, method){
  require(Rfast)
  require(philentropy)
  D <- switch(method,
              "cor"=cor(mat),
              "jaccard"=philentropy::distance(t(mat), method='jaccard'),
              "euclidean"=Rfast::Dist(t(mat), method='euclidean'))
  colnames(D) <- rownames(D) <- colnames(mat)
  return(D)
}

#' assemblePredDat
#' @description Creates a melted data-structure designed for predictions,
#' much in the same vain as the training dataset: Sample X, sample 2, value
#' 
#' @param D.vals 
#' @param known.class 
#'
#' @return
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

#' trainLogit
#' @description Trains a logistic regression
#' 
#' @param balanced 
#' @param ... 
#'
#' @return
#' @export
trainLogit <- function(balanced, ...){
  fs <- .createFormula(...)
  models <- lapply(fs, function(f){
    glm(f, family=binomial(link='logit'),data=balanced)
  })
  return(models)
}

#' mkPredictions
#' @description Uses a trained model to predict on a pred data structure
#' 
#' @param pred 
#' @param models 
#'
#' @return
#' @export
#'
#' @examples
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
  
  return(pred)
}

#### Private Functions ####
.meltDf <- function(m){
  require(reshape2)
  diag(m) <- m[upper.tri(m)] <- NA
  melt.m <- melt(m)
  melt.m[-which(is.na(melt.m$value)),]
}

.createFormula <- function(predictors){
  lapply(predictors, function(p){
    as.formula(paste0("id ~ ", paste(p, collapse="+")))
  })
}

#' .genNonmatchSnpDiff
#' @description generates a matrix of euclidean distance between probesets of 
#' sample pairs if the sample id's dont' match (NON-MATCH)
#' @param col.len 
#' @param dr.nm 
#'
#' @return
.genNonmatchSnpDiff <- function(col.len, dr.nm){
  ## Initialize a matrix of random pairs
  r.idx <- matrix(sample(1:ncol(dr.nm), col.len*2, replace=TRUE), nrow=2)
  ids.match <- TRUE
  while(ids.match){
    ## Initialize a matrix of random pairs
    id.mat <- matrix(colnames(dr.nm)[r.idx], nrow=2)
    id.mat <- gsub("^.*?_", "", id.mat)
    
    m.idx <- apply(id.mat, 2, function(i) {i[1] == i[2]})
    ids.match <- any(m.idx) # check for matching pairs
    
    ## update matching pairs to new indices
    if(ids.match){
      r.idx[,which(m.idx)] <- sample(1:ncol(dr.nm), sum(m.idx)*2, replace=TRUE)
    }
  }
  
  d <- as.data.frame(apply(r.idx, 2, function(com){
    rowDiffs(as.matrix(dr.nm[,com]))
  }))
  rownames(d) <- rownames(dr.nm)
  colnames(d) <- apply(r.idx,2, function(j) paste(colnames(dr.nm)[j],collapse=':'))
  return(d)
}

#' .genMatchSnpDiff
#' @description generates a matrix of euclidean distance between probesets of 
#' sample pairs if the sample id's MATCH
#' @param dr.nm 
#'
#' @return
.genMatchSnpDiff <- function(dr.nm){
  m.d <- sapply(meta.df$ID, function(i){
    idx <- grep(paste0("_", i, "$"), colnames(dr.nm))
    if(length(idx) > 1){
      d <- as.data.frame(apply(combn(idx,2), 2, function(com){
        rowDiffs(as.matrix(dr.nm[,com]))
      }))
      rownames(d) <- rownames(dr.nm)
      colnames(d) <- apply(combn(idx,2),2, function(j) paste(colnames(dr.nm)[j],collapse=':'))
      return(d)
    }
  })
  null.idx <- sapply(m.d, is.null)
  if(any(null.idx)){
    m.d <- do.call(cbind, m.d[-which(null.idx)])
  } 
  
  return(m.d)
}

#' splitSnpDist
#'
#' @param dm 
#' @param meta.df 
#'
#' @return
#' @export
#'
#' @examples
splitSnpDist <- function(dm, meta.df){
  dr.nm <- dm
  m.d <- abs(.genMatchSnpDiff(dr.nm))
  nm.d <- abs(.genNonmatchSnpDiff(col.len=ncol(m.d), dr.nm))
  pred.mat <- as.data.frame(t(cbind(m.d, nm.d)))
  pred.mat$id <- c(rep("M", ncol(m.d)), rep("NM", ncol(nm.d)))
  return(pred.mat)
}