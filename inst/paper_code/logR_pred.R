library(ggplot2)
library(scales)
library(reshape2)

library(CCLid)
library(Biobase)
library(taRifx)

l=1:5
datasets <- c('GDSC', 'CCLE')
pdir <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines/'
cn.dir <- file.path(pdir, 'rds')
abc.dir  <- file.path(pdir, 'CCL_paper', 'abc', 'rds')
pset.dir <- file.path(pdir, 'rds')
baf.dir <- file.path(pdir, 'CCL_paper', 'baf')


abc <- readRDS(file.path(abc.dir,'GDSC-CCLE.rds'))

lrr <- round(readRDS(file.path(baf.dir, 'dist_lrr.rds')),2)
baf <- round(readRDS(file.path(baf.dir, 'dist_baf.rds')),2)
cn.dat <- lapply(datasets, function(i){
  readRDS(file.path(cn.dir, paste0(i, "_CN.gene.RDS")))
})
names(cn.dat) <- datasets

col.ids <- list('CCLE'=c('Cell line primary name', 'SNP arrays'),
                'GDSC'=c('Sample Name', 'cel'))
meta <- lapply(datasets, function(i){
  phenoData(cn.dat[[i]])@data[,col.ids[[i]]]
})
meta.df <- Reduce(function(x,y) merge(x,y, by.x='Sample Name', 
                                      by.y='Cell line primary name', all=T), meta)
meta.df <- taRifx::remove.factors(meta.df[-which(rowSums(is.na(meta.df)) == 3),])
colnames(meta.df) <- c('ID', datasets)


######################
#### BAF Analysis ####
stopifnot(all(colnames(lrr) == colnames(baf)))
# new.ids <- assignGrpIDs(baf, meta.df)
new.ids <- sapply(colnames(baf), function(i){
  cidx <- grep(i, meta.df, ignore.case = T)
  cidx <- cidx[length(cidx)]
  ridx <- grep(paste0("^", i, "(.cel)?$"), meta.df[,cidx], ignore.case = T)
       
  if(length(ridx) > 0){
    paste0(colnames(meta.df)[cidx], "_", meta.df[ridx,]$ID)
  } else {
    i
  }
})
colnames(lrr) <- rownames(lrr) <- colnames(baf) <- rownames(baf) <- as.character(new.ids)

meltDf <- function(m){
  diag(m) <- m[upper.tri(m)] <- NA
  melt.m <- melt(m)
  melt.m[-which(is.na(melt.m$value)),]
}

dist.mats <- list("lrr"=lrr, "baf"=baf)
#D.vals <- lapply(dist.mats, splitConcordanceVals, meta.df=meta.df)
D.vals <- lapply(dist.mats, function(dm){
  dr.nm <- dm
  
  m.vals <- list()
  for(i in meta.df$ID){
    idx <- grep(paste0("_", i, "$"), colnames(dr.nm))
    m.val <- meltDf(dr.nm[idx,idx,drop=F])
    dr.nm[idx,idx] <- NA
    m.vals[[i]] <- m.val
  }
  m.vals <- do.call(rbind, m.vals)
  nm.vals <- meltDf(dr.nm)
  
  list("M"=m.vals,
       "NM"=nm.vals)
})


## Visualize matching/nonmatching distributions
plotHist <- function(D){
  nm.vals <- D$NM
  m.vals <- D$M
  
  max.x <- max(nm.vals$value) + quantile(nm.vals$value, 0.05)
  hist(m.vals$value, col=alpha("green",0.5), border='white', 
       breaks=seq(0, max.x, by=10), prob = TRUE,  
       xlab = "euclidean-distance", main = "GDSC-CCLE")
  lines(density(m.vals$value), lwd = 2, col = "darkgreen")
  
  hist(nm.vals$value, col=alpha("grey", 0.50), border='white', 
       breaks=seq(0, max.x, by=10), prob = TRUE, add=T)
  lines(density(nm.vals$value), lwd = 2, col = "black")
}

pdf(file.path(baf.dir, "baf-dist_histogram.pdf"))
par(mfrow=c(2,1), mar=c(2, 4.1, 2, 2.1))
plotHist(D.vals$lrr)
plotHist(D.vals$baf)
dev.off()

## Combine BAF-dist with LRR-dist per sample
m.vals <- Reduce(function(x,y) merge(x,y, by=c("Var1", "Var2")), lapply(D.vals, function(i) i$M))
nm.vals <- Reduce(function(x,y) merge(x,y, by=c("Var1", "Var2")), lapply(D.vals, function(i) i$NM))
colnames(nm.vals) <- colnames(m.vals) <- c("Var1", "Var2", names(D.vals))

## Train logistic regression
set.seed(12)
sample.nm <- nm.vals[sample(1:nrow(nm.vals), nrow(m.vals)),]
balanced <- data.frame("lrr"=c(m.vals$lrr, sample.nm$lrr), 
                       "baf"=c(m.vals$baf, sample.nm$baf),
                       "id.raw"=c(rep("Match", nrow(m.vals)),
                                  rep("Nonmatch", nrow(sample.nm))))
balanced$id <- as.integer(factor(balanced$id.raw)) - 1

comb.model <- glm(id ~ lrr + baf,family=binomial(link='logit'),data=balanced)
baf.model <- glm(id ~ baf,family=binomial(link='logit'),data=balanced)
lrr.model <- glm(id ~ lrr,family=binomial(link='logit'),data=balanced)

models <- list("comb"=comb.model, 
               "baf"=baf.model,
               "lrr"=lrr.model)

## Use Logistic regression to predict genotype identity
mkPredictions <- function(pred, models){
  fits <- sapply(names(models), function(m){
    model <- models[[m]]
    if(m == 'comb') m <- c('lrr', 'baf')
    
    predict(model, newdata=pred[,m,drop=F], type='response')
  })
  
  p.fits <- apply(fits, 2, function(pf){
    pf <- cut(pf, breaks = c(0,0.5,1))
    levels(pf) <- c('M', 'NM')
    as.character(pf)
  })
  colnames(fits) <- paste0(colnames(fits), ".fit")
  colnames(p.fits) <- paste0(colnames(p.fits), ".p.fit")
  
  pred <- do.call(cbind, list(pred, fits, p.fits))
  return(pred)
}

pred <- rbind(m.vals, nm.vals)
pred$class <- c(rep("M", nrow(m.vals)), rep("NM", nrow(nm.vals)))
pred <- mkPredictions(pred, models)


## Descriptive statistics
pred$match <- with(pred, class == comb.p.fit)
pred.class <- split(pred, pred$class)
sapply(pred.class, function(i) table(i$match))
#         M      NM
# FALSE   1    1821
# TRUE  463 2399051

# Show the FALSE (class != predicted fit) for each group
mismatch.ids <- lapply(pred.class, function(i){
  i[which(!i$match),]
})
mismatch.ids$M
head(mismatch.ids$NM)


s.idx <- sample(c(1:nrow(pred)), size = 2*10^5, replace = F)
with(pred[s.idx,], plot(lrr.fit, baf.fit, pch=16, col=alpha("black", 0.1)))
with(pred[s.idx,], plot(lrr, baf, pch=16, col=alpha("black", 0.1)))
ggplot(pred.class$M,aes(x=baf.fit,y=lrr.fit)) + stat_binhex()
ggplot(pred.class$NM,aes(x=baf.fit,y=lrr.fit)) + stat_binhex()


pred.class$NM[which(with(pred.class$NM, class=='NM' & p.fit == 'M')),]
pred.class$NM[order(pred.class$M$value),]
pred.class$M[order(pred.class$M$value),]
cl.ids <- c("EPLC", 'FU97', 'HCC193', 'HuCCT1', 'KP-4', 'KYSE-70',
            'MDA-MB-4361', 'MDA-MV-4681', 'NCI-H211', 'TT2', 'VM-CUB-1')
x <- lapply(cl.ids, function(cl.id){
  idx <- intersect(grep(cl.id, pred$Var1), 
                   grep(cl.id, pred$Var2))
  pred[idx,]
})

## Visualize logistic
pdf(file.path(baf.dir, "baf-dist_logistic.pdf"))
ggplot(pred,aes(x=baf,y=baf.fit)) + stat_binhex()
dev.off()

######################
#### ABC Analysis ####
drug.id <- setNames('Palbociclib', 'PD-0332991')
abc$Erlotinib
