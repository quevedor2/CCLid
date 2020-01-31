library(Rfast)
library(optparse)

## This create loads in ./eacon/SAMPLE/SAMPLE__GenomeWideSNP_6_hg19_processed.RDS
## files to extract the PROBE-LEVEL baf adn lrr values in order to aggregate it
## into a singular matrix

option_list <- list(
  make_option(c("-a", "--analysis"), type="character", default='baf',
              help="Log ratio (lrr) or B-allele frequency (baf) pre-processing", metavar="character"),
  make_option(c("-d", "--dataset"), type="character", default='GNE',
              help="Dataset must be either GDSC or CCLE", metavar="character"),
  make_option(c("-s", "--section"), type="character", default='read',
              help="Must be either 'read', 'matrix', 'dist'", metavar="character")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

###################
#### Variables ####
stopifnot(opt$dataset == 'GDSC' | opt$dataset == 'CCLE' | opt$dataset == 'GNE')
rds.id <- switch(opt$dataset,
                 "GNE"='_processed.RDS',
                 '_GenomeWideSNP_6_hg19_processed.RDS')
pdir <- file.path('/mnt/work1/users/pughlab/projects/cancer_cell_lines', opt$dataset, 'eacon')
out.dir <- "../condensed"
setwd(pdir); dir.create(out.dir)


##############################
#### Read in all the data ####
if(opt$section == 'read'){
  print("Reading in data")
  dat <- lapply(list.files(pdir), function(i){
    if(dir.exists(i)){
      print(i)
      tryCatch({
        dat <- readRDS(file.path(i, paste0(i, rds.id)))
        switch(opt$analysis,
               'baf'=dat$data$Tumor_BAF,
               'lrr'=dat$data$Tumor_LogR,
               stop("--analysis must be either lrr or baf"))
      }, error=function(e) {NULL})
    }
  })
  names(dat) <- list.files(pdir)
  if(any(sapply(dat, is.null))) dat <- dat[-which(sapply(dat, is.null))]
  saveRDS(dat, file = file.path(out.dir, paste0(opt$analysis, "s.rds"))) ## use to be bafs
  
  # BAFs generated from EaCoN which uses rawcopy package: rcnorm
  # > head(dat[[1]])
  # 201T
  # CN_473963   NA
  # CN_473964   NA
  # CN_473965   NA
  # > summary(dat[[1]][,1])
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  # -1.0     0.0     0.5     0.5     1.0     1.9 1307380 
}


#########################################
#### Assemble the data into a matrix ####
if(opt$section == 'matrix'){
  dat <- readRDS(file.path(out.dir, paste0(opt$analysis, "s.rds")))
  
  start <- seq(1, length(dat), by=50)
  end <- c(start[-1], length(dat))
  se.df <- data.frame("start"=start, "end"=end)
  
  rn1 <- rownames(dat[[1]])
  rn.equality <- sapply(dat, function(i) all(rownames(i) == rn1))
  if(all(rn.equality)){
    dat.m <- do.call(cbind, dat)
  } else {
    dat.ms <- apply(se.df, 1, function(i){
      print(paste0(i['start'], '-', i['end']))
      dat.exp <- lapply(dat[c(i['start']:i['end'])], function(b) {
        b$ID <- rownames(b)
        b
      })
      Reduce(function(x,y) merge(x,y,by="ID"), dat.exp)
    })
    save(dat.ms, file=file.path(out.dir, paste0(opt$analysis, "-ms.rda")))
    # > head(dat.ms[[1]])
    # ID 201T 22RV1 23132-87
    # 1 CN_000263   NA    NA       NA
    # 2 CN_000265   NA    NA       NA
    # 3 CN_000274   NA    NA       NA
    dat.m <- Reduce(function(x,y) merge(x,y,by="ID"), dat.ms)
  }
  rm(dat)
  
  
  save(dat.m, file=file.path(out.dir, paste0(opt$analysis, "-m.rda")))
  saveRDS(dat.m, file = file.path(out.dir, paste0(opt$analysis, "-mat.rds")))
}

#######################################
#### Combine datasets and distance ####
if(opt$section == 'dist'){
  datasets <- c('GDSC', 'CCLE')
  dat.m <- lapply(datasets, function(i){
    path <- file.path('/mnt/work1/users/pughlab/projects/cancer_cell_lines/',
                      i, 'condensed',  paste0(opt$analysis, "-mat.rds"))
    readRDS(path)
  })
  names(dat.m) <- datasets
  
  keep.idx <- switch(opt$analysis,
                     lrr=grep("CN", gsub("_.*", "", dat.m[[1]][,1])),
                     baf=grep("SNP", gsub("_.*", "", dat.m[[1]][,1])))
  dat.m <- lapply(dat.m, function(i) {
    i[,-1] <- round(i[,-1],2)
    i[keep.idx,]
  })
  na <- which(apply(dat.m[[1]][,-1], 1, function(x) all(is.na(x))))
  
  dat.r <- if(length(na) > 0){
    Reduce(function(x,y) merge(x,y,by='ID'), lapply(dat.m, function(i) i[-na,]))
  } else {
    Reduce(function(x,y) merge(x,y,by='ID'), lapply(dat.m, function(i) i))
  }
  
  #rm(dat.m); gc()
  save(dat.r, file=file.path('/mnt/work1/users/pughlab/projects/cancer_cell_lines/CCL_paper/baf', 
                                 paste0(opt$analysis, 's-matrix.rda')))
  saveRDS(dat.r, file=file.path('/mnt/work1/users/pughlab/projects/cancer_cell_lines/CCL_paper/baf', 
                                 paste0(opt$analysis, 's-matrix.rds')))
  
  #### Load RDS and Calculate Dist ####
  pdir <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines/CCL_paper/baf'
  dat.r2 <- readRDS(file = file.path(pdir, paste0(opt$analysis, 's-matrix.rds')))
  rownames(dat.r2) <- dat.r2$ID
  dat.r2 <- dat.r2[,-1]
  keep.idx <- switch(opt$analysis,
                     lrr=grep("CN", gsub("_.*", "", rownames(dat.r2))),
                     baf=grep("SNP", gsub("_.*", "", rownames(dat.r2))))
  dat.r2 <- dat.r2[keep.idx,]
  dat.r2[is.na(dat.r2)] <- median(as.matrix(dat.r2), na.rm=T)
  
  d.dat <- dista(t(dat.r2), t(dat.r2))
  rownames(d.dat) <- colnames(d.dat) <- colnames(dat.r2)
  
  save(d.dat, file=file.path(pdir, paste0("dist_", paste0(opt$analysis, ".rda"))))
  saveRDS(d.dat, file=file.path(pdir, paste0("dist_", paste0(opt$analysis, ".rds"))))
}