##################
#### snp6.dat ####
## code to prepare `snp6.dat` dataset goes here
## It's a list of Affymetrix SNP6 probeset, hg19 loci, rsIDs, and alleles
setwd("~/git/CCL_authenticator/data-raw/")
snp6.bed.raw <- read.table(unz("GenomeWideSNP_6.hg19.bed.zip", 
                               "GenomeWideSNP_6.hg19.bed"), skip = 1, 
                           sep="\t", check.names=F, header=F, 
                           stringsAsFactors = F)
snp6.bed.raw$Probe <- gsub("^AFFX-", "", snp6.bed.raw$V4) %>%
  gsub("_.*", "", .)
snp6.bed <- snp6.bed.raw[,c(4, 1, 3, 13)]
colnames(snp6.bed) <- c('Probe_Set_ID', 'chr', 'end', 'Probe')


require(GenomicRanges)
affyAnno.df <- read.csv("GenomeWideSNP_6.na35.annot.filt.csv.gz", 
                        sep="\t", header=TRUE, comment.char = "#",
                        stringsAsFactors=FALSE, check.names=FALSE)
affyAnno.df <- affyAnno.df[-grep('---', affyAnno.df$Physical_Position),]
affyAnno.df$Physical_Position <- as.numeric(affyAnno.df$Physical_Position)


affyAnno.df <- merge(affyAnno.df, snp6.bed, by="Probe_Set_ID", all.y=T)
affyAnno.df$Physical_Position <- affyAnno.df$end
affyAnno.df <- affyAnno.df[,-grep("Chromosome", colnames(affyAnno.df))]

affyAnno.gr <- makeGRangesFromDataFrame(affyAnno.df, keep.extra.columns = T,strand.field = '',
                                        start.field = 'Physical_Position', end.field = 'end')
seqlevelsStyle(affyAnno.gr) <- 'UCSC'
affyAnno.gr <- sort(affyAnno.gr)
gr.l <- as.list(split(affyAnno.gr, affyAnno.gr$Probe))

snp6.dat <- GRangesList(append(gr.l, c("All"=affyAnno.gr)))

usethis::use_data(snp6.dat, overwrite = T)

#################
#### meta.df ####
## Cell line name by filename dataframe
library(Biobase)

datasets <- c('GDSC', 'CCLE')
col.ids <- list('CCLE'=c('Cell line primary name', 'SNP arrays'),
                'GDSC'=c('Sample Name', 'cel'))

#pdir <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines/'
#cn.dir <- file.path(pdir, 'rds')
setwd("~/git/CCL_authenticator/data-raw/")
cn.dat <- lapply(datasets, function(i){
  readRDS(file.path('.', paste0(i, "_CN.gene.RDS")))
})
names(cn.dat) <- datasets

## Form metadata
meta <- lapply(datasets, function(i){
  phenoData(cn.dat[[i]])@data[,col.ids[[i]]]
})

## Reduce to dataframe
meta.df <- Reduce(function(x,y) merge(x,y, by.x='Sample Name', 
                                      by.y='Cell line primary name', all=T), meta)
meta.df <- taRifx::remove.factors(meta.df[-which(rowSums(is.na(meta.df)) == 3),])
colnames(meta.df) <- c('ID', datasets)
for(d in datasets){
  meta.df[,d] <- as.character(meta.df[,d])
}


meta.df$tmp <- tolower(gsub("[ -/]", "", meta.df$ID))
## Fix TT, T-T cell lines
tt.idx <- grep("^tt$", meta.df$tmp)
meta.df[grep("^T[-._]T$", meta.df$ID),]$tmp <- 't.t'

dup.ids <- which(table(meta.df$tmp) > 1)
for(each.dup in names(dup.ids)){
  datasets <- c('GDSC', 'CCLE')
  
  idc <- grep(paste0("^", each.dup, "$"), meta.df$tmp)
  ds.ids <- apply(meta.df[idc, datasets], 2, na.omit)
  for(each.ds in names(ds.ids)){
    meta.df[idc,each.ds] <- ds.ids[each.ds]
  }
}
meta.df <- meta.df[-which(duplicated(meta.df$tmp)),]

meta.df[16,]$CCLE <- meta.df[17,]$CCLE ## Fixes 786-0
meta.df[309,]$CCLE <- meta.df[310,]$CCLE ## Fixes G-292_Clone_A141B1
meta.df[529,]$CCLE <- meta.df[528,]$CCLE ## Fixes Ishikawa_Heraklio_02ER
meta.df[1005,]$CCLE <- meta.df[953,]$CCLE ## Fixes NIH:OVCAR-3
meta.df[1049,]$CCLE <- meta.df[1050,]$CCLE ## Fixes NIH:OVCAR-3
meta.df <- meta.df[-c(17, 309, 528, 953, 1050),]
meta.df <- meta.df[, -grep("^tmp$", colnames(meta.df))]

usethis::use_data(meta.df, overwrite = T)

#################
#### ref.dat ####
#pdir <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines/CCL_paper/baf'
pdir <- "~/git/CCL_authenticator/data-raw/"
analysis  <- 'baf'
dat.r <- readRDS(file = file.path(pdir, paste0(analysis, 's-matrix.rds')))
rownames(dat.r) <- dat.r$ID
dat.r <- dat.r[,-1]
keep.idx <- switch(analysis,
                   lrr=grep("CN", gsub("_.*", "", rownames(dat.r))),
                   baf=grep("SNP", gsub("_.*", "", rownames(dat.r))))
dat.r <- dat.r[keep.idx,]
if(any(is.na(dat.r))) dat.r[is.na(dat.r)] <- median(as.matrix(dat.r), na.rm=T)

