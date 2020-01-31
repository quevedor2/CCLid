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
rownames(meta.df) <- 1:nrow(meta.df)

meta.df[14,]$CCLE <- meta.df[15,]$CCLE ## Fixes 786-0
meta.df[278,]$CCLE <- meta.df[279,]$CCLE ## Fixes G-292_Clone_A141B1
meta.df[486,]$CCLE <- meta.df[485,]$CCLE ## Fixes Ishikawa_Heraklio_02ER
meta.df[941,]$CCLE <- meta.df[893,]$CCLE ## Fixes NIH:OVCAR-3
meta.df[978,]$CCLE <- meta.df[979,]$CCLE ## Fixes PE/CA-PJ15
meta.df[1046,]$CCLE <- meta.df[1047,]$CCLE ## Fixes RS4-11/ RS4;11
meta.df[886,]$CCLE <- meta.df[1130,]$CCLE ## Fixes NCI-SNU-1/ SNU-1

meta.df <- meta.df[-c(15, 279, 485, 893, 979, 1047, 1130),]
meta.df <- meta.df[, -grep("^tmp$", colnames(meta.df))]
meta.df$ID <- gsub(" ", "-", meta.df$ID)

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



###################
#### affy.omni ####
## Cell line name by filename dataframe

#### Functions ####
getSnp <- function(x, snp){
  if(snp=='A'){
    x %<>%
      gsub('\\[', "", .) %>%
      gsub("/.*", "", .) 
  } else if (snp =='B'){
    x %<>%
      gsub('\\[./', "", .) %>%
      gsub("\\]", "", .) 
  } else {
    stop("Specifiy snp=A or snp=B")
  }
  x
}

syncStrandedSnps <- function(snp, strand_x, strand_y){
  revcomp <- list("A"="T",
                  "T"="A",
                  "C"="G",
                  "G"="C")
  strand.idx <- which(strand_x != strand_y)
  snp[strand.idx] <- sapply(snp[strand.idx], function(x) revcomp[[x]])
  snp
}

setSnpToggleFlag <- function(snp_x1, snp_x2, snp_y1, snp_y2){
  tot.orient <- rep(NA, length(snp_x1))
  
  corr.orient <- (snp_x1 == snp_y1) & (snp_x2 == snp_y2)
  tot.orient[which(corr.orient)] <- 'CORRECT'
  
  flip.orient <- (snp_x1 == snp_y2) & (snp_x2 == snp_y1)
  tot.orient[which(flip.orient)] <- 'FLIP'
  
  tot.orient
}


#### Main  ####
require(dplyr)
require(Biobase)
#setwd("/mnt/work1/users/pughlab/projects/cancer_cell_lines/CCL_paper/reference/annotations/mapping")
setwd("~/git/CCL_authenticator/data-raw")
affy <- read.table(gzfile("GenomeWideSNP_6.na35.annot.trimmed.csv.gz"), 
                   header=TRUE, stringsAsFactors = FALSE, check.names=FALSE)
omni <- read.table(gzfile("HumanOmni2.5-4v1_H.trimmed.csv.gz"), 
                   header=TRUE, stringsAsFactors = FALSE, check.names=FALSE)

affy.id <- apply(affy[,c("Chromosome", "Physical_Position")], 1, function(x)
  paste(x, collapse=","))
omni.id <- apply(omni[,c("Chr", "MapInfo")], 1, function(x)
  paste(x, collapse=","))
omni.id <- gsub(" ", "", omni.id)

affy$mergeid <- affy.id
omni$mergeid <- omni.id

affy.omni <- merge(x=affy, y=omni, by="mergeid", all=TRUE)
affy.omni.full <- affy.omni[which(apply(affy.omni, 1, function(x) !any(is.na(x)))),]

affy.omni.full$SNP_A <- getSnp(affy.omni.full$SNP, 'A')
affy.omni.full$SNP_B <- getSnp(affy.omni.full$SNP, 'B')

affy.omni.full$scSNP_A <- with(affy.omni.full, syncStrandedSnps(SNP_A, Strand, RefStrand))
affy.omni.full$scSNP_B <- with(affy.omni.full, syncStrandedSnps(SNP_B, Strand, RefStrand))

affy.omni.full$flip <- with(affy.omni.full, setSnpToggleFlag(Allele_A, Allele_B, 
                                                             scSNP_A, scSNP_B))
metaData <- data.frame(labelDescription=c(
  "Unique merge ID ([chr],[pos])",
  "Affy6 probe set IDs",
  "Affy6 dbSNP IDs",
  "Affy6 Chromosome",
  "Affy6 Genomic Loci (GRCh37)",
  "Affy6 Strand", 
  "Affy6 A-allele",
  "Affy6 B-allele",
  "Omni2.5quad probe set IDs",
  "Omni2.5quad SNPs",
  "Omni2.5quad Chromosome",
  "Omni2.5quad Genomic Loci (GRCh37)",
  "Omni2.5quad Strand",
  "Omni2.5quad A-allele",
  "Omni2.5quad B-allele",
  "Omni2.5quad strand corrected A-allele (ref=Affy6)",
  "Omni2.5quad strand corrected b-allele (ref=Affy6)",
  "Omni2.5 order of alleles in reference to Affy6"))
affy.omni <- AnnotatedDataFrame(data=affy.omni.full, varMetadata=metaData)
usethis::use_data(affy.omni, overwrite = T)
