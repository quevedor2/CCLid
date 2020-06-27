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
#devtools::install_github("quevedor2/Rcellosaurus")
#devtools::reload(pkgload::inst('Rcellosaurus'))
library(Rcellosaurus)
library(Biobase)

datasets <- c('GDSC', 'CCLE', 'GNE')
col.ids <- list('CCLE'=c('Cell line primary name', 'SNP arrays'),
                'GDSC'=c('Sample Name', 'cel'),
                'GNE'=c('curated_ID', 'Sample_ID'))

#pdir <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines/'
#cn.dir <- file.path(pdir, 'rds')
## Form metadata
setwd("~/git/CCL_authenticator/data-raw/")
meta <- lapply(datasets, function(i){
  if(i=='GNE'){
    ds <- read.csv(file.path('.', 'GNE_mapping.csv'), header=TRUE)
    pheno <- ds[,col.ids[[i]]]
  } else {
    ds <- readRDS(file.path('.', paste0(i, "_CN.gene.RDS")))
    pheno <- phenoData(ds)@data[,col.ids[[i]]]
  }
  colnames(pheno)[1] <- "cellID"
  return(pheno)
})
names(meta) <- datasets

## Reduce to dataframe
meta.df <- Reduce(function(x,y) merge(x,y, by='cellID', all=T), meta)
meta.df <- taRifx::remove.factors(meta.df[-which(rowSums(is.na(meta.df)) == 4),])
colnames(meta.df) <- c('ID', datasets)
for(d in datasets){ meta.df[,d] <- as.character(meta.df[,d]) }


meta.df$tmp <- tolower(gsub("[ -/]", "", meta.df$ID))
## Fix TT, T-T cell lines
{
  meta.df[grep("^T[-._]T$", meta.df$ID),]$tmp <- 't.t'
  meta.df[grep("^786o$", meta.df$tmp),]$tmp <- '7860'
  meta.df[grep("^nih:ovcar3$", meta.df$tmp),]$tmp <- 'ovcar3'
  meta.df[grep("^pecapj15$", meta.df$tmp),]$tmp <- 'pe_capj15'
  meta.df[grep("^ncisnu1$", meta.df$tmp),]$tmp <- 'snu1'
  meta.df[grep("^a3kawakami$", meta.df$tmp),]$tmp <- 'a3kaw'
  meta.df[grep("^a4fukada$", meta.df$tmp),]$tmp <- 'a4fuk'
  # meta.df[grep("^bxpc3ivcc$", meta.df$tmp),]$tmp <- 'bxpc3'
  meta.df[grep("^cfpac$", meta.df$tmp),]$tmp <- 'cfpac1'
  # meta.df[grep("^colo205ivcc$", meta.df$tmp),]$tmp <- 'colo205'
  meta.df[grep("^colo320hsr$", meta.df$tmp),]$tmp <- 'colo320'
  meta.df[grep("^ef027$", meta.df$tmp),]$tmp <- 'efo27'
  meta.df[grep("^es2to$", meta.df$tmp),]$tmp <- 'es2'
  meta.df[grep("^g292clonea141b1$", meta.df$tmp),]$tmp <- 'g292_clone_a141b1'
  # meta.df[grep("^hcc2157utsw$", meta.df$tmp),]$tmp <- 'hcc2157'
  # meta.df[grep("^hct116ivcc$", meta.df$tmp),]$tmp <- 'hct116'
  # meta.df[grep("^ht29ivcc$", meta.df$tmp),]$tmp <- 'ht29'
  # meta.df[grep("^ht55ivcc$", meta.df$tmp),]$tmp <- 'ht55'
  meta.df[grep("^huh6clone5$", meta.df$tmp),]$tmp <- 'huh6'
  # meta.df[grep("^igrov1ivcc$", meta.df$tmp),]$tmp <- 'igrov1'
  meta.df[grep("^ishikawaheraklio02er$", meta.df$tmp),]$tmp <- 'ishikawa_heraklio_02er'
  meta.df[grep("^jiyoyep2003$", meta.df$tmp),]$tmp <- 'jiyoye'
  meta.df[grep("^jurkatclonee61$", meta.df$tmp),]$tmp <- 'jurkat'
  meta.df[grep("^lc1sqsf$", meta.df$tmp),]$tmp <- 'lc1sq'
  # meta.df[grep("^loximviivcc$", meta.df$tmp),]$tmp <- 'loximvi'
  # meta.df[grep("^mdamb134viatcc$", meta.df$tmp),]$tmp <- 'mdamb134vi'
  meta.df[grep("^mdamb435$", meta.df$tmp),]$tmp <- 'mdamb435s'
  # meta.df[grep("^mdamb468ivcc$", meta.df$tmp),]$tmp <- 'mdamb468'
  # meta.df[grep("^mm1sivcc$", meta.df$tmp),]$tmp <- 'mm1s'
  # meta.df[grep("^ncih1437atcc$", meta.df$tmp),]$tmp <- 'ncih1437'
  # meta.df[grep("^ncih1975ivcc$", meta.df$tmp),]$tmp <- 'ncih1975'
  meta.df[grep("^ncih510$", meta.df$tmp),]$tmp <- 'ncih510a'
  # meta.df[grep("^ncih727atcc$", meta.df$tmp),]$tmp <- 'ncih727'
  # meta.df[grep("^ov90ivcc$", meta.df$tmp),]$tmp <- 'ov90'
  meta.df[grep("^rs4;11$", meta.df$tmp),]$tmp <- 'rs411'
}

dup.ids <- which(table(meta.df$tmp) > 1)
for(each.dup in names(dup.ids)){
  idc <- grep(paste0("^", each.dup, "$"), meta.df$tmp)
  ds.ids <- apply(meta.df[idc, datasets], 2, na.omit)
  for(each.ds in names(ds.ids)){
    if(length(ds.ids[[each.ds]]) > 0) meta.df[idc,each.ds] <- ds.ids[[each.ds]][1]
  }
}
meta.df <- meta.df[-which(duplicated(meta.df$tmp)),]
rownames(meta.df) <- 1:nrow(meta.df)

# id <- 'A141B1'
# meta.df[grep(id, meta.df$ID),]
data(melt.cells)
cvcl.ids <- sapply(meta.df$ID, Rcellosaurus::getCVCL, melt.cells=melt.cells)

# cvcls == 0
{
  nil.idx <- which(sapply(cvcl.ids, length) == 0)
  cvcl.ids[[304]] <- 'CVCL_2047' # Fu-Ov-1
  cvcl.ids[[306]] <- 'CVCL_2909' # G-292_Clone_A141B1
  cvcl.ids[[464]] <- 'CVCL_0326' #   Hep_3B2_1-7   
  cvcl.ids[[513]] <- 'CVCL_0846	' # Hs 688(A).T 
  cvcl.ids[[516]] <- 'CVCL_0851' # Hs 69ST 
  cvcl.ids[[528]] <- 'CVCL_0936' # Hs 832(C).T
  cvcl.ids[[584]] <- 'CVCL_6543' # Ishikawa (Heraklio) 02 ER-
  cvcl.ids[[692]] <- 'CVCL_1338' # KP4
  cvcl.ids[[747]] <- 'CVCL_1381' # LOX-IMIV 
  cvcl.ids[[769]] <- 'CVCL_0620' # MB361.1 
  cvcl.ids[[859]] <- 'CVCL_3041' # NB_TU_1-10  
  cvcl.ids[[1015]] <- 'CVCL_0035' #NCI-PC3 
  cvcl.ids[[1042]] <- 'CVCL_1614' # OAW28NR
  cvcl.ids[[1075]] <- 'CVCL_3935' # OvCA420
  cvcl.ids[[1076]] <- 'CVCL_3936' # OvCA429
  cvcl.ids[[1077]] <- 'CVCL_0475' # OvCA433
  cvcl.ids[[1115]] <- 'CVCL_2679' # PE/CA-PJ34 (clone C12)
  cvcl.ids[[1116]] <- 'CVCL_2680' # PE/CA-PJ41 (clone D2)
  cvcl.ids[[1448]] <- 'CVCL_1779' # UACC-257
}

# cvcls > 1
{
  x <- which(sapply(cvcl.ids, length) > 1) # one cell matches to multiple names
  cvcl.ids[[328]] <- 'CVCL_N736' # G96
  cvcl.ids[[646]] <- 'CVCL_0374' # KG-1
  cvcl.ids[[718]] <- 'CVCL_8147' # L3.3
  cvcl.ids[[1106]] <- 'CVCL_0035' # PC-3
  cvcl.ids[[1426]] <- 'CVCL_1774' # T.T
}

meta.df$CVCL <- as.character(cvcl.ids)
meta.df <- meta.df[, -grep("^tmp$", colnames(meta.df))]
meta.df$ID <- gsub(" ", "-", meta.df$ID)

## Add manually curated pharmacoGX mappings
pgx.map <- read.csv("pharmacogx_mapping.csv", header=TRUE, stringsAsFactors = FALSE)
pgx.map[pgx.map=='#N/A'] <- NA
meta.df <- merge(meta.df, pgx.map, by='ID', all=TRUE)

# After manual revision of meta.df nad merging with rna
pdir <- "~/git/CCL_authenticator/data-raw/"
meta.df <- read.csv(file.path(pdir, "rna_meta_df.csv"), sep=",", header=TRUE,
                    check.names=FALSE, stringsAsFactors = FALSE)
meta.df[meta.df=='#N/A'] <- NA
# 

usethis::use_data(meta.df, overwrite = T)

###############
#### cin70 ####
library("org.Hs.eg.db")

pdir <- "~/git/CCL_authenticator/data-raw/"
cin70 <- read.csv(file = file.path(pdir, "cin70.csv"), header=TRUE, 
                  stringsAsFactors = FALSE, comment.char = "#", skip = 2)
cin70$ENS <- mapIds(org.Hs.eg.db, keys = cin70$Gene, keytype = "SYMBOL", column="ENSEMBL")

usethis::use_data(cin70, overwrite = T)

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

### Functions ###
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


### Main  ###
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

######################
#### rna.meta.dat ####
pdir <- "~/git/CCL_authenticator/data-raw/"
rna.meta.df <- read.csv(file.path(pdir, "rna_meta_df.csv"), sep=",", header=TRUE,
                        check.names=FALSE, stringsAsFactors = FALSE)
rna.meta.df[rna.meta.df=='#N/A'] <- NA
usethis::use_data(rna.meta.df, overwrite = T)

###########################
#### bigmemory ref.mat ####
# library(CCLid)
library (bigmemory)
library (biganalytics)
library (bigtabulate)

PDIR <- "/mnt/work1/users/pughlab/projects/cancer_cell_lines/CCL_paper/CCLid/CCLid"
analysis <- 'baf'
ref.mat <- downloadRefCCL(toupper(analysis), saveDir = PDIR)
saveRDS(ref.mat$ID, file=file.path(PDIR, "ref_mat_ID.rds"))

rownames(ref.mat) <- ref.mat$ID
ref.mat <- as.matrix(ref.mat[,-1] * 100)
storage.mode(ref.mat) <- 'integer'

PDIR <- '/mnt/work1/users/home2/quever/git/CCLid-web/extdata'
setwd(PDIR)
write.table(as.matrix(ref.mat), file=file.path(PDIR, "ref_mat2.csv"), sep=",", 
            row.names=FALSE, col.names=TRUE, quote=FALSE)

ref.matrix <- read.big.matrix(file.path(PDIR, "ref_mat2.csv"), 
                              type ="integer", header = TRUE, 
                              backingfile = paste0("ref_", as.integer(bin.size), ".bin"), 
                              descriptorfile = paste0("ref_", as.integer(bin.size), ".desc"), 
                              extraCols =NULL) 
desc <- describe(ref.matrix)


###################
#### gcsi.meta ####
PDIR='/mnt/work1/users/pughlab/projects/cancer_cell_lines/rnaseq_dat/vcfs/metadata'
delimDIR=file.path(PDIR, "EGAD00001000725/delimited_maps")

gcsi.meta <- read.table(file.path(delimDIR, "Run_Sample_meta_info_fixed.map"), sep=";", header=FALSE, 
                        stringsAsFactors = FALSE, check.names=FALSE, fill=FALSE, comment.char='')
colnames(gcsi.meta) <- gsub("=.*", "", gcsi.meta[1,])
gcsi.meta$Cell_line <- gsub("Cell_line=", "", gcsi.meta$Cell_line)

xml.meta <- read.table(file.path(delimDIR, "xml_meta.tsv"), sep="\t", header=FALSE, 
                       stringsAsFactors = FALSE, check.names=FALSE, fill=FALSE)
ega.meta <- read.table(file.path(delimDIR, "Study_Experiment_Run_sample.map"), sep="\t", header=FALSE, 
                       stringsAsFactors = FALSE, check.names=FALSE, fill=FALSE)
file.meta <- read.table(file.path(PDIR, "gcsi_fileList0725.txt"), header=FALSE, sep="\t",
                        stringsAsFactors = FALSE, check.names=FALSE, fill=FALSE)

m1 <- merge(gcsi.meta, xml.meta, by.x="Cell_line", by.y="V2", all=TRUE)
m2 <- merge(m1, ega.meta, by.x="V1", by.y="V15", all=TRUE)
gcsi.meta <- merge(m2, file.meta, by.x='V1', by.y='V2', all=TRUE)

gcsi.meta <- gcsi.meta[,c(1,2,28,14,24,25,30,6,3)]
saveRDS(gcsi.meta, file=file.path(PDIR, "gcsi_meta.rds"))

gcsi.meta <- gcsi.meta[-which(duplicated(gcsi.meta$Cell_line)),]
gcsi.meta$V1.y <- paste0(gcsi.meta$V1.y, "_1")
write.table(gcsi.meta, file=file.path(PDIR, "gcsi_meta.tsv"), sep="\t",
            col.names=TRUE, row.names=FALSE, quote=FALSE)

#### Generate internal memory ####
# data_dir <- "~/git/CCL_authenticator/data"
# data_files <- list.files(data_dir)
# for( i in data_files){
#   load(file.path(data_dir, i))
# }
# cat(gsub(".rda$", "", data_files), sep=", ")
# usethis::use_data(affy.omni, cin70, gne.meta, melt.cells, 
#                   meta.df, rna.meta.df, snp6.dat, 
#                   internal = TRUE, overwrite = T)

dataraw_dir <- "~/git/CCL_authenticator/data-raw"
ccl_table <- read.csv(file.path(dataraw_dir, "downloadTable.csv"), 
                      check.names = FALSE, 
                      stringsAsFactors = FALSE)
usethis::use_data(ccl_table, internal = FALSE, overwrite = T)