## code to prepare `snp6.dat` dataset goes here
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
