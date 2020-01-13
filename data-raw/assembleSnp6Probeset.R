git.dir <- '~/git/CCL_authenticator'
snp6.bed.raw <- read.table(file.path(git.dir, "data-raw", "GenomeWideSNP_6.hg19.bed"), 
                           skip = 1, sep="\t", check.names=F, header=F, stringsAsFactors = F)
snp6.bed <- snp6.bed.raw[,c(4, 1, 2, 6)]
colnames(snp6.bed) <- c('ProbesetID', 'chrom', 'pos', 'strand')
snp6.bed$chrom <- factor(snp6.bed$chrom, levels=paste0("chr", c(1:22, "X", "Y")))
snp6.bed <- snp6.bed[with(snp6.bed, order(chrom, pos)),]

snp6.dat <- split(snp6.bed, f=grepl("^CN", snp6.bed$ProbesetID))
names(snp6.dat) <- c('SNP', 'CN')
snp6.dat$All <- snp6.bed
save(snp6.dat, file=file.path(git.dir, "data", "probesets_snp6.RData"))
