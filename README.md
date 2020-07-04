# CCLid
CCLid (Cancer Cell Line identification) is designed as a toolkit to address the lack of a publicly available resource for genotype-based cell line authentication. We developed this resource to allow for genotype-matching of any given cancer cell line to the 1,497 unique cancer cell lines found in the GDSC, CCLE or gCSI datasets. Using the B-allele frequencies (BAFs) for all SNPs found in common between the input and reference datasets, this tool will allow for a genotype matching operation that trains and uses a logistic model to calculate the probability of the best cell line matches.

As BAFs can be viewed as an indirect representation of karyotype (i.e. AB has a 0.5/0.50 BAF,while AAB would have a 0.66/0.33 BAF), we developed this tool to scan for segments of the genome that are significantly different between genotypically “matching” cell lines. This function allows for the inference of karyotypically discordant regions, hence representing what is likely to be genetic drift.

## Getting Started
### Installing
The main R package can be installed using the `install_github` command:
```
devtools::install_github('bhklab/CCLid')
library(CCLid)
```

## Downloading Metadata
CCLid requires three main files that are stored on an [external server](https://zenodo.org/record/3926933):
* ref_baf.[desc/bin/rda] - Contains the raw [Probeset x Sample] BAF matrix for all CCLE files. External datasets (such as GDSC and gCSI) can be appended to this matrix to create a multi-institutional cell line matrix.
* meta.df - Contains the mapping of all cancer cell lines and their SNP array filenames used in GDSC, CCLE, and gCSI to their respective cellosaurus CVCL identifiers.
* snp6.dat - Contains probeset data and genomic locations for all SNPs used in the BAF matrix
* [1000/5000/10000/50000/100000/500000].BAF.rds - Contains pre-computed BAF variance information for genomics bins at the set bin size

The following commands need to be run at the start of each CCLid session. It ensures that all data files are downloaded and are then loaded into the Global Environment. 
```
## Load metadata into Global Env
refdir <- '/path/to/dir'
metadata <- c("meta.df",  "affy.omni", "cin70", "gne.meta", "melt.cells", "snp6.dat", "baf")
sapply(metadata, downloadRefCCL, saveDir=refdir, verbose=FALSE)

## Download and read in BAF matrices
bin.size <- 5e5
ref.dat <- CCLid::loadRef(refdir, 'baf', bin.size=bin.size, 
                          just.var=TRUE, meta.df=meta.df, verbose=verbose)
```

### Running CCLid
More in depth instruction can be found in the vignette, but this toolkit operates using a few simple commands:

  * **compareVCF()** will map your input VCF to the SNPs found in the reference BAF matrix. It will append to the existing reference matrix and subset for just SNPs that are found in both your input VCF and the BAF matrix.
```
vcfFile='/path/to/VCFfile.vcf'
vcf_mat <- compareVcf(vcfFile, 
                      var.dat=ref.dat$var, ref.mat=ref.dat$ref,         
                      max.snps=500,
                      snp6.dat=snp6.dat)
```

  * **checkForConcordance()** will do an all-by-all sample comparison of genotypes. Using cell line metadata stored in meta.df, a logistic regression is trained to identify cell lines with matching (M) annotations and nonmatching (NM). The input sample is then compared to all other cell lines and a probability is assigned using this trained model.
```
pred <- checkForConcordance(vcf_mat, sample="SampleX", meta.df=meta.df)
pred$pred$M
```

 * **bafDrift()** main purpose is to identify continuous segments of BAF that are discordant between two samples. As this analysis benefits from more SNPs to create a higher resolution estimation of genetic drift, it will require you to rerun the compareVCF() function using the maximum number of SNPs and reduced for certain samples to save on memory.

```
samples_to_compare <- c('Sample1', 'Sample2')
vcf_mat <- compareVcf(vcfFile, 
                      var.dat=ref.dat$var, ref.mat=ref.dat$ref,         
                      ids=samples_to_compare,
                      snp6.dat=snp6.dat)
bdf <- bafDrift(vcf_mat, centering='median', snp6.dat=snp6.dat)
```


## Authors

* **Rene Quevedo** - *Initial work* - [Quevedor2](https://github.com/quevedor2)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

