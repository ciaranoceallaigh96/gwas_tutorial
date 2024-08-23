# GWAS Labs Y3 Advanced Genetics & Cell Biology
## Lab 1 - Data QC and Population Stratification

Today we will be performing the first part of two-part practical lab series on conducting a genome-wide association study (GWAS). In this session we will be performing quality control (QC) on our input genotype data, followed by an invetigation of population stratification using principal component analysis (PCA). 

We will first use data and X and you will be given further data at the end of this session with which to do your own analysis and upload a final lab report. Full instructions for this assignment will be available on Moodle. 

First, we will load our genetic_matrix: 

`genetic_matrix <- read.csv("genetic_matrix.raw", sep="")` #REMINDER SNP IDS NEED TO BE CLEANED!

We should also load our BIM file which has the SNP information:

`bim_file <- read.delim("genetic_matrix.bim", header=FALSE) \n
colnames(bim_file) <- c("CHR", "SNP", "CM", "BP", "A1", "A2")`

