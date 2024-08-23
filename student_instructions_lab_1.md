# GWAS Labs Y3 Advanced Genetics & Cell Biology
## Lab 1 - Data QC and Population Stratification

Today we will be performing the first part of two-part practical lab series on conducting a genome-wide association study (GWAS). In this session we will be performing quality control (QC) on our input genotype data, followed by an invetigation of population stratification using principal component analysis (PCA). 

We will first use data and X and you will be given further data at the end of this session with which to do your own analysis and upload a final lab report. Full instructions for this assignment will be available on Moodle. 

First, we will load our genetic_matrix: 

`genetic_matrix <- read.table("genetic_matrix.raw", header=TRUE)` #REMINDER SNP IDS NEED TO BE CLEANED!

We should also load our BIM file which has the SNP information:

```
bim_file <- read.table("genetic_matrix.bim", header=FALSE)
colnames(bim_file) <- c("CHR", "SNP", "CM", "BP", "A1", "A2")
```

We can now begin our QC. Let's first calculate the missing genotype proportion per individual:

```
idv_missingness <- rowMeans(is.na(genetic_matrix))
missing_individual_df <- data.frame(Individual = 1:nrow(genetic_matrix), Missing_Proportion = idv_missingness)
```
- [ ] Task 1: How would you view the IDs of the individuals with the most missingness?
