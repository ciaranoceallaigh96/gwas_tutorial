# GWAS Labs Y3 Advanced Genetics & Cell Biology
## Lab 1 - Data QC and Population Stratification

Today we will be performing the first part of two-part practical lab series on conducting a genome-wide association study (GWAS). In this session we will be performing quality control (QC) on our input genotype matrix, followed by an investigation of population stratification using principal component analysis (PCA). 

We will first use a tutorial dataset, and you will be given another dataset at the end of this session with which to do your own analysis and upload a final lab report. Full instructions for this assignment will be available on Moodle. All genotype information in these sessions are real, but for privacy purposes, the phenotypes are simulated. 

Our tutorial dataset is a case/control cohort of X individuals of European ancestries genotyped at Y SNPs across the genome. We want to perform a GWAS on the binary disease phenotype and subsequently build a polygenic risk score (PRS) that we will test out on a Tuscan TSI test cohort (Toscani in Italia from the 1000 Genomes Project.)

The input file "genetic_matrix.raw" has a row for every individual in the dataset. The first six columns contain information on the individual (including phenotype) and the remaining columns contain the minor allele dosage of each SNP (0, 1, or 2). 

First we must define our functions that we will be using:

```
# Function to check, install, and load packages
check_and_install_packages <- function(packages) {
  # Loop over the list of packages
  for (pkg in packages) {
    # Check if the package is installed
    if (!require(pkg, character.only = TRUE)) {
      # If not installed, install the package
      install.packages(pkg, dependencies = TRUE)
      # Load the package
      library(pkg, character.only = TRUE)
    } else {
      # If already installed, just load it
      library(pkg, character.only = TRUE)
    }
  }
}

#Function to calculate Hardy-Weinberg Deviation Statistics
calculate_hwe <- function(genotypes) {
    counts <- table(factor(genotypes, levels = c(0, 1, 2)))
    n <- sum(counts)
    p <- (2 * counts[1] + counts[2]) / (2 * n)
    q <- 1 - p
    expected_counts <- c(n * p^2, 2 * n * p * q, n * q^2)
    chi_sq <- sum((counts - expected_counts)^2 / expected_counts)
    p_value <- pchisq(chi_sq, df = 1, lower.tail = FALSE)
    return(p_value)
}

# Function to calculate Minor Allele Frequencies (MAF)
calculate_maf <- function(genotypes) {
  allele_counts <- table(genotypes, useNA = "no")
  if (length(allele_counts) < 2) {
    return(NA)
  }
  p <- (2 * allele_counts[1] + allele_counts[2]) / (2 * sum(allele_counts))
  return(min(p, 1 - p))
}

#Function to replace NA values with the mean column value
replace_na_with_mean <- function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    return(x)
}
```

Now, we will load (and installl if needed) the relevant packages:

```
required_packages <- c("ggplot2", "dplyr", "qqman", "rcompanion")
check_and_install_packages(required_packages)
```

We will now load our genetic_matrix: 

`genetic_matrix_1 <- read.table("genetic_matrix_cleaned.raw", header=TRUE)` #REMINDER SNP IDS NEED TO BE CLEANED!

Let's take a look at the dataset. 

```
dim(genetic_matrix_1)
genetic_matrix_1[1:5,1:10] #Taking a peak
```

We should also load our BIM file which has the SNP information:

```
bim_file <- read.table("genetic_matrix.bim", header=FALSE)
colnames(bim_file) <- c("CHR", "SNP", "CM", "BP", "A1", "A2")
head(bim_file)
```

We can now begin our QC. Let's first calculate the missing genotype proportion per individual:

```
idv_missingness <- rowMeans(is.na(genetic_matrix_1))
missing_individual_df <- data.frame(IID = genetic_matrix_1$IID, Missing_Proportion = idv_missingness)
```
#mean(idv_missingness) or mean(missing_individual_df$Missing_Proportion)
- [ ] **Lab Task 1: Calculate the mean missingness of the individuals in this dataset (Answer: 0.19%)**

Let's take a look at the individuals with the most missingness by ordering the data:

`head(missing_individual_df[order(-missing_individual_df$Missing_Proportion), ])`

We can visualize the ammount of individuals missingness in our data using a histogram:

```
hist(missing_individual_df$Missing_Proportion, main="Histogram of Individual Missingness", xlab="Missing Proportion", ylab="Frequency")
```

We should remove individuals with a missingness of more than 1%:

```
genetic_matrix_2 <- genetic_matrix_1[idv_missingness <= 0.01,]
```

- [ ] **Lab Task 2: Calculate how many individuals were removed? (Answer: 3)**
#compare orginal dim(genetic_matrix_2) to dim(genetic_matrix) OR sum(idv_missingness>0.01)
Let's do a similar QC step except this time in terms of SNPs. We want to calcuate the amount of missingness per SNP:

```
snp_missingness <- colMeans(is.na(genetic_matrix_2[, 7:ncol(genetic_matrix_2)]))
missing_snp_df <- data.frame(SNP= colnames(genetic_matrix_2)[7:ncol(genetic_matrix_2)], Missing_Proportion = snp_missingness)
````

- [ ] **Lab Task 3: Remove SNPs with more than 2% missingness (Answer: The dimensions of your genetic_matrix_3 object should now be 107 X 49945)**
#genetic_matrix_3 <- genetic_matrix_2[,snp_missingness <= 0.02]

We've now done two important QC steps. The next step will be to select only the autosomal SNPs to simplify our analysis. 

```
autosomal_snps <- bim_file[bim_file$CHR >= 1 & bim_file$CHR <= 22, "SNP"]
autosomal_snps <- c(colnames(genetic_matrix_3)[1:6], autosomal_snps) #add the first six columns as well
genetic_matrix_4 <- genetic_matrix_3[, colnames(genetic_matrix_3) %in% autosomal_snps]
```  

We now want to use a filter for minor allele freqeuncy (MAF). We will remove variants below a 5% freqeuncy threshold:

```
maf_values <- apply(genetic_matrix_4[, 7:ncol(genetic_matrix_4)], 2, calculate_maf) #The 2 simply indicates that the function should be applied to each column. 
head(maf_values)
hist(maf_values, main="MAF Distribution", xlab = "MAF", ylab = "Minor Allele Frequency", breaks=20)
genetic_matrix_5 <- genetic_matrix_4[, c(1:6, which(maf_values >= 0.05) + 6)]
```

We need to exclude any SNPs that are completely out of Hardy-Weinberg Equilibrium (HWE), as measured by a Chi-squared test. If a SNP is out of equilibrium it could suggest a genotyping error – however for disease cases we should have a less stringent threshold as it possible selection against the disease allele could lead to deviances from HWE. 

- [ ] **Lab Task 4: Create two new objects from genotype_matrix, one for cases and one for controls (Answer: You should have 53 cases and 53 controls.)**
#controls <- genetic_matrix_5[genetic_matrix_5$PHENOTYPE == 1, ]
#cases <- genetic_matrix_5[genetic_matrix_5$PHENOTYPE == 2, ]
Let’s apply a stringent p-value threshold (1e-6) for the control SNPs as we should expect them all to be roughly in HWE.

```
hwe_p_values_controls <- apply(controls[, 7:ncol(controls)], 2, calculate_hwe)
head(hwe_p_values_controls)
genetic_matrix_6 <- genetic_matrix_5[, c(1:6, which(hwe_p_values_controls >= 1e-6) + 6)]
```
Let’s now apply a less stringent p-value threshold for cases and controls together (1e-10)

```
hwe_p_values_all <- apply(genetic_matrix_6[, 7:ncol(genetic_matrix_6)], 2, calculate_hwe)
genetic_matrix_7 <- genetic_matrix_6[, c(1:6, which(hwe_p_values_all >= 1e-10) + 6)]
```

Before, we move onto the princicpal component analysis. Let us take a look at how many SNPs and Individuals were removed during our QC processes. 

- [ ] **Lab Task 5: Calculate how many SNPs and Individuals were removed across all QC steps (Answer: 3 individuals and 509 SNPs removed))**
#dim(genetic_matrix_1) - dim(genetic_matrix_7)

We now need to perform PCA on our genotype matrix. This transformation will find the vectors (PCs) that explain the most variation across the entire matrix. 

The first few PCs tend to capture broad ancestry patterns in the data. 

We use PCA to detect genetic outliers, to ensure our sample is genetically homogenous (as a proxy for environmental homogeneity), and to then use as covariates in our GWAS regression (to further account for any environmental effects). 

We won’t need anything other than the SNP columns for PCA. For simplicity, we will also replace any NA genotype values with the mean value of the corresponsing value. 

`snp_matrix <- apply(genetic_matrix_7[,7:ncol(genetic_matrix_7)], 2, replace_na_with_mean)`

Now let us perform the PCA using the prcomp function. 

