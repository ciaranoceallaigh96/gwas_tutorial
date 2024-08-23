# GWAS Labs Y3 Advanced Genetics & Cell Biology
## Lab 1 - Data QC and Population Stratification

Today we will be performing the first part of two-part practical lab series on conducting a genome-wide association study (GWAS). In this session we will be performing quality control (QC) on our input genotype data, followed by an invetigation of population stratification using principal component analysis (PCA). 

We will first use data and X and you will be given further data at the end of this session with which to do your own analysis and upload a final lab report. Full instructions for this assignment will be available on Moodle. 

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
- [ ] **Lab Task 1: Calculate the mean missingness of the individuals in this dataset (Answer: X)**

Let's take a look at the individuals with the most missingness by ordering the data:

`head(missing_individual_df[order(-missing_individual_df$Missing_Proportion), ])`

We can visualize the ammount of individuals missingness in our data using a histogram:

```
hist(missing_individual_df$Missing_Proportion, main="Histogram of Individual Missingness", xlab="Missing Proportion", ylab="Frequency")
```

We should remove individuals with a missingness of more than 2%:

```
genetic_matrix <- genetic_matrix[idv_missingness <= 0.02,]
```

- [ ] **Lab Task 2: Calculate how many individuals were removed? (Answer: X)**

Let's do a similar QC step except this time in terms of SNPs. We want to calcuate the amount of missingness per SNP:

```
snp_missingness <- colMeans(is.na(genetic_matrix[, 7:ncol(genetic_matrix)]))
missing_snp_df <- data.frame(SNP= colnames(genetic_matrix)[7:ncol(genetic_matrix)], Missing_Proportion = snp_missingness)
````

- [ ] **Lab Task 3: Remove SNPs with more than 2% missingness (Answer: The dimensions of your genetic_matrix object should now be A X B)**

We've now done two important QC steps. The next step will be to select only the autosomal SNPs to simplify our analysis. 

```
autosomal_snps <- bim_file[bim_file$CHR >= 1 & bim_file$CHR <= 22, "SNP"]
genetic_matrix <- genetic_matrix[, colnames(genetic_matrix) %in% autosomal_snps]
```  

We now want to use a filter for minor allele freqeuncy (MAF). We will remove variants below a 5% freqeuncy threshold:

```
maf_values <- apply(genetic_matrix[, 7:ncol(raw)], 2, calculate_maf) #The 2 simply indicates that the function should be applied to each column. 
head(maf_values)
hist(maf_values, main="MAF Distribution", xlab = "MAF", ylab = "Minor Allele Frequency", breaks=20)
genetic_matrix <- genetic_matrix[, c(1:6, which(maf_values >= 0.05) + 6)]
```

We need to exclude any SNPs that are completely out of Hardy-Weinberg Equilibrium (HWE), as measured by a Chi-squared test. If a SNP is out of equilibrium it could suggest a genotyping error – however for disease cases we should have a less stringent threshold as it possible selection against the disease allele could lead to deviances from HWE. 

- [ ] **Lab Task 3: Create two new objects from genotype_matrix, one for cases and one for controls (Answer: You should have X cases and Y controls.)**

