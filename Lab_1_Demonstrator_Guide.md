# GWAS Labs Y3 Advanced Genetics & Cell Biology
## Lab 1 - Data QC and Population Stratification

Today we will be performing the first part of two-part practical lab series on conducting a genome-wide association study (GWAS). In this session we will be performing quality control (QC) on our input genotype matrix, followed by an investigation of population stratification using principal component analysis (PCA). 

We will first use a tutorial dataset, and you will be given another dataset at the end of this session with which to do your own analysis and upload a final lab report. Full instructions for this assignment will be available on Moodle. All genotype information in these sessions are real, but for privacy purposes, the phenotypes are simulated. 

Our tutorial dataset is a case/control cohort of 110 individuals of European ancestries genotyped at 10002 SNPs across the genome. In the second practical session we will perform a GWAS on the binary disease phenotype and subsequently build a polygenic risk score (PRS) to predict case/control status. 

The input file "genetic_matrix_10K_cleaned.raw" has a row for every individual in the dataset. The first six columns contain information on the individual (including phenotype) and the remaining columns contain the minor allele dosage of each SNP (0, 1, or 2). 

Data can be downloaded from https://drive.google.com/drive/folders/1nuv4UdJ7MKDOPRDTLAhj3hZYsO01sUsn?usp=drive_link

If you have downloaded the files into your Downloads folder you should switch to that directory in RSudio as follows:

```
setwd("C:/Users/johndoe/Downloads")
```

Ask your demonstrator for help if you run into issues getting into the directory your files are saved to. 

First, we must define some functions that we will be using:

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

Now, we will load (and install if needed) the relevant packages:

```
required_packages <- c("ggplot2", "dplyr", "qqman", "rcompanion")
check_and_install_packages(required_packages)
options(scipen=999) #formatting command to prevent the use of E-notation. 
```

Below is the basic workflow we will be following in this lab.<br><br>

<img src="https://github.com/ciaranoceallaigh96/gwas_tutorial/blob/main/QC%20Flowchart%20Graph.png" alt="GWAS QC" width="75%">

<br><br> As a first step, we will load our genetic_matrix: 

```
genetic_matrix_1 <- read.table("genetic_matrix_10K_cleaned.raw", header=TRUE)
```

Let's take a look at the dataset. How many SNPs and individuals does our matrix contain?

```
genetic_matrix_1[1:5,1:10] #Take a peak
dim(genetic_matrix_1)
```

We will also load our BIM file, which has the SNP information:

```
bim_file <- read.table("genetic_matrix_10K.bim", header=FALSE)
colnames(bim_file) <- c("CHR", "SNP", "CM", "BP", "A1", "A2")
head(bim_file)
```

Do you understand how the BIM file relates to the genotype matrix?

We can now begin our QC. Let's first calculate the missing genotype proportion per individual:

```
idv_missingness <- rowMeans(is.na(genetic_matrix_1)) 
missing_individual_df <- data.frame(IID = genetic_matrix_1$IID, Missing_Proportion = idv_missingness)
head(missing_individual_df)
```

- [ ] **Lab Task 1: Calculate the mean missingness of the individuals in this dataset**

(Answer: 0.19%)

mean(idv_missingness)

or

mean(missing_individual_df$Missing_Proportion)



Let's take a look at the individuals with the most missingness by ordering the data:

```
head(missing_individual_df[order(-missing_individual_df$Missing_Proportion), ])
```

We can visualize the ammount of individuals missingness in our data using a histogram:

```
hist(missing_individual_df$Missing_Proportion, main="Histogram of Individual Missingness", xlab="Missing Proportion", ylab="Frequency")
```


<img src="https://github.com/ciaranoceallaigh96/gwas_tutorial/blob/main/hist_idv_missingness.PNG" alt="GWAS QC" width="75%">


Can you interpet what this graph is showing you? 

We will remove individuals with a missingness of more than 1%:

```
genetic_matrix_2 <- genetic_matrix_1[idv_missingness <= 0.01,]
```

- [ ] **Lab Task 2: Calculate how many individuals were removed.**

Answer: 2

compare orginal dim(genetic_matrix_2) to dim(genetic_matrix)

or

sum(idv_missingness>0.01)


Let's do a similar QC step except this time in terms of SNPs. We want to calcuate the amount of missingness per SNP:

```
snp_missingness <- colMeans(is.na(genetic_matrix_2[, 7:ncol(genetic_matrix_2)]))
missing_snp_df <- data.frame(SNP= colnames(genetic_matrix_2)[7:ncol(genetic_matrix_2)], Missing_Proportion = snp_missingness)
head(missing_snp_df)
```

Do you understand the difference between missing_individual_df and missing_snp_df? 

We can also make a histrogram as before. 

```
hist(missing_snp_df$Missing_Proportion, main="Histogram of SNP Missingness", xlab="Missing Proportion", ylab="Number of SNPs")
```

<img src="https://github.com/ciaranoceallaigh96/gwas_tutorial/blob/main/hist_snp_missingness.PNG" alt="GWAS QC" width="75%">



- [ ] **Lab Task 3: Remove SNPs with more than 2% missingness**

Answer: 65 SNPs removed. The dimensions of your genetic_matrix_3 object should now be 108 X 9943. 

To remove: genetic_matrix_3 <- genetic_matrix_2[,snp_missingness <= 0.02]


We've now completed two important QC steps. 

The next step will be to select only the autosomal SNPs to simplify our analysis. 

```
autosomal_snps <- bim_file[bim_file$CHR >= 1 & bim_file$CHR <= 22, "SNP"]
autosomal_snps <- c(colnames(genetic_matrix_3)[1:6], autosomal_snps) #add the first six columns as well
genetic_matrix_4 <- genetic_matrix_3[, colnames(genetic_matrix_3) %in% autosomal_snps]
dim(genetic_matrix_4)
```  

How many non-autosomal SNPs were removed?

We now want to use a filter for minor allele freqeuncy (MAF). We can use a pre-built function called calculate_maf. 

```
maf_values <- apply(genetic_matrix_4[, 7:ncol(genetic_matrix_4)], 2, calculate_maf) #The 2 simply indicates that the function should be applied to each column. 
head(maf_values)
hist(maf_values, main="MAF Distribution", xlab = "MAF", ylab = "Minor Allele Frequency", breaks=20)
```

Can you interpret the resultant histrogram?

We will remove variants below a 5% frequency threshold:

```
genetic_matrix_5 <- genetic_matrix_4[, c(1:6, which(maf_values >= 0.05) + 6)]
```

We now need to exclude any SNPs that are completely out of Hardy-Weinberg Equilibrium (HWE), as measured by a Chi-squared test. If a SNP is out of equilibrium it could suggest a genotyping error. 

For disease cases we will use a less stringent threshold, as it possible selection against the disease allele could lead to deviances from HWE. 

- [ ] **Lab Task 4: Create two new objects from genotype_matrix, one for cases and one for controls. (You may need to consult lecture slides)**

Answer: You should have 53 cases and 54 controls. (one individual has missing phenotype information which is marked by -9 if anyone is aking). 

controls <- genetic_matrix_5[genetic_matrix_5$PHENOTYPE == 1, ]

cases <- genetic_matrix_5[genetic_matrix_5$PHENOTYPE == 2, ] 


Let’s apply a stringent p-value threshold (1e-5) for the control SNPs as we should expect them all to be roughly in HWE.

```
hwe_p_values_controls <- apply(controls[, 7:ncol(controls)], 2, calculate_hwe)
head(hwe_p_values_controls)
genetic_matrix_6 <- genetic_matrix_5[, c(1:6, which(hwe_p_values_controls >= 1e-5) + 6)]
```

How many SNPs were removed?

Let’s now apply a less stringent p-value threshold for cases and controls together (1e-10)

```
hwe_p_values_all <- apply(genetic_matrix_6[, 7:ncol(genetic_matrix_6)], 2, calculate_hwe)
genetic_matrix_7 <- genetic_matrix_6[, c(1:6, which(hwe_p_values_all >= 1e-10) + 6)]
```

How many SNPs were removed?

Before, we move onto the princicpal component analysis. Let us take a look at how many SNPs and Individuals were removed during our QC processes. 

- [ ] **Lab Task 5: Calculate how many SNPs and Individuals were removed across all QC steps**

Answer: 2 individuals and 117 SNPs removed.

dim(genetic_matrix_1) - dim(genetic_matrix_7)


We now need to perform PCA on our genotype matrix. This transformation will find the vectors (PCs) that explain the most variation across the entire matrix. 

The first few PCs tend to capture broad ancestry patterns in the data. 

We use PCA to detect genetic outliers, to ensure our sample is genetically homogenous (as a proxy for environmental homogeneity), and to then use as covariates in our futrue GWAS regression (to further account for any environmental effects). 

This is an image of what a typical PCA plot of mixed populations might look like. Clearly, the individuals do not form one homogenous population.<br><br>

<img src="https://github.com/ciaranoceallaigh96/gwas_tutorial/blob/main/pca_example.png" alt="PCA Example" width="50%">

<br><br>We won’t need anything other than the SNP columns for PCA. For simplicity, we will also replace any NA genotype values with the mean value of the corresponsing value (using a pre-built function). 

```
snp_matrix <- apply(genetic_matrix_7[,7:ncol(genetic_matrix_7)], 2, replace_na_with_mean)
```

Now let us perform the PCA using the prcomp function. We will calculate the first 2 principal components. 

```
pca_result <- prcomp(snp_matrix, center=TRUE, scale.=TRUE,rank.=2)
pca_scores <- as.data.frame(pca_result$x)
pca_scores$IID <- genetic_matrix_7$IID
head(pca_scores)
```

Now we can plot the first two principcal components of variation. 

```
plot(pca_scores$PC1, pca_scores$PC2, xlab="PC1", ylab="PC2", main="PCA of Genetic Data")
```

It looks like we might have three outliers on the PC1 axis. There are also potential outliers on the PC2 axis. 

There are multiple ways to define an outlier - one of which is any point that falls more than six standard deviations from the mean PC value. We will use this to define and remove our outliers here: 

```
pc1_upper_threshold <- mean(pca_scores$PC1) + 6 * sd(pca_scores$PC1)
pc1_lower_threshold <- mean(pca_scores$PC1) - 6 * sd(pca_scores$PC1)

pc1_outliers <- which(pca_scores$PC1 > pc1_upper_threshold | pca_scores$PC1 < pc1_lower_threshold)
pca_scores[pc1_outliers,]
```

Are you surprised by this result? Feel free to discuss with your demonstrator. 

Let's remove the outlier(s). 
```
genetic_matrix_8 <- genetic_matrix_7[-pc1_outliers,]
```

- [ ] **Lab Task 6: Remove all outliers on the PC2 axis, if any.**

Answer: There are no PC2 outliers according to the above defintion. You can discuss how abritrary some of our definitions of outlier can be.

pc2_upper_threshold <- mean(pca_scores$PC2) + 6 * sd(pca_scores$PC2)

pc2_lower_threshold <- mean(pca_scores$PC2) - 6 * sd(pca_scores$PC2) 

pc2_outliers <- which(pca_scores$PC2 > pc2_upper_threshold | pca_scores$PC2 < pc2_lower_threshold)


Let us see where these samples fall on a global sample by using labelled data from the 1000 Genomes project. 

QC has already been done on this data. 

```
global_matrix <- read.table("lab_1_1000G_cleaned.raw", header=TRUE)
dim(global_matrix)
```

Before merging, we must be sure to make sure the two matrices share the same SNPs.

```
shared_columns <- intersect(colnames(genetic_matrix_8), colnames(global_matrix))
```

We will subset both matrices to keep only the shared columns. Then we will merge them together. 

```
genetic_matrix_8_subset <- genetic_matrix_8[, shared_columns]
global_matrix <- global_matrix[, shared_columns]
merged_matrix <- merge(genetic_matrix_8_subset, global_matrix, by = shared_columns, all = TRUE)
dim(merged_matrix)
```

Let us also load in the population codes of the individuals of the 1000 Genomes project. 

They’ve give been given broad superpopulation labels: “AFR=African”, “EUR=European”, “ASN=Asian”,  “AMR=Admixed American”

```
population_codes <- read.table("population_file_lab_1.txt", header=TRUE)
head(population_codes)
```

Let’s redo our PCA with the merged matrix. We will calculate the first two PCS. We will then add back in the population info to each sample. Our own samples will remain NA so we will label it "my_sample". 

```
merged_snp_matrix <- as.matrix(merged_matrix[, 7:ncol(merged_matrix)])
merged_snp_matrix_2 <- apply(merged_snp_matrix, 2, replace_na_with_mean)
merged_pca_result <- prcomp(merged_snp_matrix_2, center = TRUE, scale. = TRUE, rank. = 2)
merged_pca_scores <- as.data.frame(merged_pca_result$x)
merged_pca_scores$IID <- merged_matrix$IID

merged_pca_scores <- merge(merged_pca_scores, population_codes, by = "IID", all.x = TRUE)
head(merged_pca_scores)
merged_pca_scores[is.na(merged_pca_scores)] <- "my_sample"
```

Let us now plot out our PCs from the 1000G samples. 

```
colours <- c("EUR" = "red", "ASN" = "blue", "AFR" = "green", "AMR" = "purple", "my_sample" = "orange")

ggplot(merged_pca_scores, aes(x=PC1, y=PC2, color=population)) +
geom_point(size=2) +
scale_colour_manual(values=colours, na.translate=FALSE) +
labs(title="PCA of Genetic Data", x="PC1", y="PC2") +
theme_minimal() +
theme(legend.title=element_blank())
```

- [ ] **Lab Task 7: Do the samples fall broadly where we expect them to in PC-space?**

Answer: Yes. These European samples fall within the European-associated cluster. It is also not suprising that there is a cline in the admixed individuals, some may have more European ancestry than others. 

