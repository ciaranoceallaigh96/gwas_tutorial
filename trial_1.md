# GWAS Labs Y3 Advanced Genetics & Cell Biology
## Lab 1 - Data QC and Population Stratification

Today we will be performing the first part of two-part practical lab series on conducting a genome-wide association study (GWAS). In this session we will be performing quality control (QC) on our input genotype matrix, followed by an investigation of population stratification using principal component analysis (PCA). 

We will first use a tutorial dataset, and you will be given another dataset at the end of this session with which to do your own analysis and upload a final lab report. Full instructions for this assignment will be available on Moodle. All genotype information in these sessions are real, but for privacy purposes, the phenotypes are simulated. 

Our tutorial dataset is a case/control cohort of X individuals of European ancestries genotyped at Y SNPs across the genome. We want to perform a GWAS on the binary disease phenotype and subsequently build a polygenic risk score (PRS). 

The input file "genetic_matrix_10K_cleaned.raw" has a row for every individual in the dataset. The first six columns contain information on the individual (including phenotype) and the remaining columns contain the minor allele dosage of each SNP (0, 1, or 2). 

Data can be downloaded from https://drive.google.com/drive/folders/1nuv4UdJ7MKDOPRDTLAhj3hZYsO01sUsn?usp=drive_link

If you have downloaded the files into your Downloads folder you should switch to that directory in RSudio as follows:

```
setwd"C:/Users/johndoe/Downloads"
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
genetic_matrix_1 <- read.table("genetic_matrix_10K_cleaned.raw", header=TRUE)` 
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

Let's take a look at the individuals with the most missingness by ordering the data:

```
head(missing_individual_df[order(-missing_individual_df$Missing_Proportion), ])`
```

We can visualize the ammount of individuals missingness in our data using a histogram:

```
hist(missing_individual_df$Missing_Proportion, main="Histogram of Individual Missingness", xlab="Missing Proportion", ylab="Frequency")
```

Can you interpet what this graph is showing you? 

We will remove individuals with a missingness of more than 1%:

```
genetic_matrix_2 <- genetic_matrix_1[idv_missingness <= 0.01,]
```

- [ ] **Lab Task 2: Calculate how many individuals were removed.**

Let's do a similar QC step except this time in terms of SNPs. We want to calcuate the amount of missingness per SNP:

```
snp_missingness <- colMeans(is.na(genetic_matrix_2[, 7:ncol(genetic_matrix_2)]))
missing_snp_df <- data.frame(SNP= colnames(genetic_matrix_2)[7:ncol(genetic_matrix_2)], Missing_Proportion = snp_missingness)
head(missing_snp_df)
```

Do you understand the difference between missing_individual_df and missing_snp_df? 

We can also make a histrogram as before. 

```
hist(missing_snp_df$Missing_Proportion, main="Histogram of SNP Missingness", xlab="Missing Proportion", ylab="Frequency")
```

- [ ] **Lab Task 3: Remove SNPs with more than 2% missingness**

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

We now need to perform PCA on our genotype matrix. This transformation will find the vectors (PCs) that explain the most variation across the entire matrix. 

The first few PCs tend to capture broad ancestry patterns in the data. 

We use PCA to detect genetic outliers, to ensure our sample is genetically homogenous (as a proxy for environmental homogeneity), and to then use as covariates in our futrue GWAS regression (to further account for any environmental effects). 

We won’t need anything other than the SNP columns for PCA. For simplicity, we will also replace any NA genotype values with the mean value of the corresponsing value (using a pre-built function). 

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

- [ ] **Lab Task 7: Do the samples fall broadly where we expect them to in PCA-space?**














# GWAS Labs Y3 Advanced Genetics & Cell Biology
## Lab 2 - GWAS and PRS
Today we will be performing the second part of two-part practical lab series on conducting a genome-wide association study (GWAS). 

Our phenotype is binary but our encoding should be 0 and 1. 

`genetic_matrix_8$PHENOTYPE <- genetic_matrix_8$PHENOTYPE - 1`

Let’s make a list of our SNPs.

`list_of_snps <- colnames(genetic_matrix_8)[7:length(colnames(genetic_matrix_8))]`
 
 We can try regression out on the first SNP. 

```
single_snp_glm <- glm(PHENOTYPE ~ genetic_matrix_8[[list_of_snps[1]]], data = genetic_matrix_8, family="binomial")
summary(single_snp_glm)
```

- [ ] **Lab Task 1: What are the effect size and p-value of the SNP (Answer: 0.19, 0.5)**

Let's redo this analysis with principcal components (PCs) included as covariates to account for enviromental effects and population stratification. We already have the pca_scores object from our previous lab. 

```
genetic_matrix_9 <- merge(genetic_matrix_8, pca_scores, by.x = "IID", by.y = "IID") 
single_snp_glm <- glm(PHENOTYPE ~ genetic_matrix_9[[list_of_snps[1]]] + PC1 + PC2, data = genetic_matrix_9, family="binomial")
summary(single_snp_glm)
```

 Let’s initialize a summary statistics object where we can keep a record of all our effect sizes and p-values. 

```
summary_stats <- data.frame(EffectSize = numeric(), PValue = numeric(), stringsAsFactors = FALSE)
summary_stats[list_of_snps[1], "PValue"] <- coef(summary(single_snp_glm))[,4][2]
summary_stats[list_of_snps[1], "EffectSize"] <- coef(summary(single_snp_glm))[,1][2]
print(summary_stats)
```     



Now we can loop through all the other SNPs and sort the summary statistics by the p-values. 

`for (i in seq_along(list_of_snps)) {
single_snp_glm <- glm(PHENOTYPE ~ genetic_matrix_9[[list_of_snps[i]]] + PC1 + PC2, data = genetic_matrix_9, family="binomial")
summary_stats[list_of_snps[i], "EffectSize"] <- coef(summary(single_snp_glm))[,1][2]
summary_stats[list_of_snps[i], "PValue"] <- coef(summary(single_snp_glm))[,4][2]
}
summary_stats$SNP <- rownames(summary_stats)
`

Let's look at our top SNPs

head(summary_stats[order(summary_stats$PValue), ])

A negative effect size ( log(OR) ) means the allele is a protective allele rather than a risk allele. 

```
alpha <- 0.05
number_of_tests <- length(list_of_snps)
```

- [ ] **Lab Task 2: Perform a bonferonni adjust on our alpha value and create the bonf_alpha variable (Answer: bonf_alpha <- alpha / number_of_tests)**


Now, filter all SNPs with p-value below this threshold
significant_snps <- summary_stats[summary_stats$PValue < bonf_alpha, ]
significant_snps
```

Let’s make a Manhattan plot. We first need to refer to our bim file to get the genomic location of each SNP. We can update our summary stats to also include the chromosomal location of each SNP. 

```
head(bim_file)
summary_stats_with_loc <- merge(summary_stats, bim_file, by.x = "SNP", by.y = "SNP")
head(summary_stats_with_loc)
```

`manhattan(summary_stats_with_loc, chr="CHR", bp="BP", p="PValue", snp="SNP", main="Manhattan Plot", col=c("blue4", "orange3"), genomewideline=-log10(bonf_alpha), ylim=c(0, -log10(bonf_alpha)+1), cex.axis=0.4, suggestiveline=FALSE)`


We can also zoom in

`manhattan(subset(summary_stats_with_loc, CHR == 3), chr="CHR", bp="BP", p="PValue", snp="SNP", xlim=c(4000000,10000000), genomewideline=-log10(bonf_alpha), annotatePval = bonf_alpha, annotateTop=FALSE, suggestiveline=FALSE)`

Let’s build a polygenic risk score. The PRS is simply a weighted sum. For each individual we will sum the effect (as estimated by the GWAS) of each SNP, depending on whether or not they have 0, 1, or 2 copies of the risk allele. 

```
snps_to_include <- summary_stats[summary_stats$PValue < 0.001, ]
prs <- data.frame(SCORE=numeric(nrow(genetic_matrix_8)), IID=character(nrow(genetic_matrix_8)), stringsAsFactors=FALSE)

#Loop through all individuals
for (i in 1:nrow(genetic_matrix_8)) {
snp_values <- genetic_matrix_8[i, 7:ncol(genetic_matrix_8)]
match_indices <- match(names(genetic_matrix_8)[7:ncol(genetic_matrix_8)], snps_to_include$SNP)
effect_sizes <- snps_to_include$EffectSize[match_indices]
prs$SCORE[i] <- sum(effect_sizes * snp_values, na.rm = TRUE)
prs$IID[i] <- genetic_matrix_8$IID[i]
}
```
- [ ] **Lab Task 3: What is the mean and standard deviation of the PRS? (Answer: mean(prs$SCORE) and sd(prs$SCORE))**

We can add our PRS results back on our genotype matrix object.

```
genetic_matrix_10 <- merge(genetic_matrix_9, prs, by.x="IID", by.y ="IID")
```

ggplot(genetic_matrix_10, aes(x = SCORE, fill = as.factor(PHENOTYPE))) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Distribution of PRS by Case (1) and Control (0) Status",
    x = "Polygenic Score (Unitless)",
    y = "Density",
    fill = "Phenotype"
  ) +
  theme_minimal()

prs_model <- glm(PHENOTYPE ~ SCORE, data = genetic_matrix_10, family=binomial)


How can we measure the success of a PRS? One metric is called R^2 (coefficient of determination) which ranges between 0 and 1, and can be interpreted as the amount of variation being explained by the model. For binary traits, we use a modified R^2 called Nagelkerke’s R^2 

```
result <- nagelkerke(prs_model)
result$Pseudo.R.squared.for.model.vs.null[3]
```

- [ ] **Lab Task 4: What percentage of the overall phenotypic variation is being explained by the PRS model? (Answer:72%))**
