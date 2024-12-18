#Function to clean up SNP names (remove underscore and allele)
clean_snp_ids <- function(ids) {
    sapply(ids, function(id) strsplit(id, "_")[[1]][1])
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


genetic_matrix <- read.csv("//wsl.localhost/Ubuntu-22.04/home/oceallc/GWA_tutorial/1_QC_GWAS/genetic_matrix.raw", sep="")
colnames(genetic_matrix) <- clean_snp_ids(colnames(genetic_matrix))

bim_file <- read.delim("//wsl.localhost/Ubuntu-22.04/home/oceallc/GWA_tutorial/1_QC_GWAS/genetic_matrix.bim", header=FALSE)
colnames(bim_file) <- c("CHR", "SNP", "CM", "BP", "A1", "A2")


# Calculate missing genotype proportion per individual
idv_missingness <- rowMeans(is.na(genetic_matrix))

# Calculate missing genotype proportion per SNP
snp_missingness <- colMeans(is.na(genetic_matrix[, 7:ncol(genetic_matrix)])) # Exclude first 6 non-SNP columns

# Data frames for visualization
missing_individual_df <- data.frame(Individual = 1:nrow(genetic_matrix), Missing_Proportion = idv_missingness)
missing_snp_df <- data.frame(SNP = colnames(genetic_matrix)[7:ncol(genetic_matrix)], Missing_Proportion = snp_missingness)

# Save histograms to PDFs
pdf("histimiss.pdf")
hist(missing_individual_df$Missing_Proportion, main="Histogram of Individual Missingness", xlab="Missing Proportion", ylab="Frequency")
dev.off()

pdf("histlmiss.pdf")
hist(missing_snp_df$Missing_Proportion, main="Histogram of SNP Missingness", xlab="Missing Proportion", ylab="Frequency")
dev.off()

# Print dimensions before filtering
print(dim(genetic_matrix))

# Delete SNPs with missingness >0.2
snp_missingness <- colMeans(is.na(genetic_matrix[, 7:ncol(genetic_matrix)])) 
genetic_matrix <- genetic_matrix[, c(1:6, which(snp_missingness <= 0.2) + 6)]
print(dim(genetic_matrix))
           
# Delete individuals with missingness >0.2
idv_missingness <- rowMeans(is.na(genetic_matrix))
genetic_matrix <- genetic_matrix[idv_missingness <= 0.2, ]
print(dim(genetic_matrix))
           
# Delete SNPs with missingness >0.02
snp_missingness <- colMeans(is.na(genetic_matrix[, 7:ncol(genetic_matrix)]))
genetic_matrix <- genetic_matrix[, c(1:6, which(snp_missingness <= 0.02) + 6)]
print(dim(genetic_matrix))
           
# Delete individuals with missingness >0.02
idv_missingness <- rowMeans(is.na(genetic_matrix))
genetic_matrix <- genetic_matrix[idv_missingness <= 0.02, ]
print(dim(genetic_matrix))

# Print dimensions after filtering
print(dim(genetic_matrix))

###################################################################
### Step 2 ####

# Select autosomal SNPs only
autosomal_snps <- bim_file[bim_file$CHR >= 1 & bim_file$CHR <= 22, "SNP"]
genetic_matrix <- genetic_matrix[, colnames(genetic_matrix) %in% autosomal_snps]





# Calculate MAF for each SNP (excluding the first 6 columns)
maf_values <- apply(genetic_matrix[, 7:ncol(genetic_matrix)], 2, calculate_maf)
maf_df <- data.frame(SNP = colnames(genetic_matrix)[7:ncol(genetic_matrix)], MAF = maf_values)
pdf("MAF_distribution.pdf")
hist(maf_df$MAF, main = "MAF distribution", xlab = "MAF", ylab = "Frequency")
dev.off()
# Filter SNPs with MAF >= 0.05
maf_threshold <- 0.05
genetic_matrix <- genetic_matrix[, c(1:6, which(maf_values >= maf_threshold) + 6)]

# Print dimensions after filtering based on MAF
print(dim(genetic_matrix))


# Separate controls and cases
controls <- genetic_matrix[genetic_matrix$PHENOTYPE == 1, ]
ctrl_snps <- controls[, 7:ncol(controls)]
cases <- genetic_matrix[genetic_matrix$PHENOTYPE == 2, ]
case_snps <- cases[, 7:ncol(cases)]

# Calculate HWE p-values for controls
hwe_p_values_controls <- apply(ctrl_snps, 2, calculate_hwe)

# Filter SNPs with HWE p-value >= 1e-6 in controls
filtered_snp_columns <- which(hwe_p_values_controls >= 1e-6)
controls_filtered <- cbind(controls[, 1:6], ctrl_snps[, filtered_snp_columns])
cases_filtered <- cbind(cases[, 1:6], case_snps[, filtered_snp_columns])

# Combine controls and cases data for the next step
filtered_genetic_matrix <- rbind(
  controls_filtered,
  cases_filtered
)

# Combine cases and controls for the final HWE check, excluding the first six columns
combined_snps <- genetic_matrix[, 7:ncol(genetic_matrix)]

# Calculate HWE p-values for each SNP in the combined dataset
hwe_p_values_combined <- apply(combined_snps, 2, calculate_hwe)

# Filter SNPs with HWE p-value >= 1e-10 in the combined dataset
genetic_matrix <- cbind(genetic_matrix[, 1:6], combined_snps[, which(hwe_p_values_combined >= 1e-10)])

snp_matrix <- as.matrix(genetic_matrix)
snp_matrix <- apply(snp_matrix, 2, replace_na_with_mean)

# Perform PCA
pca_result <- prcomp(snp_matrix, center = TRUE, scale. = TRUE, rank. = 20)
pca_scores <- as.data.frame(pca_result$x)
pca_scores$Individual <- genetic_matrix$Individual

pdf("PCA_plot.pdf")
plot(pca_scores$PC1, pca_scores$PC2, xlab = "PC1", ylab = "PC2", main = "PCA of Genetic Data")
dev.off()

###might have to skip removing relateds and removing heterozygotes and supply our own pruned SNPS
#1349   NA10854            2            1      PROBLEM         0.99 <- has a sex problem


#These people fail due to heterozygotisity:
#1330 "NA12342" 68049 67240 103571 0.02229 0.342972453679119 -3.66711854374478
#1459 "NA12874" 68802 67560 104068 0.0339 0.338874582004074 -5.04839854982741
#delete this person due to relatedness
#13291  NA07045

#These are the non-founders (could be removed manually based on relatedness from FAM file)
13281 NA12344 NA12347 NA12348 1 -9
13291 NA06995 NA07435 NA07037 1 -9
13291 NA06997 NA06986 NA07045 2 -9
13292 NA07014 NA07051 NA07031 2 -9
1330 NA12335 NA12340 NA12341 1 -9
1330 NA12336 NA12342 NA12343 2 -9
1334 NA10846 NA12144 NA12145 1 -9
1334 NA10847 NA12146 NA12239 2 -9
1340 NA07029 NA06994 NA07000 1 -9
1341 NA06991 NA06993 NA06985 2 -9
1344 NA10850 0 NA12058 2 -9
1345 NA07348 NA07357 NA07345 2 -9
1345 NA07349 NA07347 NA07346 1 -9
1346 NA10852 NA12045 0 2 -9
1347 NA10859 NA11881 NA11882 2 -9
1349 NA10853 NA11843 0 1 -9
1350 NA10855 NA11831 NA11832 2 -9
1350 NA10856 NA11829 NA11830 1 -9
1353 NA12375 0 NA12383 1 -9
1353 NA12376 NA12546 NA12489 2 -9
1354 NA12386 NA12399 NA12400 2 -9
1358 NA12707 NA12716 NA12717 1 -9
1358 NA12708 0 NA12718 2 -9
1362 NA10861 NA11994 NA11995 2 -9
1375 NA10863 NA12264 NA12234 2 -9
1377 NA10864 NA11893 NA11894 2 -9
1377 NA10865 NA11891 NA11892 1 -9
1408 NA10830 NA12154 NA12236 1 -9
1408 NA10831 NA12155 NA12156 2 -9
1416 NA10835 NA12248 NA12249 1 -9
1418 NA10836 NA12274 NA12275 2 -9
1418 NA10837 NA12272 NA12273 1 -9
1420 NA10838 NA12003 NA12004 1 -9
1420 NA10839 NA12005 NA12006 2 -9
1421 NA10840 NA12286 NA12287 2 -9
1423 NA10843 NA11919 NA11920 2 -9
1424 NA10845 NA11930 NA11931 1 -9
1444 NA12739 NA12748 NA12749 1 -9
1444 NA12740 NA12750 NA12751 2 -9
1447 NA12752 NA12760 NA12761 1 -9
1447 NA12753 NA12762 NA12763 2 -9
1451 NA12766 NA12775 NA12776 1 -9
1451 NA12767 NA12777 NA12778 2 -9
1454 NA12801 NA12812 NA12813 1 -9
1454 NA12802 NA12814 NA12815 2 -9
1456 NA12817 NA12827 NA12828 1 -9
1456 NA12818 NA12829 NA12830 2 -9
1458 NA12832 NA12842 NA12843 2 -9
1459 NA12864 NA12872 NA12873 1 -9
1459 NA12865 NA12874 NA12875 2 -9
1463 NA12877 NA12889 NA12890 1 -9
1463 NA12878 NA12891 NA12892 2 -9


# Select a single SNP for regression (e.g., the first SNP column)
snp <- genetic_matrix[, 7]

# Run linear regression without covariate
single_snp_lm <- lm(PHENOTYPE ~ snp, data = genetic_matrix)
summary(single_snp_lm)
# Run linear regression with covariate
single_snp_lm_cov <- lm(PHENOTYPE ~ snp + covariate, data = genetic_matrix)
summary(single_snp_lm_cov)


# Function to run linear regression for each SNP
run_regression <- function(snp_column) {
  snp <- genetic_matrix[, snp_column, with=FALSE]
  model <- lm(genetic_matrix$PHENOTYPE ~ snp)
  return(summary(model)$coefficients[2, ])  # Return the SNP coefficient summary
}

# Apply the function to each SNP column
snp_results <- apply(genetic_matrix[, 7:ncol(genetic_matrix)], 2, run_regression)

# Convert results to a data frame
snp_results_df <- data.frame(t(snp_results))
colnames(snp_results_df) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")

# Save results to a file
write.table(snp_results_df, file = "snp_regression_results.txt", quote = FALSE, sep = "\t", row.names = TRUE)


# Function to run linear regression for each SNP with a covariate
run_regression_with_covariate <- function(snp_column) {
  snp <- genetic_matrix[, snp_column, with=FALSE]
  model <- lm(genetic_matrix$PHENOTYPE ~ snp + genetic_matrix[, 2])
  return(summary(model)$coefficients[2, ])  # Return the SNP coefficient summary
}

# Apply the function to each SNP column
snp_results_cov <- apply(genetic_matrix[, 7:ncol(genetic_matrix)], 2, run_regression_with_covariate)

# Convert results to a data frame
snp_results_cov_df <- data.frame(t(snp_results_cov))
colnames(snp_results_cov_df) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")

# Save results to a file
write.table(snp_results_cov_df, file = "snp_regression_results_with_covariate.txt", quote = FALSE, sep = "\t", row.names = TRUE)


snp_results_df$SNP <- rownames(snp_results_df)
# Merge results with bim_file for SNP information
merged_results <- merge(snp_results_df, bim_file, by.x = "SNP", by.y = "SNP")

# Calculate -log10(p-value) for plotting
merged_results$logP <- -log10(merged_results$`Pr(>|t|)`)

manhattan_plot(merged_results, "Manhattan Plot - Without Covariate")







































           





# Load required libraries
library(qqman)

# Read the .fam file
fam_file <- read.table("//wsl.localhost/Ubuntu-22.04/home/oceallc/GWA_tutorial/1_QC_GWAS/genetic_matrix.fam", header=FALSE)
colnames(fam_file) <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")


# Read the genetic matrix
genetic_matrix <- read.csv("//wsl.localhost/Ubuntu-22.04/home/oceallc/GWA_tutorial/1_QC_GWAS/genetic_matrix.raw", sep="")
colnames(genetic_matrix) <- clean_snp_ids(colnames(genetic_matrix))



bim_file <- read.delim("//wsl.localhost/Ubuntu-22.04/home/oceallc/GWA_tutorial/1_QC_GWAS/genetic_matrix.bim", header=FALSE)
colnames(bim_file) <- c("CHR", "SNP", "CM", "BP", "A1", "A2")

# Calculate missing genotype proportion per individual
idv_missingness <- rowMeans(is.na(genetic_matrix))

# Calculate missing genotype proportion per SNP
snp_missingness <- colMeans(is.na(genetic_matrix[, 7:ncol(genetic_matrix)])) # Exclude first 6 non-SNP columns

# Data frames for visualization
missing_individual_df <- data.frame(Individual = genetic_matrix$IID, Missing_Proportion = idv_missingness)
missing_snp_df <- data.frame(SNP = colnames(genetic_matrix)[7:ncol(genetic_matrix)], Missing_Proportion = snp_missingness)

# Save histograms to PDFs
pdf("hapmap_histimiss.pdf")
hist(missing_individual_df$Missing_Proportion, main="Histogram of Individual Missingness", xlab="Missing Proportion", ylab="Frequency")
dev.off()

pdf("hapmap_histlmiss.pdf")
hist(missing_snp_df$Missing_Proportion, main="Histogram of SNP Missingness", xlab="Missing Proportion", ylab="Frequency")
dev.off()

# Print dimensions before filtering
print(dim(genetic_matrix))

# Delete SNPs with missingness >0.2
snp_missingness <- colMeans(is.na(genetic_matrix[, 7:ncol(genetic_matrix)]))
genetic_matrix <- genetic_matrix[, c(1:6, which(snp_missingness <= 0.2) + 6)]
print(dim(genetic_matrix))
           
# Delete individuals with missingness >0.2
idv_missingness <- rowMeans(is.na(genetic_matrix))
genetic_matrix <- genetic_matrix[idv_missingness <= 0.2, ]
print(dim(genetic_matrix))
           
# Delete SNPs with missingness >0.02
snp_missingness <- colMeans(is.na(genetic_matrix[, 7:ncol(genetic_matrix)]))
genetic_matrix <- genetic_matrix[, c(1:6, which(snp_missingness <= 0.02) + 6)]
print(dim(genetic_matrix))
           
# Delete individuals with missingness >0.02
idv_missingness <- rowMeans(is.na(genetic_matrix))
genetic_matrix <- genetic_matrix[idv_missingness <= 0.02, ]
print(dim(genetic_matrix))

# Print dimensions after filtering
print(dim(genetic_matrix))

###################################################################
### Step 2 ####

# Select autosomal SNPs only
autosomal_snps <- bim_file[bim_file$CHR >= 1 & bim_file$CHR <= 22, "SNP"]

# Separate the first six columns
first_six_columns <- genetic_matrix[, 1:6]

# Select columns that match the autosomal SNPs, excluding the first six columns
autosomal_snp_columns <- genetic_matrix[, (colnames(genetic_matrix) %in% autosomal_snps) & (seq_along(colnames(genetic_matrix)) > 6)]

# Combine the first six columns with the selected SNP columns
genetic_matrix <- cbind(first_six_columns, autosomal_snp_columns)


# Calculate MAF for each SNP (excluding the first 7 columns)
maf_values <- apply(genetic_matrix[, 7:ncol(genetic_matrix)], 2, calculate_maf)
maf_df <- data.frame(SNP = colnames(genetic_matrix)[7:ncol(genetic_matrix)], MAF = maf_values)
pdf("hapmap_MAF_distribution.pdf")
hist(maf_df$MAF, main = "MAF distribution", xlab = "MAF", ylab = "Frequency")
dev.off()
# Filter SNPs with MAF >= 0.05
maf_threshold <- 0.05
# Retain the first six columns and filter based on MAF
# Retain the first six columns and filter based on MAF
genetic_matrix <- genetic_matrix[, c(1:6, which(maf_values >= maf_threshold) + 6)]

# Print dimensions after filtering based on MAF
print(dim(genetic_matrix))

# Separate controls and cases
controls <- genetic_matrix[genetic_matrix$PHENOTYPE == 1, ]
ctrl_snps <- controls[, 7:ncol(controls)]
cases <- genetic_matrix[genetic_matrix$PHENOTYPE == 2, ]
case_snps <- cases[, 7:ncol(cases)]

# Calculate HWE p-values for controls
hwe_p_values_controls <- apply(ctrl_snps, 2, calculate_hwe)

# Filter SNPs with HWE p-value >= 1e-6 in controls
filtered_snp_columns <- which(hwe_p_values_controls >= 1e-6)
controls_filtered <- cbind(controls[, 1:6], ctrl_snps[, filtered_snp_columns])
cases_filtered <- cbind(cases[, 1:6], case_snps[, filtered_snp_columns])

# Combine controls and cases data for the next step
genetic_matrix_filtered <- rbind(
  controls_filtered,
  cases_filtered
)

# Combine cases and controls for the final HWE check, excluding the first six columns
combined_snps <- genetic_matrix_filtered[, 7:ncol(genetic_matrix_filtered)]

# Calculate HWE p-values for each SNP in the combined dataset
hwe_p_values_combined <- apply(combined_snps, 2, calculate_hwe)

# Filter SNPs with HWE p-value >= 1e-10 in the combined dataset
genetic_matrix <- cbind(genetic_matrix_filtered[, 1:6], combined_snps[, which(hwe_p_values_combined >= 1e-10)])

# Print dimensions after final filtering based on HWE
print(dim(genetic_matrix))



############################################
#1000G
####################################################
# Read the .fam file
global_fam_file <- read.table("//wsl.localhost/Ubuntu-22.04/home/oceallc/GWA_tutorial/2_Population_stratification/1000G_tutorial_data.fam", header=FALSE)
colnames(fam_file) <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")


# Read the genetic matrix
global_matrix <- read.csv("//wsl.localhost/Ubuntu-22.04/home/oceallc/GWA_tutorial/2_Population_stratification/1000G_tutorial_data.raw", sep="")
colnames(global_matrix) <- clean_snp_ids(colnames(global_matrix))



global_bim_file <- read.delim("//wsl.localhost/Ubuntu-22.04/home/oceallc/GWA_tutorial/2_Population_stratification/1000G_tutorial_data.bim", header=FALSE)
colnames(global_bim_file) <- c("CHR", "SNP", "CM", "BP", "A1", "A2")

# Calculate missing genotype proportion per individual
idv_missingness <- rowMeans(is.na(global_matrix))

# Calculate missing genotype proportion per SNP
snp_missingness <- colMeans(is.na(global_matrix[, 7:ncol(global_matrix)])) # Exclude first 6 non-SNP columns

# Data frames for visualization
missing_individual_df <- data.frame(Individual = global_matrix$IID, Missing_Proportion = idv_missingness)
missing_snp_df <- data.frame(SNP = colnames(global_matrix)[7:ncol(global_matrix)], Missing_Proportion = snp_missingness)

# Save histograms to PDFs
pdf("hapmap_histimiss.pdf")
hist(missing_individual_df$Missing_Proportion, main="Histogram of Individual Missingness", xlab="Missing Proportion", ylab="Frequency")
dev.off()

pdf("hapmap_histlmiss.pdf")
hist(missing_snp_df$Missing_Proportion, main="Histogram of SNP Missingness", xlab="Missing Proportion", ylab="Frequency")
dev.off()

# Print dimensions before filtering
print(dim(global_matrix))

# Delete SNPs with missingness >0.2
snp_missingness <- colMeans(is.na(global_matrix[, 7:ncol(global_matrix)]))
global_matrix <- global_matrix[, c(1:6, which(snp_missingness <= 0.2) + 6)]
print(dim(global_matrix))
           
# Delete individuals with missingness >0.2
idv_missingness <- rowMeans(is.na(global_matrix))
global_matrix <- global_matrix[idv_missingness <= 0.2, ]
print(dim(global_matrix))
           
# Delete SNPs with missingness >0.02
snp_missingness <- colMeans(is.na(global_matrix[, 7:ncol(global_matrix)]))
global_matrix <- global_matrix[, c(1:6, which(snp_missingness <= 0.02) + 6)]
print(dim(global_matrix))
           
# Delete individuals with missingness >0.02
idv_missingness <- rowMeans(is.na(global_matrix))
global_matrix <- global_matrix[idv_missingness <= 0.02, ]
print(dim(global_matrix))

# Print dimensions after filtering
print(dim(global_matrix))


###################################################################
### Step 2 ####

# Select autosomal SNPs only
autosomal_snps <- global_bim_file[global_bim_file$CHR >= 1 & global_bim_file$CHR <= 22, "SNP"]

# Separate the first six columns
first_six_columns <- global_matrix[, 1:6]

# Select columns that match the autosomal SNPs, excluding the first six columns
autosomal_snp_columns <- global_matrix[, (colnames(global_matrix) %in% autosomal_snps) & (seq_along(colnames(global_matrix)) > 6)]

# Combine the first six columns with the selected SNP columns
global_matrix <- cbind(first_six_columns, autosomal_snp_columns)


# Calculate MAF for each SNP (excluding the first 7 columns)
maf_values <- apply(global_matrix[, 7:ncol(global_matrix)], 2, calculate_maf)
maf_df <- data.frame(SNP = colnames(global_matrix)[7:ncol(global_matrix)], MAF = maf_values)
pdf("hapmap_MAF_distribution.pdf")
hist(maf_df$MAF, main = "MAF distribution", xlab = "MAF", ylab = "Frequency")
dev.off()
# Filter SNPs with MAF >= 0.05
maf_threshold <- 0.05
global_matrix <- global_matrix[, c(1:6, which(maf_values >= maf_threshold) + 6)]

# Print dimensions after filtering based on MAF
print(dim(global_matrix))

#No Need to Separate controls and cases

# No need for HWE check


shared_columns <- intersect(colnames(genetic_matrix), colnames(global_matrix))

# Subset both matrices to keep only the shared columns
genetic_matrix_subset <- genetic_matrix[, shared_columns]
global_matrix_subset <- global_matrix[, shared_columns]

# Merge the matrices based on the shared columns
merged_matrix <- merge(genetic_matrix_subset, global_matrix_subset, by = shared_columns, all = TRUE)

# Print the dimensions of the merged matrix
print(dim(merged_matrix))





#########################
# Load the superpopulation file
# Perform PCA on genetic_matrix (excluding the first six columns and the superpopulation column)
snp_matrix <- as.matrix(merged_matrix[, 7:ncol(merged_matrix)])
snp_matrix <- apply(snp_matrix, 2, replace_na_with_mean)
pca_result <- prcomp(snp_matrix, center = TRUE, scale. = TRUE, rank. = 20)
pca_scores <- as.data.frame(pca_result$x)
pca_scores$IID <- merged_matrix$IID

library(ggplot2)

#information on the definitions of these populations can be found here -> https://www.coriell.org/1/NHGRI/Collections/1000-Genomes-Project-Collection/1000-Genomes-Project
colors <- c("ASW" = "red", "CEU" = "blue", "CHB" = "green", "OWN" = "purple", "CHS" = "orange", "FIN" = "black", "GBR" = "pink", "JPT" = "gray", "LWK" = "yellow", "MXL" = "brown", "PUR" = "aquamarine", "TSI" = "darkgreen", "YRI" = "deeppink")
population <- read.table("//wsl.localhost/Ubuntu-22.04/home/oceallc/GWA_tutorial/2_Population_stratification/pop_file.txt", header = TRUE)
pca_scores <- as.data.frame(pca_result$x)
pca_scores$IID <- merged_matrix$IID
pca_scores <- merge(pca_scores, population, by = "IID", all.x = TRUE)

           
pdf("PCA_pop_plot.pdf")
ggplot(pca_scores, aes(x = PC1, y = PC2, color = population)) +
  geom_point(size = 2) +
  scale_color_manual(values = colors, na.translate = FALSE) +
  labs(title = "PCA of Genetic Data", x = "PC1", y = "PC2") +
  theme_minimal() +
  theme(legend.title = element_blank())
dev.off()


# Ensure genetic_matrix is a data frame
genetic_matrix <- as.data.frame(genetic_matrix)

run_regression <- function(snp_column) {
    snp <- as.numeric(snp_column)
    model <- lm(genetic_matrix$PHENOTYPE ~ snp)
    return(summary(model)$coefficients[2, ])  # Return the SNP coefficient summary
}

# Apply the function over the specified columns
snp_results <- lapply(genetic_matrix[, 7:ncol(genetic_matrix)], run_regression)

# Optionally, convert the results to a matrix or data frame if needed
snp_results_df <- do.call(rbind, snp_results)
snp_results_df <- as.data.frame(snp_results_df)
# Ensure the row names of snp_results_df match the SNP column in bim_file
snp_results_df$SNP <- rownames(snp_results_df)

# Merge the two data frames
merged_data <- merge(snp_results_df, bim_file, by.x = "SNP", by.y = "SNP")

# Rename columns to match qqman expectations
colnames(merged_data) <- c("SNP", "Estimate", "Std.Error", "t.value", "P", "CHR", "CM", "BP", "A1", "A2")

# Ensure CHR and BP are numeric
merged_data$CHR <- as.numeric(merged_data$CHR)
merged_data$BP <- as.numeric(merged_data$BP)
merged_data$P <- as.numeric(merged_data$P)

library(qqman)
manhattan(merged_data, chr = "CHR", bp = "BP", p = "P", snp = "SNP", main = "Manhattan Plot", col = c("blue4", "orange3"))  #suggestiveline = -log10(1e-03


