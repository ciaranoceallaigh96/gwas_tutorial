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

#Note need to fix the raw file and the bim file prior to the lab



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
filtered_genetic_matrix <- genetic_matrix[, c(1:6, which(snp_missingness <= 0.2) + 6)]
# Delete individuals with missingness >0.2
filtered_genetic_matrix <- filtered_genetic_matrix[idv_missingness <= 0.2, ]
# Delete SNPs with missingness >0.02
snp_missingness <- colMeans(is.na(filtered_genetic_matrix[, 7:ncol(filtered_genetic_matrix)]))
filtered_genetic_matrix <- filtered_genetic_matrix[, c(1:6, which(snp_missingness <= 0.02) + 6)]
# Delete individuals with missingness >0.02
idv_missingness <- rowMeans(is.na(filtered_genetic_matrix))
filtered_genetic_matrix <- filtered_genetic_matrix[idv_missingness <= 0.02, ]

# Print dimensions after filtering
print(dim(filtered_genetic_matrix))

###################################################################
### Step 2 ####

# Skip sex discrepancy check for now

# Select autosomal SNPs only
autosomal_snps <- bim_file[bim_file$CHR >= 1 & bim_file$CHR <= 22, "SNP"]
genetic_matrix <- filtered_genetic_matrix[, colnames(filtered_genetic_matrix) %in% autosomal_snps]





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
controls <- genetic_matrix[filtered_genetic_matrix$PHENOTYPE == 1, ]
ctrl_snps <- controls[, 7:ncol(controls)]
cases <- genetic_matrix[filtered_genetic_matrix$PHENOTYPE == 2, ]
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
combined_snps <- filtered_genetic_matrix[, 7:ncol(filtered_genetic_matrix)]

# Calculate HWE p-values for each SNP in the combined dataset
hwe_p_values_combined <- apply(combined_snps, 2, calculate_hwe)

# Filter SNPs with HWE p-value >= 1e-10 in the combined dataset
filtered_genetic_matrix <- cbind(filtered_genetic_matrix[, 1:6], combined_snps[, which(hwe_p_values_combined >= 1e-10)])



snp_matrix <- apply(filtered_genetic_matrix, 2, replace_na_with_mean)
pca_result <- prcomp(snp_matrix, center = TRUE, scale. = TRUE, rank=20)
pdf("PCA_plot.pdf")
plot(pca_scores$PC1, pca_scores$PC2, xlab = "PC1", ylab = "PC2", main = "PCA of Genetic Data")
dev.off()

###might have to skip removing relateds and removing heterozygotes and supply our own pruned SNPS
#1349   NA10854            2            1      PROBLEM         0.99 <- has a sex problem


#These people fail due to heterozygotisity:
#1330 "NA12342" 68049 67240 103571 0.02229 0.342972453679119 -3.66711854374478
#1459 "NA12874" 68802 67560 104068 0.0339 0.338874582004074 -5.04839854982741


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
