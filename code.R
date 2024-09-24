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

required_packages <- c("ggplot2", "dplyr", "qqman", "rcompanion")
check_and_install_packages(required_packages)

genetic_matrix_1 <- read.table("genetic_matrix_cleaned.raw", header=TRUE)

bim_file <- read.table("genetic_matrix.bim", header=FALSE)
colnames(bim_file) <- c("CHR", "SNP", "CM", "BP", "A1", "A2")
head(bim_file)

idv_missingness <- rowMeans(is.na(genetic_matrix_1))
missing_individual_df <- data.frame(IID = genetic_matrix_1$IID, Missing_Proportion = idv_missingness)

hist(missing_individual_df$Missing_Proportion, main="Histogram of Individual Missingness", xlab="Missing Proportion", ylab="Frequency")

genetic_matrix_2 <- genetic_matrix_1[idv_missingness <= 0.01,]

snp_missingness <- colMeans(is.na(genetic_matrix_2[, 7:ncol(genetic_matrix_2)]))
missing_snp_df <- data.frame(SNP= colnames(genetic_matrix_2)[7:ncol(genetic_matrix_2)], Missing_Proportion = snp_missingness)

genetic_matrix_3 <- genetic_matrix_2[,snp_missingness <= 0.02]

autosomal_snps <- bim_file[bim_file$CHR >= 1 & bim_file$CHR <= 22, "SNP"]
autosomal_snps <- c(colnames(genetic_matrix_3)[1:6], autosomal_snps) #add the first six columns as well
genetic_matrix_4 <- genetic_matrix_3[, colnames(genetic_matrix_3) %in% autosomal_snps]

maf_values <- apply(genetic_matrix_4[, 7:ncol(genetic_matrix_4)], 2, calculate_maf) #The 2 simply indicates that the function should be applied to each column. 
head(maf_values)
hist(maf_values, main="MAF Distribution", xlab = "MAF", ylab = "Minor Allele Frequency", breaks=20)
genetic_matrix_5 <- genetic_matrix_4[, c(1:6, which(maf_values >= 0.05) + 6)]

controls <- genetic_matrix_5[genetic_matrix_5$PHENOTYPE == 1, ]
cases <- genetic_matrix_5[genetic_matrix_5$PHENOTYPE == 2, ] 

hwe_p_values_controls <- apply(controls[, 7:ncol(controls)], 2, calculate_hwe)
head(hwe_p_values_controls)
genetic_matrix_6 <- genetic_matrix_5[, c(1:6, which(hwe_p_values_controls >= 1e-6) + 6)]

hwe_p_values_all <- apply(genetic_matrix_6[, 7:ncol(genetic_matrix_6)], 2, calculate_hwe)
genetic_matrix_7 <- genetic_matrix_6[, c(1:6, which(hwe_p_values_all >= 1e-10) + 6)]

snp_matrix <- apply(genetic_matrix_7[,7:ncol(genetic_matrix_7)], 2, replace_na_with_mean)

pca_result <- prcomp(snp_matrix, center=TRUE, scale.=TRUE,rank.=2)
pca_scores <- as.data.frame(pca_result$x)
pca_scores$IID <- genetic_matrix_7$IID
head(pca_scores)

plot(pca_scores$PC1, pca_scores$PC2, xlab="PC1", ylab="PC2", main="PCA of Genetic Data")

pc1_upper_threshold <- mean(pca_scores$PC1) + 6 * sd(pca_scores$PC1)
pc1_lower_threshold <- mean(pca_scores$PC1) - 6 * sd(pca_scores$PC1)

pc1_outliers <- which(pca_scores$PC1 > pc1_upper_threshold | pca_scores$PC1 < pc1_lower_threshold)
pca_scores[pc1_outliers,]
genetic_matrix_8 <- genetic_matrix_7[-pc1_outliers,]

pc2_upper_threshold <- mean(pca_scores$PC2) + 6 * sd(pca_scores$PC2)
pc2_lower_threshold <- mean(pca_scores$PC2) - 6 * sd(pca_scores$PC2)
pc2_outliers <- which(pca_scores$PC2 > pc2_upper_threshold | pca_scores$PC2 < pc2_lower_threshold)

global_matrix <- read.table("../lab_1_1000G_cleaned.raw", header=TRUE)
dim(global_matrix)
shared_columns <- intersect(colnames(genetic_matrix_8), colnames(global_matrix))

genetic_matrix_8_subset <- genetic_matrix_8[, shared_columns]
global_matrix <- global_matrix[, shared_columns]
merged_matrix <- merge(genetic_matrix_8_subset, global_matrix, by = shared_columns, all = TRUE)
dim(merged_matrix)

population_codes <- read.table("../population_file_lab_1.txt", header=TRUE)
head(population_codes)

merged_snp_matrix <- as.matrix(merged_matrix[, 7:ncol(merged_matrix)])
merged_snp_matrix_2 <- apply(merged_snp_matrix, 2, replace_na_with_mean)
merged_pca_result <- prcomp(merged_snp_matrix_2, center = TRUE, scale. = TRUE, rank. = 2)
merged_pca_scores <- as.data.frame(merged_pca_result$x)
merged_pca_scores$IID <- merged_matrix$IID

merged_pca_scores <- merge(merged_pca_scores, population_codes, by = "IID", all.x = TRUE)
head(merged_pca_scores)
merged_pca_scores[is.na(merged_pca_scores)] <- "my_sample"

colours <- c("EUR" = "red", "ASN" = "blue", "AFR" = "green", "AMR" = "purple", "my_sample" = "orange")

ggplot(merged_pca_scores, aes(x=PC1, y=PC2, color=population)) +
geom_point(size=2) +
scale_colour_manual(values=colours, na.translate=FALSE) +
labs(title="PCA of Genetic Data", x="PC1", y="PC2") +
theme_minimal() +
theme(legend.title=element_blank())

genetic_matrix_8$PHENOTYPE <- genetic_matrix_8$PHENOTYPE - 1
list_of_snps <- colnames(genetic_matrix_8)[7:length(colnames(genetic_matrix_8))]
single_snp_glm <- glm(PHENOTYPE ~ genetic_matrix_8[[list_of_snps[1]]], data = genetic_matrix_8, family="binomial")
summary(single_snp_glm)

genetic_matrix_9 <- merge(genetic_matrix_8, pca_scores, by.x = "IID", by.y = "IID") 
single_snp_glm <- glm(PHENOTYPE ~ genetic_matrix_9[[list_of_snps[1]]] + PC1 + PC2, data = genetic_matrix_9, family="binomial")
summary(single_snp_glm)

summary_stats <- data.frame(EffectSize = numeric(), PValue = numeric(), stringsAsFactors = FALSE)
summary_stats[list_of_snps[1], "PValue"] <- coef(summary(single_snp_glm))[,4][2]
summary_stats[list_of_snps[1], "EffectSize"] <- coef(summary(single_snp_glm))[,1][2]
print(summary_stats)


for (i in seq_along(list_of_snps)) { 
single_snp_glm <- glm(PHENOTYPE ~ genetic_matrix_9[[list_of_snps[i]]] + PC1 + PC2, data = genetic_matrix_9, 
family="binomial") 
summary_stats[list_of_snps[i], "EffectSize"] <- coef(summary(single_snp_glm))[,1][2] 
summary_stats[list_of_snps[i], "PValue"] <- coef(summary(single_snp_glm))[,4][2] 
} 

summary_stats$SNP <- rownames(summary_stats)


alpha <- 0.05
number_of_tests <- length(list_of_snps)
bonf_alpha <- alpha / number_of_tests

summary_stats_with_loc <- merge(summary_stats, bim_file, by.x = "SNP", by.y = "SNP")


manhattan(summary_stats_with_loc, chr="CHR", bp="BP", p="PValue", snp="SNP", main="Manhattan Plot", col=c("blue4", "orange3"), genomewideline=-log10(bonf_alpha), ylim=c(0, -log10(bonf_alpha)+1), cex.axis=0.4, suggestiveline=FALSE)

manhattan(subset(summary_stats_with_loc, CHR == 3), chr="CHR", bp="BP", p="PValue", snp="SNP", xlim=c(7000000,9000000), genomewideline=-log10(bonf_alpha), annotatePval = bonf_alpha, annotateTop=FALSE)

snps_to_include <- summary_stats[summary_stats$PValue < 0.001, ] 
prs <- data.frame(SCORE=numeric(nrow(genetic_matrix_8)), IID=character(nrow(genetic_matrix_8)), stringsAsFactors=FALSE)

for (i in 1:nrow(genetic_matrix_8)) {
snp_values <- genetic_matrix_8[i, 7:ncol(genetic_matrix_8)] 
match_indices <- match(names(genetic_matrix_8)[7:ncol(genetic_matrix_8)], snps_to_include$SNP) 
effect_sizes <- snps_to_include$EffectSize[match_indices] 
prs$SCORE[i] <- sum(effect_sizes * snp_values, na.rm = TRUE) 
prs$IID[i] <- genetic_matrix_8$IID[i] 
}


genetic_matrix_10 <- merge(genetic_matrix_9, prs, by.x="IID", by.y ="IID")


ggplot(genetic_matrix_10, aes(x = SCORE, fill = as.factor(PHENOTYPE))) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Distribution of PRS by Case (1) and Control (0) Status",
    x = "PRS",
    y = "Density",
    fill = "Phenotype"
  ) +
  theme_minimal()

prs_model <- glm(PHENOTYPE ~ SCORE, data = genetic_matrix_10, family=binomial)

result <- nagelkerke(prs_model) 

result$Pseudo.R.squared.for.model.vs.null[3]


