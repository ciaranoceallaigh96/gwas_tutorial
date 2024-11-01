# GWAS Labs Y3 Advanced Genetics & Cell Biology
## Lab 2 - GWAS and PRS
Today we will be performing the second part of two-part practical lab series on conducting a genome-wide association study (GWAS). In the last session, we performed quality control (QC) on our data and computed the princicpal compoents of the genotype matrix that we can use as GWAS covariates. In this session, we will perform the actual GWAS on ~10,000 SNPs and then conduct a polygenic risk score (PRS) that can be used to discriminate between cases and controls. 

All genotype information in these sessions are real, but for privacy purposes, the phenotypes are simulated.

Our tutorial dataset is a case/control cohort of 107 individuals of European ancestries genotyped at ~10,000 SNPs across the genome.

Our phenotype is binary but our encoding should be 0 and 1. 

```
genetic_matrix_8$PHENOTYPE <- genetic_matrix_8$PHENOTYPE - 1
```

Let’s make a list of our SNPs.

```
list_of_snps <- colnames(genetic_matrix_8)[7:length(colnames(genetic_matrix_8))]
 ```

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

```
for (i in seq_along(list_of_snps)) {
single_snp_glm <- glm(PHENOTYPE ~ genetic_matrix_9[[list_of_snps[i]]] + PC1 + PC2, data = genetic_matrix_9, family="binomial")
summary_stats[list_of_snps[i], "EffectSize"] <- coef(summary(single_snp_glm))[,1][2]
summary_stats[list_of_snps[i], "PValue"] <- coef(summary(single_snp_glm))[,4][2]
}
summary_stats$SNP <- rownames(summary_stats)
```

Let's look at our top SNPs:

```
head(summary_stats[order(summary_stats$PValue), ])
```

A negative effect size ( log(OR) ) means the allele is a protective allele rather than a risk allele. 

```
alpha <- 0.05
number_of_tests <- length(list_of_snps)
```

- [ ] **Lab Task 2: Perform a bonferonni adjust on our alpha value and create the bonf_alpha variable (Answer: bonf_alpha <- alpha / number_of_tests)**


Now, filter all SNPs with p-value below this threshold
```
significant_snps <- summary_stats[summary_stats$PValue < bonf_alpha, ]
significant_snps
```

Let’s make a Manhattan plot. We first need to refer to our bim file to get the genomic location of each SNP. We can update our summary stats to also include the chromosomal location of each SNP. 

```
head(bim_file)
summary_stats_with_loc <- merge(summary_stats, bim_file, by.x = "SNP", by.y = "SNP")
head(summary_stats_with_loc)
```

```
manhattan(summary_stats_with_loc, chr="CHR", bp="BP", p="PValue", snp="SNP", main="Manhattan Plot", col=c("blue4", "orange3"), genomewideline=-log10(bonf_alpha), ylim=c(0, -log10(bonf_alpha)+1), cex.axis=0.4, suggestiveline=FALSE)`
```

We can also zoom in on a particular region (Chromosome 3). 

```
manhattan(subset(summary_stats_with_loc, CHR == 3), chr="CHR", bp="BP", p="PValue", snp="SNP", xlim=c(4000000,10000000), genomewideline=-log10(bonf_alpha), annotatePval = bonf_alpha, annotateTop=FALSE, suggestiveline=FALSE)`
```

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
```

Le us now calculate how well our PRS model is discriminating between cases and controls.  
```
prs_model <- glm(PHENOTYPE ~ SCORE, data = genetic_matrix_10, family=binomial)
```

How can we measure the success of a PRS? One metric is called R^2 (coefficient of determination) which ranges between 0 and 1, and can be interpreted as the amount of variation being explained by the model. For binary traits, we use a modified R^2 called Nagelkerke’s R^2 

```
result <- nagelkerke(prs_model)
result$Pseudo.R.squared.for.model.vs.null[3]
```

- [ ] **Lab Task 4: What percentage of the overall phenotypic variation is being explained by the PRS model? (Answer:72%))**
