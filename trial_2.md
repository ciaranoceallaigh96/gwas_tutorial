# GWAS Labs Y3 Advanced Genetics & Cell Biology

<img src="https://github.com/att-y3-gcb/gwas_tutorial/blob/main/manhattan.png" alt="GWAS QC" width="75%">


## Lab 2 - GWAS and PRS
Today we will be performing the second part of two-part practical lab series on conducting a genome-wide association study (GWAS). In the last session, we performed quality control (QC) on our data and computed the princicpal compoents of the genotype matrix that we can use as GWAS covariates. In this session, we will perform the actual GWAS on ~10,000 SNPs and then conduct a polygenic risk score (PRS) that can be used to discriminate between cases and controls. 

All genotype information in these sessions are real, but for privacy purposes, the phenotypes are simulated.

Our tutorial dataset is a case/control cohort of 107 individuals of European ancestries genotyped at ~10,000 SNPs across the genome.

Data can be downloaded from [https://drive.google.com/drive/folders/1nuv4UdJ7MKDOPRDTLAhj3hZYsO01sUsn?usp=drive_link](https://drive.google.com/drive/folders/1HIN9Q742nz-Ffd5TFZ7n4nt0PUT6LMXn?usp=drive_link)

If you have downloaded the files into your Downloads folder you should switch to that directory in RStudio as follows:

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
```

Now, we will load (and install if needed) the relevant packages:

```
required_packages <- c("ggplot2", "dplyr", "qqman", "rcompanion")
check_and_install_packages(required_packages)
options(scipen=999) #formatting command to prevent the use of E-notation. 
```

As a first step, we will load our post-QC genetic_matrix (this was our genetic_matrix_8 object from the previous practical.): 

```
genetic_matrix_8 <- read.table("genetic_matrix_8.raw", header=TRUE)
```

Let's take a look at the dataset. How many SNPs and individuals does our matrix contain?

```
genetic_matrix_8[1:5,1:10] #Take a peak
dim(genetic_matrix_8)
```

We will also load our BIM file, which has the SNP information:

```
bim_file <- read.table("genetic_matrix_10K.bim", header=FALSE)
colnames(bim_file) <- c("CHR", "SNP", "CM", "BP", "A1", "A2")
head(bim_file)
```

Hopefully, all this is familiar with you. Feel free to ask if you need clarity on what these files represent. 

We can see the phenotype in our PHENOTYPE column. Our phenotype is binary (2 and 1) but our encoding should be 0 and 1 to make the regression work properly. 

```
genetic_matrix_8$PHENOTYPE <- genetic_matrix_8$PHENOTYPE - 1
```

Let’s now make a list of our SNPs.

```
list_of_snps <- colnames(genetic_matrix_8)[7:length(colnames(genetic_matrix_8))]
 ```

 Remember, a GWAS is just a regression for everyone SNP. We can try regression out on the first SNP. 

```
single_snp_glm <- glm(PHENOTYPE ~ genetic_matrix_8[[list_of_snps[1]]], data = genetic_matrix_8, family="binomial")
summary(single_snp_glm)
```

- [ ] **Lab Task 1: What are the effect size and p-value of the SNP (Answer: 0.19, 0.5)**

<img src="https://github.com/att-y3-gcb/gwas_tutorial/blob/main/summary1.PNG" alt="GWAS QC" width="50%">


The effect from a logistic regresion is reported as the log of the odds ratio, or log(OR). To convert this to a more interpretable number, you can find the exponent. 

```
exp(single_snp_glm$coefficients[2])
```

This gives us the odds ratio. For a rare outcome this can be interpretated as the relative risk; i.e. a carrier of one copy of the risk allele is 1.2 times as likely to be a case rather than a control. 

But, is this the true effect of the SNP on the phenotype? Remember, any ancestry-specific SNP correlated with an enviromental effect might become spuriously associated with the phenotype. We can discount this ancestry-mediated effect by making use of our princicpal components we calcucated in the last lab. Remember, these PCs tend to capture broad ancestry patterns in the genetic data. 


```
pca_scores <- read.table("pca_scores.txt", header=TRUE)
head(pca_scores)
```

We can merge our pca_scores with our genotype matrix.

```
genetic_matrix_9 <- merge(genetic_matrix_8, pca_scores, by.x = "IID", by.y = "IID") 
```

Let's redo this analysis with principcal components (PCs) included as covariates to account for enviromental effects and population stratification. This is a simple modification of the regression equation.  

```
single_snp_glm <- glm(PHENOTYPE ~ genetic_matrix_9[[list_of_snps[1]]] + PC1 + PC2, data = genetic_matrix_9, family="binomial")
summary(single_snp_glm)
```

- [ ] **Lab Task 2: What is the new effect size and p-value of the SNP once ancestry/enviroment is taken into account? (Answer: 0.21, 0.45)**

<img src="https://github.com/att-y3-gcb/gwas_tutorial/blob/main/summary2.PNG" alt="GWAS QC" width="50%">


 Let’s initialize a summary statistics object where we can keep a record of all our effect sizes and p-values. 

```
summary_stats <- data.frame(EffectSize = numeric(), PValue = numeric(), stringsAsFactors = FALSE)
summary_stats[list_of_snps[1], "PValue"] <- coef(summary(single_snp_glm))[,4][2]
summary_stats[list_of_snps[1], "EffectSize"] <- coef(summary(single_snp_glm))[,1][2]
print(summary_stats)
```     


Now we can loop through all the other SNPs and sort the summary statistics by the p-values. This is a far more efficient way to perform the GWAS than taking each SNP one-by-one. 

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

A negative effect size ( log(OR) ) means the minor allele is a protective allele rather than a risk allele. A GWAS is always reported in reference to the minor (less common allele). Another way of seeing this is that the major (common) allele of the SNP confers risk on the carrier. 

Remember, for a normal GWAS the significance p-value theshold is set to 10^-8. This is beacuse you generally test 1,000,000 independant regions of the genome. Here, we are testing only ~9,000 SNPs so our p-value threshold doesn't need to be as stringent. Let's calculate a threshold appropriate for our analysis. 

```
alpha <- 0.05
number_of_tests <- length(list_of_snps)
```

- [ ] **Lab Task 3: Perform a bonferonni adjust on our alpha value and create the bonf_alpha variable (Answer: bonf_alpha <- alpha / number_of_tests)**


Now, let us filter all SNPs with p-value below this threshold:

```
significant_snps <- summary_stats[summary_stats$PValue < bonf_alpha, ]
significant_snps
```

How many significant SNPs do we have?


<img src="https://github.com/att-y3-gcb/gwas_tutorial/blob/main/manhattan.png" alt="GWAS QC" width="75%">


Let’s now make a Manhattan plot. We first need to refer to our bim file to get the genomic location of each SNP. We can update our summary stats to also include the chromosomal location of each SNP. 

```
head(bim_file)
summary_stats_with_loc <- merge(summary_stats, bim_file, by.x = "SNP", by.y = "SNP")
head(summary_stats_with_loc)
```

A manhattan plot shows the genomic location along the x-axis and the p-value (not the effect size!) along the y-axis. We have a package function that can take summary statistics and generate a manhattan plot. 

```
manhattan(summary_stats_with_loc, chr="CHR", bp="BP", p="PValue", snp="SNP", main="Manhattan Plot", col=c("blue4", "orange3"), genomewideline=-log10(bonf_alpha), ylim=c(0, -log10(bonf_alpha)+1), cex.axis=0.4, suggestiveline=FALSE)
```

<img src="https://github.com/att-y3-gcb/gwas_tutorial/blob/main/manhattan3.PNG" alt="GWAS QC" width="50%">


Can you make sense of this graph? Look at the peak(s). 


We can also zoom in on a particular peak to really get a sense of what is going on (Chromosome 3). 

```
manhattan(subset(summary_stats_with_loc, CHR == 3), chr="CHR", bp="BP", p="PValue", snp="SNP", xlim=c(4000000,10000000), genomewideline=-log10(bonf_alpha), annotatePval = bonf_alpha, annotateTop=FALSE, suggestiveline=FALSE)
```

<img src="https://github.com/att-y3-gcb/gwas_tutorial/blob/main/manhattan32.PNG" alt="GWAS QC" width="50%">


Congratulations! You have performed a mini-GWAS from start to finish. You might choose to focus on figuring out which gene(s) are implicated in your phenotype (and how) but this involves a lot of follow-up work. Instead, one can build a polygenic risk score from the results of the GWAS. This is because now we know how much each SNP contributes to the phenotype through the effect size. We can calculate each individuals sum of risk alleles. 

<img src="https://github.com/att-y3-gcb/gwas_tutorial/blob/main/polygenic.PNG" alt="GWAS QC" width="50%">


So let’s build the polygenic risk score. The PRS is simply a weighted sum. For each individual we will sum the effect (as estimated by the GWAS) of each SNP, depending on whether or not they have 0, 1, or 2 copies of the minor allele. 

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

Let's take a look at our calculated PRS.

```
head(prs)
```

- [ ] **Lab Task 4: What is the mean and standard deviation of the PRS? (Answer: mean(prs$SCORE) and sd(prs$SCORE))**

We can add our PRS results back on our genotype matrix object.

```
genetic_matrix_10 <- merge(genetic_matrix_9, prs, by.x="IID", by.y ="IID")
```

The PRS is a little hard to interpret...maybe if we plot out the PRS distribution in cases vs. controls we can get a better feel for how well we can predict case status. 

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

<img src="https://github.com/att-y3-gcb/gwas_tutorial/blob/main/prs.PNG" alt="GWAS QC" width="50%">


There seems to be very good discrimination between cases and controls. Can we put a number on this?

Remember, from our lectures that the maxmimum amount of variation of a PRS can be expected to explain is equal to the heritability of a trait i.e. if 80% of the variation in height is due to genetic factos, then a predictive tool built solely from genetic factors (PRS) cannot exceed 80% of variation explained. For the simulated disease phenotype, the heritability has also been set to 80%. 

SO how do we measure the success of a PRS? One metric is called R^2 (coefficient of determination) which ranges between 0 and 1, and can be interpreted as the amount of variation being explained by the model. For binary traits, we use a modified R^2 called Nagelkerke’s R^2. 

To calculate the R^2 we must first try to predict the phenotype using the PRS as the predictor. 

```
prs_model <- glm(PHENOTYPE ~ SCORE, data = genetic_matrix_10, family=binomial)
```

We can know check how much variation in the phenotype our PRS model explains. 

```
result <- nagelkerke(prs_model)
result$Pseudo.R.squared.for.model.vs.null[3]
```

- [ ] **Lab Task 5: What percentage of the overall phenotypic variation is being explained by the PRS model? (Answer:72%))**

That's capturing almost all of the genetic variation in the trait! Unfortunatly, one should always be cautious in testing a model on data it was trained on i.e. testing a PRS on the same population that the GWAS was performed on. It's far better to test a model on an independant dataset. 

We will use a Japanese-ancestry population to test the predictive performance of our PRS. 

The GWAS was performed on individuals of European ancestry and so we can expect this PRS to perform worse for several reasons. 

Perhaps most importantly, remember that our SNPs may not be causal - the may simply be tagging the causal variant through LD structure. The LD structure differs between the European and Japanese population. Therefore, the SNP we have included in the PRS model may not be predictive at all in the Japanese cohort. 


Let's load in our Japanese genotype matrix. 
```
genetic_matrix_JPT <- read.table("genetic_matrix_JPT.raw", header=TRUE)
```

Now let's set up the PRS calcuation again (using our European GWAS effect sizes from before). 

```
JPT_prs <- data.frame(SCORE=numeric(nrow(genetic_matrix_JPT)), IID=character(nrow(genetic_matrix_JPT)), stringsAsFactors=FALSE)

for (i in 1:nrow(genetic_matrix_JPT)) {
snp_values <- genetic_matrix_JPT[i, 7:ncol(genetic_matrix_JPT)]
match_indices <- match(names(genetic_matrix_JPT)[7:ncol(genetic_matrix_JPT)], snps_to_include$SNP)
effect_sizes <- snps_to_include$EffectSize[match_indices]
JPT_prs$SCORE[i] <- sum(effect_sizes * snp_values, na.rm = TRUE)
JPT_prs$IID[i] <- genetic_matrix_JPT$IID[i]
}

head(JPT_prs)
```

Let us plot out the discrimination this time. 

```
genetic_matrix_JPT_with_prs <- merge(genetic_matrix_JPT, JPT_prs, by.x="IID", by.y ="IID")

ggplot(genetic_matrix_JPT_with_prs, aes(x = SCORE, fill = as.factor(PHENOTYPE))) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Distribution of PRS by Case (1) and Control (0) Status",
    x = "Polygenic Score (Unitless)",
    y = "Density",
    fill = "Phenotype"
  ) +
  theme_minimal()
```

<img src="https://github.com/att-y3-gcb/gwas_tutorial/blob/main/jptprs.PNG" alt="GWAS QC" width="50%">


It looks like our discrimnation is a lot worse, in the more distant population. This is because our SNPs are no longer properly tagging the underlying causal variation through the differences in LD between variants. 

- [ ] **Lab Task 6: What percentage of the overall phenotypic variation is being explained by the PRS model? (Answer:16%))**

Answer:
JPT_prs_model <- glm(PHENOTYPE ~ SCORE, data = genetic_matrix_JPT_with_prs, family=binomial)
result <- nagelkerke(JPT_prs_model)
result$Pseudo.R.squared.for.model.vs.null[3]


It is thus vitally important for the future of precision medicine that we ensure genomics studies (such as GWAS) are done on diverse populations, so that the benefits (PRS prediction) can reach everyone. 

<img src="https://github.com/att-y3-gcb/gwas_tutorial/blob/main/diversity.PNG" alt="GWAS QC" width="75%">
