#Start with HapMap_3_r3_1 files from 1_QC_GWAS.zip.

#Remove individuals not suitable for this tutorial and subset to smaller number of SNPs (make sure to include significant SNPs!)
plink --extract subset_snps.txt --remove idv_remove.txt --bfile HapMap_3_r3_1 --make-bed --out genetic_matrix

mv genetic_matrix.log  genetic_matrix.old.log

#generates genotype matrix which can be loaded into R
plink --bfile genetic_matrix --recode A --out genetic_matrix

