#Start with HapMap_3_r3_1 files from 1_QC_GWAS.zip.
#cd ~/GWA_tutorial/1_QC_GWAS/

shuf --random-source=<(yes 42) indepSNP.prune.in | head -n 3000 > subset_indep_snps.txt
sort -g -k 9,9 ../3_Association_GWAS/logistic_results.assoc_2.logistic | awk '{print $2}' | head -n 50 > subset_top_snps.txt
diff HapMap_3_r3_3.bim HapMap_3_r3_4.bim | grep 'rs' | awk '{print $3}' | shuf --random-source=<(yes 42) | head -n 37 > subset_geno_0.02.txt
diff HapMap_3_r3_8.bim HapMap_3_r3_7.bim | grep 'rs' | awk '{print $3}' | shuf --random-source=<(yes 42) | head -n 89 > subset_maf_0.05.txt

cat subset_indep_snps.txt subset_top_snps.txt subset_geno_0.02.txt subset_maf_0.05.txt > subset_snps.txt

#Remove individuals not suitable for this tutorial and subset to smaller number of SNPs (make sure to include significant SNPs!)
plink --extract subset_snps.txt --remove idv_remove.txt --bfile HapMap_3_r3_1 --make-bed --out genetic_matrix

mv genetic_matrix.log  genetic_matrix.old.log

#generates genotype matrix which can be loaded into R
plink --bfile genetic_matrix --recode A --out genetic_matrix

