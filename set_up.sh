#Start with HapMap_3_r3_1 files from 1_QC_GWAS.zip.

cd ~/GWA_tutorial/1_QC_GWAS/

shuf --random-source=<(yes 42) indepSNP.prune.in | head -n 50000 > subset_indep_snps.txt
sort -g -k 9,9 ../3_Association_GWAS/logistic_results.assoc_2.logistic | awk '{print $2}' | head -n 50 > subset_top_snps.txt
diff HapMap_3_r3_3.bim HapMap_3_r3_4.bim | grep 'rs' | awk '{print $3}' | shuf --random-source=<(yes 42) | head -n 37 > subset_geno_0.02.txt
diff HapMap_3_r3_8.bim HapMap_3_r3_7.bim | grep 'rs' | awk '{print $3}' | shuf --random-source=<(yes 42) | head -n 89 > subset_maf_0.05.txt
awk '{if ($1 == 23) print}' HapMap_3_r3_2.bim | shuf --random-source=<(yes 42) | head -n 20 > subset_x.txt
awk '{if ($1 == 25) print}' HapMap_3_r3_2.bim | shuf --random-source=<(yes 42) | head -n 2 > subset_y.txt

cat subset_indep_snps.txt subset_top_snps.txt subset_geno_0.02.txt subset_maf_0.05.txt subset_x.txt subset_y.txt > subset_snps.txt

#Remove individuals not suitable for this tutorial and subset to smaller number of SNPs (make sure to include significant SNPs!)
plink --extract subset_snps.txt --remove idv_remove.txt --bfile HapMap_3_r3_1 --make-bed --out genetic_matrix

mv genetic_matrix.log  genetic_matrix.old.log

#generates genotype matrix which can be loaded into R
plink --bfile genetic_matrix --recode A --out genetic_matrix

#Add in ebough missingness for an individuals to remove
sed -i 's/1458 NA12843 0 0 2 2 0 0 0 1 1 1 2 2 1 0 0 1 0 0 1 1 0 0/1458 NA12843 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA/g' genetic_matrix.raw

###############################################################################################
###########################Prepare 1000G Data####################################################
cp subset_indep_snps.txt ../2_Population_stratification
cp genetic_matrix.bim ../2_Population_stratification
cd ../2_Population_stratification
diff 1kG_MDS2.bim 1kG_MDS3.bim | grep 'rs' | awk '{print $3}' | shuf --random-source=<(yes 42) | head -n 70 > subset_geno_0.02.txt
diff 1kG_MDS4.bim 1kG_MDS5.bim | grep 'rs' | awk '{print $3}' | shuf --random-source=<(yes 42) | head -n 92 > subset_maf_0.05.txt

cat subset_indep_snps.txt subset_geno_0.02.txt subset_maf_0.05.txt > subset_snps.txt #note 1000G data is only chr1:22

#Remove individuals not suitable for this tutorial and subset to smaller number of SNPs 
plink --extract subset_snps.txt --bfile 1kG_MDS --make-bed --out bad_snps
plink --extract genetic_matrix.bim --bfile 1kG_MDS --make-bed --out good_snps
plink --bmerge good_snps --bfile bad_snps --make-bed --out 1000G_tutorial_data

mv 1000G_tutorial_data.log  1000G_tutorial_data.old.log

plink --bfile 1000G_tutorial_data --recode A --out 1000G_tutorial_data

#contents of idv_remove.txt (contains idvs failing due to heterozygotisity, relatedness and the non-founders)
'''
1330 NA12342
1459 NA12874
13291 NA07045
13281 NA12344
13291 NA06995
13291 NA06997
13292 NA07014
1330 NA12335
1330 NA12336
1334 NA10846
1334 NA10847
1340 NA07029
1341 NA06991
1344 NA10850
1345 NA07348
1345 NA07349
1346 NA10852
1347 NA10859
1349 NA10853
1350 NA10855
1350 NA10856
1353 NA12375
1353 NA12376
1354 NA12386
1358 NA12707
1358 NA12708
1362 NA10861
1375 NA10863
1377 NA10864
1377 NA10865
1408 NA10830
1408 NA10831
1416 NA10835
1418 NA10836
1418 NA10837
1420 NA10838
1420 NA10839
1421 NA10840
1423 NA10843
1424 NA10845
1444 NA12739
1444 NA12740
1447 NA12752
1447 NA12753
1451 NA12766
1451 NA12767
1454 NA12801
1454 NA12802
1456 NA12817
1456 NA12818
1458 NA12832
1459 NA12864
1459 NA12865
1463 NA12877
1463 NA12878
'''
