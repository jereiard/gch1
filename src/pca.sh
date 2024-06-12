#!/bin/bash
wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip
unzip ./plink_linux_x86_64_20231211.zip
wget https://s3.amazonaws.com/plink2-assets/plink2_linux_avx2_20240609.zip
unzip ./plink2_linux_amd_avx2_20240609.zip
rm ./plink_linux_x86_64_20231211.zip
rm ./plink2_linux_avx2_20240609.zip

parallel ::: \
"./plink2 --pfile chr14_AAC_release7 --make-bed --out AAC" \
"./plink2 --pfile chr14_AFR_release7 --make-bed --out AFR" \
"./plink2 --pfile chr14_AJ_release7 --make-bed --out AJ" \
"./plink2 --pfile chr14_AMR_release7 --make-bed --out AMR" \
"./plink2 --pfile chr14_CAH_release7 --make-bed --out CAH" \
"./plink2 --pfile chr14_CAS_release7 --make-bed --out CAS" \
"./plink2 --pfile chr14_EAS_release7 --make-bed --out EAS" \
"./plink2 --pfile chr14_EUR_release7 --make-bed --out EUR" \
"./plink2 --pfile chr14_FIN_release7 --make-bed --out FIN" \
"./plink2 --pfile chr14_MDE_release7 --make-bed --out MDE" \
"./plink2 --pfile chr14_SAS_release7 --make-bed --out SAS"

rm -rf chr14_???_release7.*

parallel ::: \
"./plink --bfile AAC --bmerge AFR --make-bed --out combined1" \
"./plink --bfile AJ --bmerge AMR --make-bed --out combined2" \
"./plink --bfile CAH --bmerge CAS --make-bed --out combined3" \
"./plink --bfile EAS --bmerge SAS --make-bed --out combined4" \
"./plink --bfile FIN --bmerge MDE --make-bed --out combined5" \
"./plink --bfile MDE --bmerge MDE --make-bed --out combined6"

parallel ::: \
"./plink --bfile combined1 --bmerge combined2 --make-bed --out combined11" \
"./plink --bfile combined3 --bmerge combined4 --make-bed --out combined12" \
"./plink --bfile combined5 --bmerge combined6 --make-bed --out combined13"

rm -rf combined?.*

parallel ::: \
"./plink --bfile combined11 --bmerge combined12 --make-bed --out combined21" \
"./plink --bfile combined13 --bmerge EUR --make-bed --out combined22"

rm -rf combined1?.*

./plink --bfile combined21 --bmerge combined22 --make-bed --out merged

rm -rf combined2?.*

plink2 --threads 92 --bfile merged --make-pgen --out chr14_all_ancestry_merged

# GRCh38, GCH1 flanking 10bp
./plink2 --threads 92 --pfile chr14_all_ancestry_merged --chr 14 --from-bp 54841998 --to-bp 54902836 --make-pgen --out gch1_all_ancestry_merged
./plink2 --threads 92 --pfile gch1_all_ancestry_merged --pca 10 --out gch1_pca_result

# Filter with MAF 0.01%, Hardy-Weinberg, Linkage disequilibrium
./plink2 --threads 92 --pfile gch1_all_ancestry_merged --maf 0.0001 --max-maf 0.05 --hwe 1e-6 --make-pgen --out gch1_all_ancestry_filtered
./plink2 --pfile gch1_all_ancestry_filtered --indep-pairwise 50 5 0.2 --out pruned_data
./plink2 --pfile gch1_all_ancestry_filtered --extract pruned_data.prune.in --make-pgen --out pruned_data_filtered
./plink2 --threads 92 --pfile pruned_data_filtered --pca 3 --out gch1_filtered_pca_result

# Filter with MAF 0.01%, Hardy-Weinberg, Linkage disequilibriumz
./plink2 --threads 92 --pfile gch1_all_ancestry_merged --max-maf 0.01 --hwe 1e-6 --make-pgen --out gch1_all_ancestry_filtered_re
./plink2 --threads 92 --pfile gch1_all_ancestry_filtered_re --indep-pairwise 50 5 0.2 --out pruned_data_re
./plink2 --threads 92 --pfile gch1_all_ancestry_filtered_re --extract pruned_data_re.prune.in --make-pgen --out pruned_data_filtered_re
./plink2 --threads 92 --pfile pruned_data_filtered_re --pca 3 --out gch1_filtered_pca_result_re

# Remaped Case-Control by Age of onset 50. Late-onset is case (remapping was done separately).
./plink2 --threads 92 --pfile gch1_all_ancestry_merged --keep samples_to_keep_onset.list--make-pgen --out gch1_all_ancestry_aos
./plink2 --threads 92 --pfile gch1_all_ancestry_aos --max-maf 0.01 --hwe 1e-6 --make-pgen --out gch1_all_ancestry_aos_filtered
./plink2 --threads 92 --pfile gch1_all_ancestry_aos_filtered --indep-pairwise 50 5 0.2 --out pruned_data_aos
./plink2 --threads 92 --pfile gch1_all_ancestry_aos_filtered --extract pruned_data_aos.prune.in --make-pgen --out pruned_data_filtered_aos
./plink2 --threads 92 --pfile pruned_data_filtered_aos --pca 3 --out gch1_filtered_pca_result_aos

#./plink2 -pfile pruned_data_filtered_aos --export vcf --out pruned_data_filtered_aos
#./plink2 --vcf pruned_data_filtered_aos.vcf --make-pgen --out pruned_data_filtered_aos_recon
#cp pruned_data_filtered_aos.psam pruned_data_filtered_aos_recon.psam
#./plink2 --threads 92 --pfile pruned_data_filtered_aos_recon --pca 3 --out gch1_filtered_pca_result_aos_recon
#./plink2 --threads 92 --pfile pruned_data_filtered_aos_recon --make-bed --out pruned_data_filtered_aos_recon
